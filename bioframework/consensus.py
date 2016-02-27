"""
Usage:
     consensus --ref <ref> --vcf <vcf> [--mind <mind>] [--majority <majority>] [-o <output>]

Options:
    --ref=<ref>             Reference fasta file
    --vcf=<vcf>             VCF output
    --majority=<majority>   Percentage required [default: 80]
    --mind=<mind>           minimum depth to call base non-N [default: 10]
    --output,-o=<output>       output file [default: ]
"""
from toolz.itertoolz import first, rest, peek
from operator import itemgetter as get
from toolz import compose, curry
from functools import partial
from toolz.dicttoolz import merge, dissoc, merge_with, valfilter
from itertools import ifilter, imap, groupby, takewhile, repeat, starmap, izip_longest
#from contracts import contract
import os, sys
from docopt import docopt
from schema import Schema, Use
from Bio import SeqIO, SeqRecord
import vcf

#############
# Constants #
#############

AMBIGUITY_TABLE = { 'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'N': 'N',
                       'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT':
                       'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D',
                       'CGT': 'B', 'ACGT': 'N' }

MAJORITY_PERCENTAGE = 80
MIN_DEPTH = 10

###########
# Reducer #
###########
#@contract(reference=str, muts='seq(tuple(str, str, int))'   )
def make_consensus(reference, muts):
    ''' Actually builds a consensus string by recursively applying
          the mutations.'''
    def _do_build((accString, string, lastPos), (x, y, bigPos)):
        pos = bigPos - lastPos
        return (accString + (string[:pos] + y), string[pos+len(x):],  bigPos+len(x))
    return reduce(_do_build, muts, ('', reference, 0))[0]


##############
#  Mappers   #
##############

#TODO: Failing Case:
# a = {'ref': u'AT', 'pos': 2, 'AO': (51771, 41537, 42398, 9342), 'alt': [u'A',
# u'TT', u'AATTG', u'AAGAA'], 'chrom': u'o', 'DP': 87288}
#  bioframework.consensus.call_many(10, 80, a)

#@contract(dp=int,ref=str,alts='dict(string,int)')
def call_base_multi_alts(min_depth, majority_percentage, dp, alts, ref):
    """when there are multiple alts, zip through with each of them
    zip(*alts), character by character. compare the percentages, and
    sum the percentages for each base. (groupby, sum) pick each character
    (call each base) based on the given rules (using call_base)."""
    #TODO: majority_percentage gets ignored, so replace constants
    #TODO: behavior is undefined if sum(AO) > dp.
    #import ipdb; ipdb.set_trace()
    if dp < min_depth: #could call REF here sometimes
        return 'N'
    total_ao = lambda: sum(alts.values()) # avoid evaluating unless necessary

    if ref is None: # this is an insert
        if total_ao()/float(dp) < .50:
            return ref
        # if the insert is above threshold, keep going and call the insert like a normal base
    if '-' in alts:
        if alts['-']/float(dp) > .50: # a deletion
            return ''
        dp -= alts['-']
        alts_without_insert = dissoc(alts, '-')
    else:
        alts_without_insert = alts
    over_depth = lambda x: lambda depth: depth/float(dp) > x
    picked_alt = valfilter(over_depth(0.8), alts_without_insert)
    if picked_alt:
        return picked_alt.keys()[0]
    #add ref so that it will be considered in creating ambiguous base
    alts_with_ref = merge(alts_without_insert, ({ref : (dp - total_ao()) } if ref else {}))
    over20 = valfilter(over_depth(0.2), alts_with_ref)
    as_ambiguous = ''.join(sorted(over20.keys()))
    # this could return a single base, (including the reference), becuase i.e.  A => A in the ambiguity table
    return AMBIGUITY_TABLE[as_ambiguous]


def call_many(min_depth, majority_percentage, rec):
    #TODO: switch to generators
    muts = zip(rec['AO'], rec['alt'])
    ref, dp, pos = rec['ref'], rec['DP'], rec['pos']
    longest_len = max(map(lambda x: len(x[-1]), muts))
    longest_len = max(longest_len, len(ref))
    def fill_gap(r):
        ao, s = r
        return (ao, str(s) + (longest_len - len(s)) * '-')
    xs = map(fill_gap, muts) # fill in the shorter alts with '-'.
    def merge_sum(x,y):
        return x if y is None else (y if x is None else merge_with(sum, x, y))
    def seq_count(acc, ao_and_nts):
        ao, nts = ao_and_nts
        return map(merge_sum, acc, [{nt:ao} for nt in nts])
    # create a list of {base : count}, where the index matches the position
    mut_dicts = reduce(seq_count, xs, [{}])
    base_caller = partial(call_base_multi_alts, min_depth, majority_percentage, dp)
    res = map(base_caller, mut_dicts, ref)
    # trim None values at the end, (which indicate deletion)
    result = takewhile(bool, res)
    return (ref, ''.join(result), pos)

def flatten_vcf_record(rec):
    _rec = merge({
  'alt' : rec.ALT, 'ref' : rec.REF,
  'pos' : rec.POS, 'chrom' : rec.CHROM},
        rec.INFO)
    if not hasattr(_rec['alt'], '__iter__'): #TODO: put this somewhere else
        return merge(_rec, dict(alt=[_rec['alt']], AO=[_rec['AO']]))
    else: return _rec

##############
# Group By   #
##############
def group_muts_by_refs(references, muts):
    '''group and sort the mutations so that they match the order of the references.'''
    mut_groups = groupby(muts, get('chrom'))
    def index_of_ref(key):
        return sum(1 for _ in takewhile(lambda x: x.id != key, references))
    muts_by_ref = sorted(mut_groups, key=index_of_ref)
    return muts_by_ref



###############
# Runner      #
###############

#@contract(references='SeqRecord', muts='seq(dict)', mind=int, majority=int)
def all_consensuses(references, muts, mind, majority):
    ''' generates conesnsuses, including for flu and other mult-reference VCFs.
    applies filters and base callers to the mutations.
    then builds the consensus using these calls and `make_consensus`'''
    muts_by_ref = group_muts_by_refs(references, muts)
    mut_groups = map(compose(list, get(1)), muts_by_ref)
    def single_consensus(muts, ref):
        the_muts = map(partial(call_many, mind, majority), muts)
        ref_and_alt_differ = lambda x: x[0] != x[1]
        real_muts = filter(ref_and_alt_differ, the_muts)
        return make_consensus(str(ref.seq), real_muts)
    return references, imap(single_consensus, mut_groups, references)


##########
# I/O    #
##########
def consensus_str(ref, consensus):
    return ">{0}:Consensus\n{1}".format(ref.id, consensus)


#@contract(ref_fasta=str, vcf=str, mind=int, majority=int)
def run(ref_fasta, freebayes_vcf, outfile, mind, majority):
    refs = SeqIO.parse(ref_fasta, 'fasta')
    with open(freebayes_vcf, 'r') as vcf_handle:
        muts = imap(flatten_vcf_record, vcf.Reader(vcf_handle))
        refs, muts = list(refs), list(muts)
        refs, seqs = all_consensuses(refs, muts, mind, majority)
        strings = imap(consensus_str, refs, seqs)
        result = '\n'.join(strings)
        outfile.write(result)
        outfile.close()
    return 0

def main():
    scheme = Schema(
        { '--vcf' : os.path.isfile,
          '--ref' : os.path.isfile,
          '--majority' : Use(int),
          '--mind' : Use(int),
          '--output' : Use(lambda x: sys.stdout if not x else open(x, 'w'))})
    raw_args = docopt(__doc__, version='Version 1.0')
    args = scheme.validate(raw_args)
    run(args['--ref'], args['--vcf'], args['--output'],
        args['--mind'], args['--output'])

if __name__ == '__main__':
    main()
