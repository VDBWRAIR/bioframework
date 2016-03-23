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
#stdlib
#TO type-check: (requires python3)
# $ MYPYPATH=$HOME/bioframework/mypy mypy --py2 --disallow-untyped-calls  bioframework/consensus.py

from operator import itemgetter as get
from functools import partial
from itertools import ifilter, ifilterfalse, imap, groupby, takewhile, repeat, starmap, izip_longest
import os, sys
import collections

from typing import Tuple, Dict, List, Iterator, Iterable, Any, Callable, NamedTuple, BinaryIO

from Bio import SeqIO #done
from Bio.SeqRecord import SeqRecord #done
import vcf #done
from vcf.model import _Record
import sh #todo
#from toolz import compose
from toolz.dicttoolz import merge, dissoc, merge_with, valfilter, keyfilter #done
#from docopt import docopt #ignore
#from schema import Schema, Use #ignore
#from contracts import contract, new_contract #can ignore
#from mypy.types import VCFRow
#############
# Constants #
#############
VCFRow = NamedTuple("VCFRow",
                    [('ref', str),
                     ('AO', List[int]),
                     ('DP', int),
                     ('QA', List[int]),
                     ('QR', int),
                     ('chrom',str),
                     ('pos', int),
                     ('alt', List[str])])
AMBIGUITY_TABLE = { 'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'N': 'N', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N' }

MAJORITY_PERCENTAGE = 80
MIN_DEPTH = 10
PMut = NamedTuple('PMut', [('ALT', str), ('AO', int), ('DP', int), ('QA', int), ('QR', int)])
Mut = Tuple[str, str, int]
###########
# Reducer #
###########
#@contract(reference='string', muts='list(tuple(string, string, int))'   )
def make_consensus(reference, muts):
    # type: (str, List[Mut]) -> Tuple[str, List[Mut]]
    ''' Actually builds a consensus string by recursively applying
          the mutations.'''
    def _do_build(t1, t2): # type: (Tuple[str,str,int], Tuple[str,str,int]) -> Tuple[str,str,int]
        (accString, string, lastPos), (x, y, bigPos) = t1, t2
        pos = bigPos - lastPos
        return (accString + (string[:pos] + y), string[pos+len(x):],  bigPos+len(x))
    result, remaining, _ = reduce(_do_build, muts, ('', reference, 0))
    return result + remaining, muts


##############
#  Mappers   #
##############

#TODO: Failing Case:
# a = {'ref': u'AT', 'pos': 2, 'AO': (51771, 41537, 42398, 9342), 'alt': [u'A',
# u'TT', u'AATTG', u'AAGAA'], 'chrom': u'o', 'DP': 87288}
#  bioframework.consensus.call_many(10, 80, a)

#@contract(min_depth='number,>=0', majority_percentage='number,>=0,<=100',dp='number,>=0', ref='string|None',alts='dict(string: number)')
def call_base_multi_alts(min_depth, majority_percentage, dp, alts, ref):
    # type: (int, int, int, List[PMut], str) -> str
    """when there are multiple alts, zip through with each of them
    zip(*alts), character by character. compare the percentages, and
    sum the percentages for each base. (groupby, sum) pick each character
    (call each base) based on the given rules (using call_base)."""
    #TODO: majority_percentage gets ignored, so replace constants
    #TODO: behavior is undefined if sum(AO) > dp.
#    if dp < min_depth: #could call REF here sometimes
#        return 'N'

#    """
#    1) "if the total quality of the alternates is < 200 (25 * 8), don't call it an alternate."
#    2) "if the total quality of the references is < 200, don't call the reference."
#    3) "if neither the reference nor the alternate can be called, call an N."
#    """
#    altsTooLow = sum(map(lambda x: x.QA)) < 200
#    refTooLow = alts[0].QR < 200
#
#    if altsTooLow and refTooLow:
#        return 'N'
#
#    elif altsTooLow:
#        return ref

    get_ao = lambda x: x.AO
    total_ao = sum(map(get_ao, alts))

    if ref is None: # this is an insert
        if total_ao/float(dp) < .50:
            return ref
        # if the insert is above threshold, keep going and call the insert like a normal base
    is_insert = lambda x: x.ALT == '-' # type: Callable[[PMut],bool]
    inserts = list(filter(is_insert, alts))
    if inserts:
        assert len(inserts) == 1
        ins = inserts[0]
        if ins.AO/float(dp) > .50: # a deletion
            return ''
        dp -= ins.AO
        alts_without_insert = list(ifilterfalse(is_insert, alts))
    else:
        alts_without_insert = alts
    over_depth = lambda x: lambda y: y.AO/float(dp) > x
    picked_alt = filter(over_depth(0.8), alts_without_insert)
    if picked_alt:
        return picked_alt[0].ALT
    #add ref so that it will be considered in creating ambiguous base
    if ref:
        alts.append(PMut(ALT=ref, AO=(dp-total_ao), DP=dp, QA=-1, QR=-1))
    over20 = filter(over_depth(0.2), alts)
    as_ambiguous = ''.join(sorted(map(lambda x: x.ALT, over20)))
    # this could return a single base, (including the reference), becuase i.e.  A => A in the ambiguity table
    return AMBIGUITY_TABLE[as_ambiguous] if as_ambiguous != '' else ''

def pmuts_from_row(row):
    # type: (VCFRow) -> List[PMut]
    muts = zip(row.QA, row.AO, row.alt)
    longest_len = max(map(lambda x: len(x[-1]), muts))
    longest_len = max(longest_len, len(row.ref))
    def fill_gap(r):
        qa, ao, s = r
        return (qa, ao, str(s) + (longest_len - len(s)) * '-')
    def merge_sum(x,y):
        def combine(tups):
            return sum(map(get(0), tups)), sum(map(get(1), tups))
        #return x if y is None else (y if x is None else merge_with(sum, x, y))
        return x if y is None else (y if x is None else merge_with(combine, x, y))
    def seq_count(acc, ao_and_nts):
        qa, ao, nts = ao_and_nts
        #return map(merge_sum, acc, [{nt:ao, 'qa': qa} for nt in nts])
        return map(merge_sum, acc, [{nt:(ao, qa)} for nt in nts])
    xs = map(fill_gap, muts) # fill in the shorter alts with '-'.
    # create a list of {base : count}, where the index matches the position
    mut_dicts = reduce(seq_count, xs, [{}]) # type: Iterable[Dict[str,int]]
    def to_pmuts(d): # type: (Dict[str,int]) -> List[PMut]
        #qa = d.pop('qa')
        #assert len(d.values()) == 1
        make = lambda k,v: PMut(k, AO=v[0], QA=v[1], DP=row.DP,QR=row.QR)
        return list(starmap(make, d.items()))
#        res = PMut(d.keys()[0], d.values()[0], row.DP, qa, row.QR)
#        return res
    res = list(map(to_pmuts, mut_dicts))
    return res

def call_many(min_depth, majority_percentage, rec):
    # type: (int, int, VCFRow) -> Mut
    muts_by_column = pmuts_from_row(rec)
    #TODO: switch to generators
    base_caller = lambda m,r: call_base_multi_alts(
        min_depth, majority_percentage, rec.DP, m, r) #   # # ?Callable[[Dict[Any,Any], str], str]
    res = map(base_caller, muts_by_column, rec.ref)
    # trim None values at the end, (which indicate deletion)
    result = takewhile(bool, res)
    return (rec.ref, ''.join(result), rec.pos)

#@contract(rec='dict',returns='dict')
def flatten_vcf_record(rec):
    # type: (_Record) -> VCFRow
    fields = ['ref', 'AO', 'DP', 'QA', 'QR', 'chrom' 'pos', 'alt']
    _rec = merge({
  'alt' : rec.ALT, 'ref' : rec.REF,
  'pos' : rec.POS, 'chrom' : rec.CHROM},
        rec.INFO)
    if not hasattr(_rec['alt'], '__iter__'): #TODO: put this somewhere else
        d = merge(_rec, dict(alt=[_rec['alt']],
                             AO=[_rec['AO']],
                             QA=[_rec['QA']]))
    else: d = _rec
    d = keyfilter(fields.__contains__, d)
    return VCFRow(**d)

##############
# Group By   #
##############
#NOTE: could possibly drop lists, use fn.Stream all the time,
# and write a Stream instance for contracts like:
# https://github.com/AndreaCensi/contracts/blob/831ec7a5260ceb8960540ba0cb6cc26370cf2d82/src/contracts/library/lists.py
#@contract(references='list[N]($SeqRecord),N>0', muts='list(dict)',returns='tuple(list(dict))')
def group_muts_by_refs(references, muts):
    # type: (List[SeqRecord], List[VCFRow]) -> List[List[VCFRow]]
    '''group and sort the mutations so that they match the order of the references.'''
    #NOTE: muts will already be "sorted" in that they are grouped together in the vcf
    #fix the groupby so it doesn't incidentally drain the first object of the group
    unzip = lambda x: zip(*x)
    chroms, groups = unzip(map(lambda kv: (kv[0], list(kv[1])), groupby(muts, lambda x: x.chrom)))
    #@contract(key='tuple(string,list)')
    def index_of_ref(key): # type: (Tuple[str, List[SeqRecord]]) -> int
        chrom=key[0]
        index_of_chrom =  map(lambda x: x.id, references).index(chrom)
        return index_of_chrom
    _, muts_by_ref = unzip(sorted(zip(chroms, groups), key=index_of_ref))
    return muts_by_ref



###############
# Runner      #
###############

#@contract(references='SeqRecord', muts='seq(dict)', mind=int, majority=int)
def all_consensuses(references, muts, mind, majority):
    # type: (List[SeqRecord], List[VCFRow], int, int) -> Tuple[List[SeqRecord], Iterable[Tuple[str, List[Mut]]]]
    ''' generates conesnsuses, including for flu and other mult-reference VCFs.
    applies filters and base callers to the mutations.
    then builds the consensus using these calls and `make_consensus`'''
    muts_by_ref = group_muts_by_refs(references, muts)
    def single_consensus(muts, ref):
        # type: (List[VCFRow], SeqRecord) -> Tuple[str, List[Mut]]
        #the_muts = map(partial(call_many, mind, majority), muts)
        the_muts = map(lambda x: call_many(mind, majority, x), muts)
        ref_and_alt_differ = lambda x: x[0] != x[1]
        # vcf is index-starting-at-1
        #real_muts = map(lambda (a,b,pos): (a,b,pos-1), filter(ref_and_alt_differ, the_muts))
        real_muts = map(lambda x: (x[0], x[1], x[2] - 1), filter(ref_and_alt_differ, the_muts))
        return make_consensus(str(ref.seq), real_muts)
    return references, imap(single_consensus, muts_by_ref, references)


##########
# I/O    #
##########
def consensus_str(ref, consensus): # type: (SeqRecord, str) -> str
    return ">{0}:Consensus\n{1}".format(ref.id, consensus)

def zero_coverage_positions(bam_file, ref_file): # type: (str, str) -> Iterable[int]
    pileup = sh.Command('mpileup')(bam_file, f=ref_file, _iter=True)
    get_pos = lambda x: int(x.split()[1]) # type: Callable[[str],int]
    return imap(get_pos, pileup)

#TODO: is pileup 0-based or 1-based index?
def trim_ref(ref, positions): # type: (str, Iterator[int]) -> str
    start, end = next(positions), collections.deque(positions, 1)[0]
    return '-'*start + ref[:start:end] + '-'*(len(ref) - end)



#@contract(ref_fasta=str, vcf=str, mind=int, majority=int)
def run(ref_fasta, freebayes_vcf, outfile, mind, majority):
    # type: (str, str, BinaryIO, int, int) -> int
    _refs = SeqIO.parse(ref_fasta, 'fasta')
    with open(freebayes_vcf, 'r') as vcf_handle:
        _muts = map(flatten_vcf_record, vcf.Reader(vcf_handle))
        refs, muts = list(_refs), list(_muts)
        the_refs, seqs_and_muts = all_consensuses(refs, muts, mind, majority)
        strings = imap(consensus_str, the_refs, imap(get(0), seqs_and_muts))
        result = '\n'.join(strings)
        outfile.write(result)
        outfile.close()
    return 0





#def main(): # type: () -> None
#    scheme = Schema(
#        { '--vcf' : os.path.isfile,
#          '--ref' : os.path.isfile,
#          '--majority' : Use(int),
#          '--mind' : Use(int),
#          '--output' : Use(lambda x: sys.stdout if not x else open(x, 'w'))})
#    raw_args = docopt(__doc__, version='Version 1.0')
#    args = scheme.validate(raw_args)
#    run(args['--ref'], args['--vcf'], args['--output'],
#        args['--mind'], args['--output'])

#if __name__ == '__main__':
#    main()
