from toolz.itertoolz import first, rest, peek
from operator import itemgetter as get
from toolz import compose, curry
from functools import partial
from toolz.dicttoolz import merge, dissoc
from itertools import ifilter, imap, groupby, takewhile, repeat, starmap, izip_longest
from contracts import contract
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

########################
# Function combintaors #
########################

''' just a way to compose filters and mapping functions together '''
def compose_transormer(transformer, combinator, *funcs):
    return partial(transformer, reduce(combinator, funcs))
_and = lambda f,g: lambda x: f(x) and g(x)
compose_filters = partial(compose_transormer, ifilter, _and)
compose_mappers = partial(compose_transormer, imap, compose)


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
    return reduce(_do_build, muts, ('', reference, 0))


##############
#  Mappers   #
##############

#when there are multiple alts, zip through with each of them
#zip(*alts), character by character. compare the percentages, and
#sum the percentages for each base. (groupby, sum) pick each character (call each base) based on the given rules (using call_base).
from toolz.dicttoolz import merge_with
#@contract(dp=int,ref=str,alts='dict(string,int)')
from toolz.dicttoolz import valfilter
def find_val(pred, d):
    res = valfilter(pred, d)
    return res if res else None
def find(pred, xs):
    res = ifilter(pred, xs)
    try:
        return next(res)
    except:
        return None
#TODO: Failing Case:
# a = {'ref': u'AT', 'pos': 2, 'AO': (51771, 41537, 42398, 9342), 'alt': [u'A',
# u'TT', u'AATTG', u'AAGAA'], 'chrom': u'o', 'DP': 87288}
#  bioframework.consensus.call_many(10, 80, a)


def call_base_multi_alts(min_depth, majority_percentage, dp, alts, ref):
    #TODO: majority_percentage gets ignored, so replace constants
    #TODO: behavior is undefined if sum(AO) > dp.
    import ipdb; ipdb.set_trace()
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


def call_many(min_depth, majority_percentage, _rec):
    #TODO: switch to generators
    if not hasattr(_rec['alt'], '__iter__'): #TODO: put this somewhere else
        rec = merge(_rec, dict(alt=[_rec['alt']], AO=[_rec['AO']]))
    else: rec = _rec
    muts = zip(rec['AO'], rec['alt'])
    ref, dp, pos = rec['ref'], rec['DP'], rec['pos']
    longest_len = max(map(lambda x: len(x[-1]), muts))
    longest_len = max(longest_len, len(ref))
    def fill_gap(r):
        ao, s = r
        return (ao, s + (longest_len - len(s)) * '-')
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
    #    return call_base(min_depth, majority_percentage, rec)


#@contract(min_depth=int, majoirty_percentage=int, rec=dict, returns='tuple(str, str, int)')
#DEPRECATED
def call_base(min_depth, majoirty_percentage, rec):
#DEPRECATED
    '''if a base is under the min depth, it gets called as an `N`.
    If the alternate base is >= `majoirty_percentage`, it will be
    converted to an ambiguous base if it includes multiple bases.
    Otherwise, the reference base replaces itself (gets ignored).
    Returns a tuple of the form:
        (referene_base, alternate_base, position)
        where reference_base == alternate_base if the reference was preferred.'''
    get_ambiguous = compose(lambda x: AMBIGUITY_TABLE.get(x, x), ''.join, sorted)
    if rec['DP'] < min_depth: alt = 'N'# (rec['ref'], 'N', rec['POS'])
    elif alt_over_percent(majoirty_percentage)(rec):
        alt = get_ambiguous(rec['alt'])
    else: alt = rec['ref']
    return (rec['ref'], alt, rec['pos'])

def flatten_vcf_record(rec):
    return merge({
  'alt' : rec.ALT, 'ref' : rec.REF,
  'pos' : rec.POS, 'chrom' : rec.CHROM},
        rec.INFO)

##############
#  Filters   #
##############
ref_and_alt_differ = lambda x: x[0] != x[1]
alt_over_percent =  lambda n: lambda rec: rec['AO']/float(rec['DP']) > (n/float(100))


##############
# Group By   #
##############
def group_muts_by_refs(references, muts):
    '''group and sort the mutations so that they match the order of the references.'''
    mut_groups = groupby(muts, get('CHROM'))
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
    mappers = [partial(call_many, mind, majority)]
    filters = [ref_and_alt_differ]
    muts_by_ref = group_muts_by_refs(references, muts)
    muts = imap(get(0), muts_by_ref)
    def single_ref_consensus(recs, ref):
        fix_recs = compose(compose_filters(filters), compose_mappers(mappers))
        ref_str = str(ref.seq)
        return make_consensus(ref_str, fix_recs(recs))
    return references, imap(single_ref_consensus, muts, references)



##########
# I/O    #
##########

def consensus_str(ref, consensus):
    return ">{0}:Consensus\n{1}".format(ref.id, consensus)

make_fasta_str = compose('\n'.join,
                         partial(imap, consensus_str), all_consensuses)

#@contract(ref_fasta=str, vcf=str, mind=int, majority=int)
def run(ref_fasta, freebayes_vcf, outpath, mind, majority):
    refs = SeqIO.parse(ref_fasta, 'fasta')
    with open(freebayes_vcf, 'r') as vcf_handle, open(outpath, 'w') as out:
        muts = imap(flatten_vcf_record, vcf.Reader(vcf_handle))
        result = make_fasta_str(refs, muts, mind, majority)
        out.write(result)
    return 0

