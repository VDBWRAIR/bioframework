from toolz.itertoolz import first, rest, peek
from operator import itemgetter as get
from toolz import compose, curry
from functools import partial
from toolz.dicttoolz import merge
from itertools import ifilter, imap, groupby, takewhile
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
@contract(reference=str, muts='seq(tuple(str, str, int))'   )
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
@contract(min_depth=int, majoirty_percentage=int, rec=dict, returns='tuple(str, str, int)')
def call_base(min_depth, majoirty_percentage, rec):
    '''if a base is under the min depth, it gets called as an `N`.
    If the alternate base is >= `majoirty_percentage`, it will be
    converted to an ambiguous base if it includes multiple bases.
    Otherwise, the reference base replaces itself (gets ignored).
    Returns a tuple of the form:
        (referene_base, alternate_base, position)
        where reference_base == alternate_base if the reference was preferred.'''
    get_ambiguous = compose(AMBIGUITY_TABLE.__getitem__, ''.join, sorted)
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
alt_over_percent =  lambda n: lambda rec: (n/float(100)) > rec['AO']/float(rec['DP'])


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
@contract(references='SeqRecord', muts='seq(dict)', mind=int, majority=int)
def all_consensuses(references, muts, mind, majority):
    ''' generates conesnsuses, including for flu and other mult-reference VCFs.
    applies filters and base callers to the mutations.
    then builds the consensus using these calls and `make_consensus`'''
    mappers = [partial(call_base, mind, majority)]
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

@contract(ref_fasta=str, vcf=str, mind=int, majority=int)
def run(ref_fasta, freebayes_vcf, outpath, mind, majority):
    refs = SeqIO.parse(ref_fasta, 'fasta')
    with open(freebayes_vcf, 'r') as vcf_handle, open(outpath, 'w') as out:
        muts = imap(flatten_vcf_record, vcf.Reader(vcf_handle))
        result = make_fasta_str(refs, muts, mind, majority)
        out.write(result)
    return 0

