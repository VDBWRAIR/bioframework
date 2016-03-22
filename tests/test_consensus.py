from __future__ import division
from biotest.biohypothesis import ref_with_vcf_dicts_strategy_factory, \
    make_seqrec, vcf_dict_strategy_factory
from hypothesis import given
#from fn import _
from hypothesis import strategies as st
from hypothesis import given, assume
from operator import itemgetter as get
from bioframework.consensus import call_many, all_consensuses, make_consensus, VCFRow
import string
import itertools
import unittest

simple_vcf_dict_strategy = st.tuples(st.text(string.ascii_letters),
                     st.integers(min_value=1),
                     st.text(alphabet='ACTGN', min_size=1, max_size=6)) \
                     .flatmap(lambda tup:\
                          vcf_dict_strategy_factory(*tup)).map(lambda d: VCFRow(**d))
pos_int = st.integers(min_value=0)

#TODO: these 10, 80 for trhesh and majority_percentage should be factored out and possibly be strategies themselves
def just_ref(*args):
    return next(all_consensuses(*args)[1])[0]
class CallBaseHypothesisTest(unittest.TestCase):
    @given(simple_vcf_dict_strategy, pos_int)
    def test_under_mind_is_N(self, mut, mind):
        assume(mut.DP < mind)
        result = call_many(mind, 80, mut)[1]
        self.assertTrue(all(map(lambda x: x == 'N', result)))

    @given(simple_vcf_dict_strategy)
    def test_ao_under_minority_is_ref(self, mut):
        assume(sum(mut.AO) / mut.DP < 0.2)
        result = call_many(0, 80, mut)[1]
        self.assertEquals(result, mut.ref)

    @given(simple_vcf_dict_strategy)
    def test_over_majority_is_alt(self, mut):
        #TODO: this is slow
        assume(sum(mut.AO) / mut.DP > 0.8)
        assume(len(mut.alt) == 1)
        result = call_many(0, 80, mut)[1]
        self.assertEquals(result, mut.alt[0])

#Commented out because it's not actually always true,
# e.g. mut={'ref': u'AA', 'pos': 1, 'AO': [784313725491], 'alt': [u'A'],
# 'chrom': u'', 'DP': 3921568627454})
# should result in AA
#    @given(simple_vcf_dict_strategy)
#    def test_over_minoriy_is_not_ref(self, mut):
#        assume(sum(mut.AO) / mut.DP > 0.2)
#        result = call_many(0, 80, mut)[1]
#        self.assertNotEquals(result, mut.ref)

class ConsesusExampleTest(unittest.TestCase):
    def test_make_consensus_example(self):
        muts = [('CG', 'TT', 1),
                ('C', 'A', 5)]
        ref = 'ACGTACGT'
        expected = 'ATTTAAGT'
        actual = make_consensus(ref, muts)[0]
        self.assertEquals(expected, actual)

    def test_single_example(self):
        raw_muts = [{
            'pos' : 2,
            'ref' : 'CG',
            'alt' : ['TT'],
            'AO'  : [99],
            'DP'  : 100,
            'chrom' : 'X'
        },
        {
            'pos' : 6,
            'ref' : 'C',
            'alt' : ['G', 'A'],
            'AO' : [15, 120],
            'DP' : 150,
            'chrom' : 'X'
        }]
        muts = map(lambda d: VCFRow(**d), raw_muts)
        ref = make_seqrec('X', 'ACGTACGT')
        expected = 'ATTTAAGT'
        result = just_ref([ref], muts, 10, 80)
        self.assertEquals(expected, result)
ref_with_vcf_dicts_strategy = ref_with_vcf_dicts_strategy_factory().map(
    lambda (r, muts): (make_seqrec(muts[0]['chrom'], r), map(lambda d: VCFRow(**d), muts)))
from collections import Counter
countof = lambda c: lambda x: Counter(x).get(c, 0)
def run_cons(*args):
    _, alt_and_cons = all_consensuses(*args)
    cons, alts  = zip(*alt_and_cons)
    return  cons[0], alts[0]
class ConsensusHypothesisTest(unittest.TestCase):
    #ref_and_muts=(SeqRecord(seq=Seq(u'AAAAAAAAAA', IUPACAmbiguousDNA()), id=u'', name='<unknown name>', description='', dbxrefs=[]), [
#    {'ref': u'A', 'pos': 1, 'AO': [479, 777, 119, 604], 'alt': [u'G', u'C', u'G', u'TG'], 'chrom': u'', 'DP': 2635},
#    {'ref': u'A', 'pos': 3, 'AO': [291, 241, 583, 420], 'alt': [u'CTG', u'C', u'G', u'C'], 'chrom': u'', 'DP': 1627}]), rand=random.seed(0))
#  AssertionError: 1 != 0

    @given(ref_with_vcf_dicts_strategy, st.random_module())
    def test_n_count(self, ref_and_muts, rand):
        ref, muts = ref_and_muts
        originalNs = countof('N')(ref)
        alts = map(lambda x: x.alt, muts)
        refs = map(lambda x: x.ref, muts)
        assume(not filter(lambda x: 'N' in x, itertools.chain(*alts)))
        assume(not filter(lambda x: len(x) > 1, itertools.chain(*alts)))
        assume(not filter(lambda x: len(x) > 1, refs))
        # needed because  ACGT -> N
        expectedNs = len(filter(lambda x: x.DP < 10, muts))  + originalNs
        result = just_ref([ref], muts, 10, 80)
        self.assertEquals(countof('N')(result), expectedNs)

    @given(ref_with_vcf_dicts_strategy)
    def test_less_or_equal_length_when_no_inserts(self, ref_and_muts):
        ref, muts = ref_and_muts
        cons, alts = run_cons([ref], muts, 10, 80)
        assume(not any(map(lambda x: len(x[0]) < len(x[1]), alts)))
        self.assertGreaterEqual(len(ref), len(cons))

    @given(ref_with_vcf_dicts_strategy)
    def assume_greater_or_equal_length_when_no_deletions(self, ref_and_muts):
        ref, muts = ref_and_muts
        def has_deletion(mut):
             filter(lambda x: len(x) < mut.ref, mut.alt)
        assume(not any(map(has_deletion, muts)))
        result = just_ref([ref], muts, 10, 80)
        self.assertLesserEqual(len(ref), len(result))

class ConsensusMetamorphicTests(unittest.TestCase):
    '''paper on metamorphic testing: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-24 '''
    @given(ref_with_vcf_dicts_strategy, st.integers(), st.integers())
    def test_more_or_equal_ns_with_lower_threshold(self, ref_and_muts, n1, n2):
        ref, muts = ref_and_muts
        assume(n1 < n2)
        cons1 = just_ref([ref], muts, n1, 80)
        cons2 = just_ref([ref], muts, n2, 80)
        nsCount1, nsCount2 = countof('N')(cons1), countof('N')(cons2)
        self.assertLessEqual(nsCount1, nsCount2)

    @given(ref_with_vcf_dicts_strategy)
    def test_consensus_from_consensus_contains_more_alts(self, ref_and_muts):
        ref, muts = ref_and_muts
        assume(not any(map(lambda x: len(x.alt) > 1, muts)))
        n1 = 10
        cons1, alts = run_cons([ref], muts, n1, 80)
        assume(not any(map(lambda x: len(x[0]) > len(x[1]), alts)))
        cons2, _ = run_cons([make_seqrec(muts[0].chrom, cons1)], muts, n1, 80)
        picked_alts = map(get(1), alts)
        altCounts1 = sum(map(lambda f: f(cons1),  map(countof, picked_alts)))
        altCounts2 = sum(map(lambda f: f(cons2),  map(countof, picked_alts)))
        self.assertLessEqual(altCounts1, altCounts2)


        #NOTE: the below test appears to be meaningless,
        # ass the values are always equal
    percent = st.integers(min_value=0, max_value=100)
    @given(ref_with_vcf_dicts_strategy, percent, percent)
    def test_lower_majority_required_contains_more_alts(self, ref_and_muts, p1, p2):
        ref, muts = ref_and_muts
        assume(p1 < p2)
        assume(not any(map(lambda x: len(x.alt) > 1, muts)))
        n1 = 10
        cons1, alts = run_cons([ref], muts, n1, p1)
        assume(not any(map(lambda x: len(x[0]) > len(x[1]), alts)))
        cons2, _ = run_cons([ref], muts, n1, p2)
        picked_alts = map(get(1), alts)
        altCounts1 = sum(map(lambda f: f(cons1),  map(countof, picked_alts)))
        altCounts2 = sum(map(lambda f: f(cons2),  map(countof, picked_alts)))
        self.assertLessEqual(altCounts1, altCounts2)

        #TODO: could test for exceptions when e.g.,
        # 1. POS is repeated
        # 2. POS is greater than the size of sequence
