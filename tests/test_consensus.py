from __future__ import division
from biotest.biohypothesis import ref_with_vcf_dicts_strategy_factory, \
    make_seqrec, vcf_dict_strategy_factory
from hypothesis import given
from fn import _
from hypothesis import strategies as st
from hypothesis import given, assume
from bioframework.consensus import call_many, all_consensuses, make_consensus
import string
import unittest

simple_vcf_dict_strategy = st.tuples(st.text(string.ascii_letters),
                     st.integers(min_value=1),
                     st.text(alphabet='ACTGN', min_size=1, max_size=6)) \
                     .flatmap(lambda tup:\
                          vcf_dict_strategy_factory(*tup))
pos_int = st.integers(min_value=0)
def just_ref(*args):
    return next(all_consensuses(*args)[1])
class CallBaseHypothesisTest(unittest.TestCase):
    @given(simple_vcf_dict_strategy, pos_int)
    def test_under_mind_is_N(self, mut, mind):
        assume(mut['DP'] < mind)
        result = call_many(mind, 80, mut)[1]
        self.assertTrue(all(map(lambda x: x == 'N', result)))

    @given(simple_vcf_dict_strategy)
    def test_ao_under_minority_is_ref(self, mut):
        assume(sum(mut['AO']) / mut['DP'] < 0.2)
        result = call_many(0, 80, mut)[1]
        self.assertEquals(result, mut['ref'])

    @given(simple_vcf_dict_strategy)
    def test_over_majority_is_alt(self, mut):
        #TODO: this is slow
        assume(sum(mut['AO']) / mut['DP'] > 0.8)
        assume(len(mut['alt']) == 1)
        result = call_many(0, 80, mut)[1]
        self.assertEquals(result, mut['alt'][0])

#Commented out because it's not actually always true,
# e.g. mut={'ref': u'AA', 'pos': 1, 'AO': [784313725491], 'alt': [u'A'],
# 'chrom': u'', 'DP': 3921568627454})
# should result in AA
#    @given(simple_vcf_dict_strategy)
#    def test_over_minoriy_is_not_ref(self, mut):
#        assume(sum(mut['AO']) / mut['DP'] > 0.2)
#        result = call_many(0, 80, mut)[1]
#        self.assertNotEquals(result, mut['ref'])

class ConsesusExampleTest(unittest.TestCase):
    def test_make_consensus_example(self):
        muts = [('CG', 'TT', 1),
                ('C', 'A', 5)]
        ref = 'ACGTACGT'
        expected = 'ATTTAAGT'
        actual = make_consensus(ref, muts)
        self.assertEquals(expected, actual)

    def test_single_example(self):
        muts = [{
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
        ref = make_seqrec('X', 'ACGTACGT')
        expected = 'ATTTAAGT'
        result = just_ref([ref], muts, 10, 80)
        self.assertEquals(expected, result)

from collections import Counter
countof = lambda c: lambda x: Counter(x).get(c, 0)
class ConsensusHypothesisTest(unittest.TestCase):
    #@st.random_module
    @given(ref_with_vcf_dicts_strategy_factory(), st.random_module())
    def test_n_count(self, ref_and_muts, rand):
        ref, muts = ref_and_muts
        originalNs = countof('N')(ref)
        assume(not filter(lambda x: 'N' in x, muts))
        expectedNs = len(filter(_['DP'] < 10, muts))  + originalNs
        result = just_ref([ref], muts, 10, 80)
        self.assertEquals(countof('N')(result), expectedNs)




