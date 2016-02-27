from biottest.biohypothesis import ref_with_vcf_dicts_strategy_factory
from hypothesis import given
from hypothesis import strategies as st
from hypothesis import given, assume
from bioframework.consensus import call_many, all_consensuses
import unittest

simple_vcf_dict_strategy = st.tuples(st.text(string.ascii_letters),
                     st.integers(min_value=1),
                     st.text(string.ascii_letters)) \
                     .flatmap(lambda tup:
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
        self.assertEquals(result == mut['ref'])

    @given(simple_vcf_dict_strategy)
    def test_over_minoriy_is_alt(self, mut):
        assume(sum(mut['AO']) / mut['DP'] > 0.2)
        assume(len(mut['alt']) == 1)
        result = call_many(0, 80, mut)[1]
        self.assertEquals(result == mut['alt'])

class ConsesusExampleTest(unittest.TestCase):
    def test_single_example(self):
        muts = [{
            'pos' : 2,
            'ref' : 'CG',
            'alt' : 'TT',
            'AO'  : 99,
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
    ref = 'ACGTACGT'
    expected = 'ATTTAAGT'
    result = just_ref([ref], muts, 10, 80)
    self.assertEquals(expected, result)

class ConsensusHypothesisTest(unittest.TestCase):
    from collections import Counter


