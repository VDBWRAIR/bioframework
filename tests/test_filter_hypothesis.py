#TODO: make io test use commandline jip
#NOTE: this duplicates  a lot of the tested code's logic (see IO, above_min, imports) and API dependencies (see imports)
from Bio.SeqRecord import SeqRecord
import sh
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from functools import partial
import re
from bioframework.seqio import interleave
import unittest 
import operator
from toolz.itertoolz import partition
import hypothesis
from hypothesis import strategies as st
from hypothesis import find, note, given, assume, example
from bioframework.index_filter import qual_filter, filter_on_index_quality_interleaved

above_min = lambda minimum: lambda x: min(x.letter_annotations['phred_quality']) >= minimum
make_seqrec = lambda id, seq, quals: \
                SeqRecord(Seq(seq, IUPAC.ambiguous_dna), id=str(id), description='', letter_annotations={'phred_quality':quals})

rec = st.integers(min_value=1, max_value=10).flatmap(
    lambda n:
        st.builds(
            make_seqrec,
            st.integers(),
            st.text(alphabet='ATGCN', min_size=n, max_size=n),
            st.lists(st.integers(min_value=35, max_value=40), min_size=n, max_size=n)
        )
)

reads_and_indices = st.integers(min_value=1,max_value=5).flatmap(
   lambda n:
   st.tuples(
      st.integers(min_value=0, max_value=40),
      st.lists(rec, min_size=n, max_size=n).map(lambda x: list(interleave(x, x))),
      st.lists(rec, min_size=n, max_size=n),
      st.lists(rec, min_size=n, max_size=n),
   ))

def ilen(seq): return sum(1 for _ in seq)

class TestIndexQualityFilter(unittest.TestCase):
#    @given(reads_and_indices)
#    def test_idempotent(self, seqs):
#        self.assertSequenceEqual(list(qual_filter(*seqs)), list(qual_filter(list(qual_filter(*seqs)), *seqs[1:])))

    @given(reads_and_indices)
    def test_always_less_or_equal_size(self, seqs):
        self.assertLessEqual(ilen(qual_filter(*seqs)), ilen(seqs[1]))

    @given(reads_and_indices)
    def test_result_is_even(self, seqs):
        self.assertTrue(ilen(qual_filter(*seqs)) % 2 == 0)

    @given(reads_and_indices)
    def test_result_is_interleaved(self, seqs):
        result = qual_filter(*seqs)
        pairs = partition(2, result)
        matches = map(lambda (x,y): x.id == y.id, pairs)
        self.assertTrue(all(matches)) 

    @given(reads_and_indices)
    def test_drops_all_if_all_lower(self, seqs):
        (_min, reads, i1, i2) = seqs
        assume(not any(map(above_min(_min), i1)))
        result = qual_filter(*seqs)
        self.assertEqual([], list(result))

    @given(reads_and_indices)
    def test_drops_none_if_all_above(self, seqs):
        (_min, reads, i1, i2) = seqs
        assume(all(map(above_min(_min), i1)))
        assume(all(map(above_min(_min), i2)))
        result = qual_filter(*seqs)
        self.assertEqual(reads, list(result))

    @given(reads_and_indices)
    def test_drops_some_if_any_lower(self, seqs):
        (_min, reads, i1, i2) = seqs
        assume(not all(map(above_min(_min), i1 + i2)))
        result = qual_filter(*seqs)
        self.assertLess(ilen(result), ilen(reads))

    @given(reads_and_indices)
    #@hypothesis.settings(max_examples=3)
    def test_drops_some_if_any_lower_with_IO(self, seqs):
        (_min, reads, i1, i2) = seqs
        assume(not all(map(above_min(_min), i1 + i2)))
        names =  ['r', 'i1', 'i2']
        map(partial(SeqIO.write, format='fastq'), seqs[1:], names)
        outfile = 'foo'
        filter_on_index_quality_interleaved(names[0], names[1], names[2], outfile, _min)
        # version using the jip module itself not working
        #cmd = sh.Command('jip_modules/interleaved_index_filter.jip')
        #cmd(i=names[0], index1=names[1], index2=names[2], minimum=_min, o=outfile) # could use stdout
        result = SeqIO.parse(outfile, 'fastq')
        self.assertLess(ilen(result), ilen(reads))


if __name__ == '__main__':
    unittest.main()
# index files should be unchanged
# the number of reads left should be original - number of bad reads in each index which don't overlap
# can use find(strategy, boolean) 
#_not = lambda f: lambda x: not f(x)
