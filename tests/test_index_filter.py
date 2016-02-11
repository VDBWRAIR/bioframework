from . import unittest
from testfixtures import tempdir, TempDirectory

import common
from os.path import *
import os
import shutil

testname = splitext(basename(__file__))[0]
testdirbasepath = join(common.TESTDIR, 'testoutput', testname)

indexes = [
    [40]*10,
    [0]*10,
    [20]*10,
    [40]*4 + [0] + [40]*5,
]
num_indexes = len(indexes)
num_values = len(indexes[0])

check_index = lambda: common.seqrec_stream(range(1,num_indexes+1), range(1,num_indexes+1), ['A'*num_values]*num_indexes, indexes)
good_index = lambda: common.seqrec_stream(range(1,num_indexes+1), range(1,num_indexes+1), ['A'*num_values]*num_indexes, [[40]*num_values]*num_indexes)

class RunsInFixtureDir(object):
    def setUp(self):
        # Start new tempdir
        if exists(testdirbasepath):
            shutil.rmtree(testdirbasepath)
        self.tdir = TempDirectory(path=testdirbasepath)
        p = os.makedirs(testdirbasepath)

from bioframework import index_filter
class TestFiltersUsingIndex(RunsInFixtureDir, unittest.TestCase):
    def test_filters_both_pairs_when_forward_is_lowqual(self):
        i_read = common.seqrec_stream(range(1,9), range(1,9), ['A']*8, [[40]]*8)
        f_index = check_index()
        r_index = good_index()
        ifastq = self.tdir.write('i_read.fastq', '')
        rfastq = self.tdir.write('r_index.fastq', '')
        ffastq = self.tdir.write('f_index.fastq', '')
        ofastq = self.tdir.write('output.fastq', '')
        common.write_seq_stream(f_index, ffastq)
        common.write_seq_stream(r_index, rfastq)
        common.write_seq_stream(i_read, ifastq)
        c = index_filter.filter_on_index_quality_interleaved(ifastq, ffastq, rfastq, ofastq, 20)
        self.tdir.compare(sorted(['f_index.fastq','r_index.fastq','i_read.fastq', 'output.fastq']))
        self.assertEqual(4, c)

    def test_filters_both_pairs_when_reverse_is_lowqual(self):
        i_read = common.seqrec_stream(range(1,9), range(1,9), ['A']*8, [[40]]*8)
        f_index = good_index()
        r_index = check_index()
        ifastq = self.tdir.write('i_read.fastq', '')
        rfastq = self.tdir.write('r_index.fastq', '')
        ffastq = self.tdir.write('f_index.fastq', '')
        ofastq = self.tdir.write('output.fastq', '')
        common.write_seq_stream(f_index, ffastq)
        common.write_seq_stream(r_index, rfastq)
        common.write_seq_stream(i_read, ifastq)
        c = index_filter.filter_on_index_quality_interleaved(ifastq, ffastq, rfastq, ofastq, 20)
        self.tdir.compare(sorted(['f_index.fastq','r_index.fastq','i_read.fastq', 'output.fastq']))
        self.assertEqual(4, c)
