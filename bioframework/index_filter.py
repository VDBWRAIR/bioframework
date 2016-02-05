from Bio import SeqIO
from toolz.itertoolz import second, first, partition
import re
import os
from itertools import ifilter, imap, izip
from toolz import compose
from seqio import write_zip_results
from functools import partial

def filter_on_index(predicate, seqs, index_seqs):
    pred = compose(predicate, second)
    return imap(first, ifilter(pred, izip(seqs, index_seqs)))

def write_index_filter(input, output, predicate):
    def get_index(path):
        index_name = re.sub(r"_R([12])_", r"_I\1_", path)
        if os.path.exists(index_name): return index_name
    index = get_index(input)
    if index is None:
        raise ValueError("Index %s for file %s not found" % (index, input))
    return write_zip_results(partial(filter_on_index, predicate),
                      output, 'fastq', input, index)

def filter_on_index_quality(input, output, minimum):
    """Removes reads from a paired read file if its associated
    index file has a quality lower than minimum for that read index."""
    below_qual = lambda x: min(x.quality) < minimum
    write_index_filter(input, output, below_qual)
    return 0


def filter_on_index_quality_interleaved(interleaved, index1, index2, output, minimum):
    '''enpair the interleaved read file,zip that enpaired with the two indexes
    drop pairs from the interleaved file if either *index* is below the minimum'''
    above_min = lambda x: min(x.quality) >= minimum
    def indexes_above_min(seqs, idx1,idx2):
        return above_min(idx1) and above_min(idx2)
    def qual_filter(interleaved, idx1, idx2):
        pairs = partial(partition, 2)
        interleaved = pairs(interleaved)
        zipped = izip(interleaved, idx1, idx2)
        filtered = imap(first, ifilter(indexes_above_min, zipped))
        return izip(*filtered)
    return write_zip_results(qual_filter, output, 'fastq')
