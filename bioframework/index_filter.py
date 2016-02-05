from Bio import SeqIO
import re
import os
from itertools import ifilter, imap, izip
from toolz import compose
from toolz.itertoolz import second, first
from seqio import write_zip_results
from functools import partial
def filter_on_index(seqs, index_seqs, predicate):
    pred = compose(predicate, second)
    return imap(first, ifilter(pred, izip(seqs, index_seqs)))

def write_index_filter(input, output, predicate):
    def get_index(path):
        index_name = re.sub(r"_R([12])_", r"_I\1_", path)
        if os.path.exists(index_name): return index_name
    index = get_index(input)
    if index is None:
        raise ValueError("Index %s for file %s not found" % (index, input))
    return write_zip_results(partial(filter_on_index, predicate=predicate),
                      output, 'fastq', input, index)

def filter_on_index_quality(input, output, minimum):
    """Removes reads from a paired read file if its associated
    index file has a quality lower than minimum for that read index."""
    below_qual = lambda x: min(x.quality) < minimum
    write_index_filter(input, output, below_qual)
    return 0

