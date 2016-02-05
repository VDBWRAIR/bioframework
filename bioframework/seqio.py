import itertools
from Bio import SeqIO
from functools import partial
def interleave(iter1, iter2):
    for forward, reverse in itertools.izip(iter1, iter2):
        assert forward.id == reverse.id, "%s did not match %s" % \
            (forward.id, reverse.id)
        yield forward
        yield reverse

def write_zip_results(func, output, out_format, *files):
    parse_fastq = partial(SeqIO.parse, format='fastq')
    inputs = map(parse_fastq, files)
    results = func(*inputs)
    count = SeqIO.write(results, output, out_format)

def paired_to_interleave(forward, reverse, output, out_format):
    return write_zip_results(interleave, output, out_format, forward, reverse)




