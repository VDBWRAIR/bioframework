import itertools
from Bio import SeqIO

def interleave(iter1, iter2):
    for forward, reverse in itertools.izip(iter1, iter2):
        assert forward.id == reverse.id, "%s did not match %s" % \
            (forward.id, reverse.id)
        yield forward
        yield reverse

def write_zip_results(func, file1, file2, output, out_format):
    f1, f2 = open(file1), open(file2)
    records = func(SeqIO.parse(f1, 'fastq'), SeqIO.parse(f2, 'fastq'))
    outfile = open(output, 'w')
    count = SeqIO.write(records, outfile, out_format)
    outfile.close()

def paired_to_interleave(forward, reverse, output, out_format):
    return write_zip_results(interleave, forward, reverse, output, out_format)




