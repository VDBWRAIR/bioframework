from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import SeqIO
from StringIO import StringIO

import os
from os.path import *

TESTDIR = dirname(abspath(__file__))
PROJDIR = dirname(TESTDIR)

def make_dna_seqrec(id, desc, seq, qual=None):
    rec = SeqRecord(
        Seq(seq, ambiguous_dna),
        id=id,
        description=desc
    )
    if qual:
        rec.letter_annotations['phred_quality'] = qual
    return rec

def seqrec_stream(ids, descs, seqs, quals):
    for id, desc, seq, qual in zip(ids, descs, seqs, quals):
        yield make_dna_seqrec(str(id), str(desc), seq, qual)

def write_seq_stream(stream, output, format='fastq'):
    ''' Write seqrecords to output handle '''
    return SeqIO.write(stream, output, format)
