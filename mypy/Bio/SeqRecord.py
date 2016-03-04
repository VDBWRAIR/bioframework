# from Bio.SeqIO import SeqIO
from typing import NamedTuple
class Stringable(object):
    def __str__(self): # type: () -> str
        pass
SeqRecord = NamedTuple('SeqRecord', [('id', str), ('seq', Stringable)])
