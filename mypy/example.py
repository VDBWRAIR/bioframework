from typing import List, Dict, Generator, Iterator, Iterable, Tuple
from Bio import SeqIO
from itertools import imap
from Bio.SeqRecord import SeqRecord
def test_long(): # type: () -> int
    return 11999999L
def test_seqIO_map_fails(s): # type: (str) -> List[SeqRecord]
    return map(lambda x: x.id, SeqIO.parse(s))

#def test_seqIO_map_fails2(s): # type: (str) -> Iterator[SeqRecord]
#    return map(lambda x: x.id, SeqIO.parse(s))
def test_seqIO_map_passes(s): # type: (str) -> Iterable[str]
    return imap(lambda x: x.id, SeqIO.parse(s))

def test_seqIO(s): # type: (str) -> Iterator[SeqRecord]
    return SeqIO.parse(s)
def test_list_seqIO(s): # type: (str) -> List[SeqRecord]
    return list(SeqIO.parse(s))
def test_seqIO_fails(s): # type: (str) -> List[str]
    return SeqIO.parse(s)
def test_should_pass(s): # type: (SeqRecord) -> str
    return s.id
def test_should_fail(s): # type: (SeqRecord) -> int
    return s.id
#def test_should_fail(): # type: () -> List[SeqRecord]
#    return 3

#a = test_should_fail()
def test_ordered_dict(od): # type: (Dict[str,int]) -> Dict[str,int]
    return 1   #type error 1
#
#a = test_ordered_dict(1)   #type error 2
#
#def test_me():
#    a = test_ordered_dict(1)  # type error 3 is not reported

####def test_ordered_dict(od: typing.Dict[str,int]) -> typing.Dict[str,int]:
####    return 1   #type error 1
####
####a = test_ordered_dict(1)   #type error 2
####
####def test_me():
####    a = test_ordered_dict(1)  # type error 3 is not reported
###
