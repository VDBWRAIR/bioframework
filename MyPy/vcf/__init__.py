from typing import Union, Dict, List, NamedTuple, Iterator

#fields = [("ALT", Union[str, List[str]]), ("REF", str), ("POS", int), ("CHROM", str), ("INFO", Dict[str, Union[int, List[int]]])]
#
#VCFRecord = NamedTuple('VCFRecord', fields)

VCFRecord = NamedTuple('VCFRecord', [("ALT", Union[str, List[str]]), ("REF", str), ("POS", int), ("CHROM", str), ("INFO", Dict[str, Union[int, List[int]]])]
)
def Reader(s): # type: (str) -> Iterator[VCFRecord]
    pass
