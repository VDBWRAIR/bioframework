from typing import Union, Dict, List, NamedTuple, Iterator
_Record = NamedTuple('VCFRecord', [("ALT", Union[str, List[str]]), ("REF", str), ("POS", int), ("CHROM", str), ("INFO", Dict[str, Union[int, List[int]]])]
)
