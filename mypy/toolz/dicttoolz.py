from typing import Dict, Any, Callable, TypeVar
K = TypeVar('K')
V = TypeVar('V')
V2 = TypeVar('V2')
V3 = TypeVar('V3')
def merge(d1, d2): # type: (Dict[K,V], Dict[K,V]) -> Dict[K,V]
    pass

def dissoc(d, k): # type: (Dict[K,V], K) -> Dict[K,V]
  pass

def merge_with(f, d1, d2): # type: (Callable[[V,V2], V3], Dict[K,V], Dict[K,V2]) -> Dict[K,V3]
    pass

def valfilter(f, d): # type: (Callable[[V], bool], Dict[K,V]) -> Dict[K,V]
  pass



#from typing import Dict, Any, Callable, TypeVar
#T = TypeVar('T')
#def merge(d1, d2): # type: (Dict[Any,Any], Dict[Any,Any]) -> Dict[Any,Any]
#    pass
#
#def dissoc(d, k): # type: (Dict[Any,Any], Any) -> Dict[Any,Any]
#  pass
#
#def merge_with(f, d1, d2): # type: (Callable, Dict[Any,Any], Dict[Any,Any]) -> Dict[Any,Any]
#    pass
#
#def valfilter(f, d): # type: (Callable, Dict[Any,Any]) -> Dict[Any,Any]
#  pass
