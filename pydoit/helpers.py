import os
from glob import glob# ls =  partial(glob.glob1, ".")
from functools import wraps


class F (str):
  """represents a file"""
  def __div__(self, other):
      return F(os.path.join(self, other))

  def __rdiv__(self, other):
      return F(os.path.join(other, self))

  def __getattr__(self, ext):
    try:
      return super(self.__class__, self).__getattr__(ext)
    except:
      return F('{0}.{1}'.format(self, ext))

def assoc_all(d, **kwargs):
    dict = d.copy()
    for k,v in kwargs.items():
        dict[k] = v
    return dict

def update_all(d, **kwargs):
    dict = d.copy()
    for k,f in kwargs.items():
        dict[k] = f(dict[k])
    return dict

class D (dict):
  """dictionary with extra methods"""

  def assoc(self, **kwargs):
    return assoc_all(self, **kwargs)

  def apply(self, **kwargs):
    return update_all(self, **kwargs)

def to_fifo(func, *argnames):  #works as composition
    code = func.func_code
    names = code.co_varnames[:code.co_argcount]
    @wraps(func)
    def decorated(*args,**kwargs):
        for argname in argnames:
            argval = kwargs.get(argname, args[names.index(argname)])
#            try:
#                argval = args[names.index(argname)]
#            except ValueError:
#                argval = kwargs[argname]
            if hasattr(argval, '__iter__'):
                for f in argval:
                    os.mkfifo(f)
            else:
                os.mkfifo(argval)

        return func(*args, **kwargs)
    return decorated

#decorator
def fifo(*argnames):
    def decorated(func):
        return to_fifo(func, *argnames)
    return decorated

fwd, rev = glob("*_R1_*.fastq"), glob("*_R2_*.fastq")
ext = lambda s: s.split('.')[-1]
swap_ext = lambda s, ext: '.'.join(s.split('.')[:-1] + [ext])


#{
#    cutadapt_up : filter_up,
#    filter_up : get_input_up,
#    cutadapt_paired : filter_paired,
#    filter_paired : get_input_paired,
#    bwa_up : cutadapt_up,
#    bwa_paired : cutadapt_paired,
#    samtools_merge : [bwa_up, bwa_paired],
#    samtools_index : samtools_merge,
#    freebayes : [samtools_index
#
