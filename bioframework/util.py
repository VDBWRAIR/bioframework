import sys

import os.path

def normalize_handle(option, default=sys.stdin, mode='r'):
    '''
    Get file handle for file option
    '''
    # Normalize input to handle
    if option.get():
        return open(option.get(), mode)
    else:
        return default

def format_from_ext(option_or_path):
    '''
    Return file format from option or path

    >>> format_from_ext('/path/file.format')
    'format'
    >>> format_from_ext('file.format')
    'format'
    >>> format_from_ext('file')
    Traceback (most recent call last):
    ...
    ValueError: no extension
    '''
    pth = option_or_path
    if hasattr(option_or_path, 'get'):
        pth = option_or_path.get()
    path, ext = os.path.splitext(pth)
    if not ext.startswith('.'):
        raise ValueError('no extension')
    return ext[1:]
