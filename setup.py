from glob import glob
from os.path import join
import sys

from setuptools import setup, find_packages

import bioframework

def jip_modules(path=bioframework.JIP_PATH):
    return glob(join(path, '*.jip'))

setup(
    name = bioframework.__projectname__,
    version = bioframework.__release__,
    packages = find_packages(),
    author = bioframework.__authors__,
    author_email = bioframework.__authoremails__,
    description = bioframework.__description__,
    license = "GPLv2",
    keywords = bioframework.__keywords__,
    entry_points = {
    },
    install_requires = [
        'pyjip',
        'toolz',
    ],
    data_files=[
        (join(sys.prefix,'bin'), jip_modules()),
    ]
)
