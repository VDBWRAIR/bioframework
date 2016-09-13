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
    # scripts = glob('scripts/*'),
    scripts = ['scripts/consensus.sh'],
    entry_points = {
        'console_scripts': [
            'fb_tagreads = bioframework.tagreads:main',
            'fb_consensus = bioframework.consensus:main'
        ]
    },
    install_requires = [
        'pyjip',
        'toolz',
        'docopt',
        'schema',
        'pyvcf',
        'typing'
    ],
    data_files=[
        (join(sys.prefix,'bin'), jip_modules()),
    ]
)

