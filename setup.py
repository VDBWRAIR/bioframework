from setuptools import setup, find_packages

import bioframework

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
)
