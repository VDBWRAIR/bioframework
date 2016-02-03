# BioFramework

[![Build Status](https://travis-ci.org/VDBWRAIR/bioframework.svg?branch=master)](https://travis-ci.org/VDBWRAIR/bioframework)

# Setup/Install

1. Install/Setup miniconda

   If you do not already have miniconda installed
   
   ```
   wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O conda.sh
   bash conda.sh -p bioframework
   export PATH=${PWD}/myconda/bin:$PATH
   conda config --add channels r
   conda config --add channels bioconda
   ```

2. Install project

   ```
   conda create -n bioframework pip
   source activate bioframework
   conda install --file requirements-conda.txt
   python setup.py install
   ```
   
Any time you open a new terminal you will need to activate the `bioconda` environment with: 

`source activate bioframework`
