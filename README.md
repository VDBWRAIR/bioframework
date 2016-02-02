# BioFramework

# Setup/Install

1. Install miniconda

   ```
   wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O conda.sh
   bash conda.sh -p myconda
   export PATH=${PWD}/myconda/bin:$PATH
   conda config --add channels r
   conda config --add channels bioconda
   ```

2. Install required software

   ```
   conda install --file requirements.txt
   ```
