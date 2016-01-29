# Simple jip pipeline

[Jip Docs](http://pyjip.readthedocs.org/en/latest)

# Dirty install

Would be fun later to maybe make a task to install executables inside a jip task
using conda.

For now do it manually
```
conda install bwa samtools cutadapt pyzmq
pip install pyjip
wget https://github.com/gnuplot/gnuplot/archive/Release_4_6_0.tar.gz -O- | tar xzvf -
cd gnuplot*
./prepare
./configure --prefix=$(cd $(dirname $(dirname $(which conda))) && pwd)
make && make install
```

# Execute pipeline

For the examples/tests I ran I used the 947 sample that comes with the ngs_mapper pipeline

For now you have to supply `-i` at the end or it tries to submit multiple jobs for each fastq

```
./simplepipe.jip -r ../../functional/947.ref.fasta -f ../../functional/947/947_S32_L001_R{1,2}_001_2013_12_17.fastq -i
```

## Dry run

```
./simplepipe.jip -r ../../functional/947.ref.fasta -f ../../functional/947/947_S32_L001_R{1,2}_001_2013_12_17.fastq -- --show --dry
```
