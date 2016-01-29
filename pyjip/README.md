# Simple jip pipeline

[Jip Docs](http://pyjip.readthedocs.org/en/latest)

# Dirty install

Would be fun later to maybe make a task to install executables inside a jip task
using conda.

For now do it manually
```
conda install bwa samtools cutadapt pyzmq
pip install pyjip
```

# Execute pipeline

For the examples/tests I ran I used the 947 sample that comes with the ngs_mapper pipeline


## Single threaded(non-piped)

```
./simplepipe.jip -r ../../functional/947.ref.fasta -f ../../functional/947/947_S32_L001_R{1,2}_001_2013_12_17.fastq
```

### Dry run

```
./simplepipe.jip -r ../../functional/947.ref.fasta -f ../../functional/947/947_S32_L001_R{1,2}_001_2013_12_17.fastq -- --show --dry
```

## Multy threaded(named pipe)

```
jip server &
cat <<EOF > ~/.jip/jip.json
    {
        "cluster": "jip.grids.JIP",
        "jip_cluster":{
            "port": 5556
        }
    }
EOF

jip submit simplepipe.jip -r ../../functional/947.ref.fasta -f ../../functional/947/947_S32_L001_R{1,2}_001_2013_12_17.fastq --pipe
```
