#!/usr/bin/env bash

export JIP_PATH=$PWD/jip_modules

# Super simple test to make sure pieces all work nicely together
tests/jip_modules/map_pipe.jip tests/testinput/f1.fastq tests/testinput/f1.fastq -r tests/testinput/ref.fasta --indexqual 21 --indexes tests/testinput/f1.index.fastq tests/testinput/f2.index.fastq --removebases 1
