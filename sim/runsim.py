
ss="MS" # miseq platform
l=250 #length of reads
p=True # paired end
m=500 #mean size of DNA fragments of PE
s=10 stdev of fragment size (see m)
-f   --fcov     the fold of read coverage to be simulated or number of reads/read pairs generated for each amplicon
-c toltal # reads to be generated (excludes -f)
-f # create a zero-error sam file
-ir = insertsion rate
-dr = deletion rate # ther is also ir2 and dr2
mp = mate-pair simulation
rs = random seed!
sp = separate qual profiles for diff. bases
Linux-amd64/bin/pbsim --data-type CLR --depth 20 --model_qc data/model_qc_clr sample/sample.fasta
	  -qs2 --qShift2  the amount to shift every second-read quality score by
	                  NOTE: For -qs/-qs2 option, a positive number will shift up quality scores (the max is 93)

* ART by default selects a built-in quality score profile according to the read length specified for the run.
"-ss", "HS25", "-l", "150", "-f", "10", "-p", "-m", "500", "-s", "10", "-sam"]
