
BAM=$(NGSDIR)/$(notdir $(NGSDIR)).bam
REF=$(NGSDIR)/$(notdir $(NGSDIR)).ref.fasta
TAGGED=$(NGSDIR)/tagged.bam
FBVCF=$(NGSDIR)/freebayes.vcf
FBCONSENSUS=$(NGSDIR)/freebayes_consensus.fasta

all: $(FBCONSENSUS)

$(TAGGED): $(BAM)
	python tagreads.py $< | samtools view  - -Sbh | samtools sort - $(TAGGED:.bam=)

$(TAGGED).bai: $(TAGGED)
	samtools index $<

$(FBVCF): $(TAGGED) $(NGSDIR)/tagged.bam.bai $(REF)
	freebayes -f $(REF) $< > $@

$(FBCONSENSUS): $(BAM) $(REF) $(FBVCF) 
	python consensus.py --vcf $< --ref $(REF) > $@
	#cd ~/bioframes && python bioframes/consensus.py $^ > $@
