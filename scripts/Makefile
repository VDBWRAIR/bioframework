
SAMPLE=$(notdir $(NGSDIR))
BAM=$(NGSDIR)/$(SAMPLE).bam
REF=$(NGSDIR)/$(SAMPLE).ref.fasta
TAGGED=$(NGSDIR)/tagged.bam
FBVCF=$(NGSDIR)/freebayes.vcf
FBCONSENSUS=$(NGSDIR)/$(SAMPLE).freebayes_consensus.fasta

all: $(FBCONSENSUS)

$(TAGGED): $(BAM)
	fb_tagreads $< | samtools view  - -Sbh | samtools sort - $(TAGGED:.bam=)

$(TAGGED).bai: $(TAGGED)
	samtools index $<

$(FBVCF): $(TAGGED) $(NGSDIR)/tagged.bam.bai $(REF)
	freebayes -f $(REF) $< > $@

$(FBCONSENSUS):  $(FBVCF) $(REF)
	fb_consensus --vcf $< --ref $(REF) --sample $(SAMPLE) > $@
#cd ~/bioframes && python bioframes/consensus.py $^ > $@
