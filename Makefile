ILLUMINA_URL = "http://www.niehs.nih.gov/research/resources/assets/docs/artbinchocolatecherrycake031915linux64tgz.tgz" PBSIM_URL = "https://pbsim.googlecode.com/files/pbsim-1.0.3-Linux-amd64.tar.gz"
define wget_and_untar
$(basename 2): 
	wget $(1)$(2)
	tar xf $(2)

endef
bin/art_illumina: 
	$(call wget_and_untar, ILLUMINA_URL, 


bin/pbsim:

