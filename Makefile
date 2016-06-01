include gmsl

ILLUMINA_URL := http://www.niehs.nih.gov/research/resources/assets/docs/artbinchocolatecherrycake031915linux64tgz.tgz
PBSIM_URL := https://pbsim.googlecode.com/files/pbsim-1.0.3-Linux-amd64.tar.gz
DIR_FROM_URL := $(basename $(call last, $(call split, "/", $(1))))
#ILLUMINA_DIR := $(call DIR_FROM_URL, ILLUMINA_URL)
#PBSIM_DIR := $(call DIR_FROM_URL, PBSIM_URL)

pbsim-1.0.3-Linux-amd64.tar.gz:
	wget $(PBSIM_URL)

pbsim-1.0.3-Linux-amd64: pbsim-1.0.3-Linux-amd64.tar.gz
	tar -xf $<
artbinchocolatecherrycake031915linux64tgz.tgz: 
	wget $(ILLUMINA_URL)


artbinchocolatecherrycake031915linux64tgz:
define wget_and_untar

$(basename $(call last, $(call split,/, $(1)))):  $(call last, $(call split,/, $(1))) 
	tar -xf $(call last, $(call split,/, $(1)))

$(call last, $(call split,/, $(1))):
	wget $(1)
	tar -xf $(call last, $(call split,/, $(1)))
endef
s := /
#define wget_and_untar
#$(basename 2): 
#	wget $(1)$(2)
#	tar xf $(2) 
#endef

bin/art_illumina: 
	echo $(call wget_and_untar, $(ILLUMINA_URL))

#bin/pbsim:
#	$(call wget_and_untar, PBSIM_URL)

