#!/bin/bash

cd ../HMMs/

# Align the exons preceding the B30.2 exon from the reference genome

mafft --auto B30.2.special_exon.fasta > B30.2.special_exon.aln.fa

# Create an HMM and press it

hmmbuild B30.2.special_exon.hmm B30.2.special_exon.aln.fa

hmmpress B30.2.special_exon.hmm


# Search the NLR contigs of the fish with the model

cd ../fish_information/

cat sample_names.txt | \
while read fish; do \

	/home/yschaefer/scripts/bin/hmmsearch \
	--nobias --nonull2 \
	-o hifiasm/${fish}/${fish}.special_exon_search.output \
	--domtblout hifiasm/${fish}/${fish}.special_exon_search.domtbl \
	../HMMs/B30.2.special_exon.hmm \
	hifiasm/${fish}/${fish}.fasta; \

done
