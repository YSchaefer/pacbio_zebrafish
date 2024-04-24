#!/bin/bash

cat /import/yschaefer/pacbio_zebrafish/fish_information/sample_names.txt | \
while read fish; do \

	echo Map the CCS reads of ${fish} to the reference genome; \

	pbmm2 align \
        --preset ccs \
        --sort \
        /import/yschaefer/pacbio_zebrafish/WGA/GRCz11.fasta \
        /home/yschaefer/unmapped_reads/fastq_files/${fish}.ccs.nochim.fastq > \
        ../mappings/ccs_GRCz11/${fish}.ccs.GRCz11.bam; \


	echo Assign the reads of ${fish} to the targeted genes; \

	samtools view ${fish}.ccs.GRCz11.bam | awk '$2 == 0 || $2 == 16' | \
	while read line; do \

		read=$(echo ${line} | awk '{print $1}'); \

		chr_identifier=$(echo ${line} | awk '{print $3}'); \
		chromosome=$(grep 'GRCz11' /import/yschaefer/pacbio_zebrafish/WGA/assembly_chromosome_header.tbl | grep "${chr_identifier}" | awk '{print $2}'); \
		# assembly_chromosome_header.tbl is a file with chromosome identifiers and what chromosome they refer to

		position=$(echo ${line} | awk '{print $4}'); \

		gene=$(grep "^${chromosome} " /import/yschaefer/pacbio_zebrafish/targets/GRCz11_targets_and_NLRs.bed | \
		awk -v position="${position}" '$2 < position && $3 > position {print $4}' | tr "\n" "," | sed 's/,$//'); \
		# The bed file GRCz11_targets_and_NLRs.bed contains the genomic positions of all targeted genes and all reference NLR genes from biomart

		if [[ -z ${gene} ]]; then \
			gene=$(echo not_in_a_targeted_gene); \
		fi; \

		echo ${read} ${chromosome} ${position} ${gene}; \

	done | sed '1s/^/Read Chromosome Position Gene\n/' | column -t -s ' ' \
	> ${fish}.reads_and_genes.tbl; \

done
