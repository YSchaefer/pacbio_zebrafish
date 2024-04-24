#!/bin/bash

cat sample_names.txt | while read fish; do \

	echo Assemble the CCS reads of fish ${fish}; \

	mkdir hifiasm/${fish}; \

	hifiasm -o hifiasm/${fish}/${fish} /home/yschaefer/unmapped_reads/fastq_files/${fish}.ccs.nochim.fastq; \

	echo The CCS reads from fish ${fish} were assembled; \


	cd hifiasm/${fish}/

	echo Extract the contigs of haplotype 1 from fish ${fish}; \

	grep '^S' ${fish}.bp.hap1.p_ctg.gfa | awk '{print $2, "\t" $3}' | sed 's/h1/>h1/' | sed 's/\t/\n/' \
	> ${fish}.hap1.fasta; \


	echo Extract the contigs of haplotype 2 from fish ${fish}; \

	grep '^S' ${fish}.bp.hap2.p_ctg.gfa | awk '{print $2, "\t" $3}' | sed 's/h2/>h2/' | sed 's/\t/\n/' \
	> ${fish}.hap2.fasta; \

	cat ${fish}.hap1.fasta ${fish}.hap2.fasta > ${fish}.fasta; \

done
