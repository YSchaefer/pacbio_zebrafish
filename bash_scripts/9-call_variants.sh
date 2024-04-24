#!/bin/bash


## Define the working directories and the bin version


BIN_VERSION="latest"
INPUT_DIR="/import/yschaefer/pacbio_zebrafish/fish_information/hifiasm/pan_NLRome_analysis/clusters_with_contigs/blast/mapping"
OUTPUT_DIR="/import/yschaefer/pacbio_zebrafish/fish_information/hifiasm/pan_NLRome_analysis/clusters_with_contigs/blast/mapping/variant_calls"


cat /import/yschaefer/pacbio_zebrafish/fish_information/sample_names.txt | \
while read fish; do \

	echo Index the bam file of ${fish}; \

	samtools index ${fish}.ccs.filtered_clusters.bam; \

done


samtools faidx filtered_clusters.fasta; \


echo Index files for all reference fastas and mapping files were created
echo Next is the variant calling


cat /import/yschaefer/pacbio_zebrafish/fish_information/sample_names.txt | \
while read fish; do \

	docker run --user "$(id -u):$(id -g)" \
	-v "${INPUT_DIR}:/input" \
	-v "${OUTPUT_DIR}:/output" \
	google/deepvariant:"${BIN_VERSION}" /opt/deepvariant/bin/run_deepvariant \
	--model_type=PACBIO \
	--ref=/input/filtered_clusters.fasta \
	--reads=/input/${fish}.ccs.filtered_clusters.bam \
	--output_vcf=/output/${fish}.ccs.FISNACHT_clusters.vcf.gz \
	--output_gvcf=/output/${fish}.ccs.FISNACHT_clusters.g.vcf.gz \
	--num_shards=10 \
	2> log.txt; \

done
