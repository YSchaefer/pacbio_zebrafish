#!/bin/bash


cat sample_names.txt | \
while read fish; do \

        echo Map the HiFi reads from ${fish} to the filtered cluster representatives; \

        pbmm2 align \
        --sort --preset HiFi --sample ${fish} -x 95 \
        ../filtered_clusters.fasta \
        /home/yschaefer/unmapped_reads/bam_files/${fish}.ccs.nochim.bam \
        > ${fish}.ccs.filtered_clusters.bam; \


	# Filters:
	# 95% Identity (set in the mapping step)
	# Primary alignment (SAM flag 0 or 16)
	# Perfect mapping quality of 60
	# No more than 9 soft-clipping bases at the ends of the aligned reads
	# No reads with length less than 1 kb allowed

        echo Calculate the read depths for the matched clusters of ${fish}; \

        samtools view ${fish}.ccs.filtered_clusters.bam | \
        awk '($2 == 0 || $2 == 16) && $5 == 60' | grep -v '[0-9][0-9]S' | \
        while read line; do \

                read=$(echo ${line} | awk '{print $1}'); \
                read_length=$(grep "${fish}" /home/yschaefer/unmapped_reads/fasta_files/read_lengths.tbl | grep "${read}" | awk '{print $3}'); \

                if [[ ${read_length} -lt 1000 ]]; then \
                        :; \
                else \
                        echo ${line}; \
                fi; \

        done | \
        awk '{print $3}' | uniq -c | awk '{print $2, $1}' | sed '1s/^/Cluster Reads\n/' | sed 's/cluster_//' | column -t -s ' ' | sort -k1,1n \
        > ${fish}.absolute_depths.tbl; \

        echo Add the relative NLR read depths to the depth table of ${fish}; \

        NLR_reads=$(grep "${fish}" /import/yschaefer/pacbio_zebrafish/mappings/ccs_GRCz11/reads_mapped_to_ref_NLRs.tbl | awk '{print $2}'); \

        sed '1d' ${fish}.absolute_depths.tbl | awk -v NLR_reads="${NLR_reads}" '{print $1, $2, $2/NLR_reads}' | \
        sed '1s/^/Cluster Absolute_Depth Relative_Depth\n/' | column -t -s ' ' \
        > ${fish}.read_depths.tbl; \

        rm ${fish}.absolute_depths.tbl; \

done


# Filter for clusters with reads mapped to them

cat *.read_depths.tbl | \
sed '/Cluster/d' | awk '{print $1}' | sort -n | uniq | \
while read cluster_number; do \

	grep -A1 "^>cluster_${cluster_number}$" filtered_clusters.fasta; \

done > confirmed_clusters.fasta
