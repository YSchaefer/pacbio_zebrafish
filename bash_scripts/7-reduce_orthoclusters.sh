#!/bin/bash


# Blast the representative sequences of the clusters against each other

blastn \
 -query all_representatives.fasta \
 -subject all_representatives.fasta \
 -outfmt 6 \
> reps_vs_reps.blast.tsv


# First, remove any alignments with less than 98% identity and less than 3 kb alignment length
# We are looking for very similar sequences and at least 1 kb of intronic match


awk '$1 != $2' reps_vs_reps.blast.tsv | awk '$3 >= 98 && $4 >= 3000 {print $1, $2, $3, $4, $11, $12}' | \
sed '1s/^/Cluster_1 Cluster_2 Identity Alignment_Length E_Value Bit_Score\n/' | column -t -s ' ' \
> reps_vs_reps.blast.tbl


# Next, find out how many percent of the representative sequences of the clusters are matched in the alignment
# If less than 98% are matched, remove the alignment
# This way, we make sure that even if representatives have similar parts but are not from the same gene copy,
# we keep them separate

cat ../representative_lengths.txt | \
while read line; do \

        cluster=$(echo ${line} | awk '{print $1}'); \
        rep_length=$(echo ${line} | awk '{print $2}'); \

        grep "^${cluster} " reps_vs_reps.blast.tbl | \
        while read blast_result; do \

                second_cluster=$(echo ${blast_result} | awk '{print $2}'); \
                alignment_length=$(echo ${blast_result} | awk '{print $4}'); \

                percent=$(echo "scale=2; " ${alignment_length} / ${rep_length} " * 100" | bc | sed 's/\.00$//'); \
                if [[ ${percent} -gt 97 ]]; then \
                        echo ${blast_result} ${rep_length} ${percent}; \
                fi; \

        done; \

done | \
sed '1s/^/Cluster_1 Cluster_2 Identity Alignment_Length E_Value Bit_Score Rep_Length Alignment_Percentage\n/' | column -t -s ' ' \
> filtered_reps_vs_reps.blast.tbl

# In the resulting file, there are
# a) clusters with very long representative sequences which have one or few other clusters fitting to them
# b) clusters with rather short representatives sequences which have one or few other clusters fitting to them
# c) clusters with rather short representative sequences which have a bunch of matches
# We should combine clusters from a) and b) into single clusters while clusters from group c) should be taken out of the dataset
# If a cluster representative fits one or few other clusters very well, they probably originate from the same gene copy
# If a cluster representative is short and matches lots of other clusters well, it probably belongs to one of the other clusters but
# a decision cannot be made because of the short length. It is therefore redundant and can be removed.


# Find the groups of clusters

sed '1d' filtered_reps_vs_reps.blast.tbl | awk '{print $1}' | uniq | \
while read cluster; do \

	combine_with=$(grep "^${cluster} " filtered_reps_vs_reps.blast.tbl | awk '{print $1, $2}' | sed "s/${cluster}[ $]//"); \

	echo ${cluster} ${combine_with}; \

done > groups_of_clusters.txt


# Combine the contigs of the clusters in the groups into a fasta file and create new consensus sequences

cat groups_of_clusters.txt | \
while read line; do \

	sum=$(echo ${line} | sed 's/ cluster_/+/g' | sed 's/^cluster_//' | bc); \
	rm cluster_${sum}n_contigs.fa; \

	echo ${line} | tr " " "\n" | sed 's/cluster_//' | \
	while read cluster_number; do \

		cat ../${cluster_number}_contig*.fa* >> cluster_${sum}n_contigs.fa; \

	done; \

done


ls cluster_*n_contigs.fa | sed 's/_contigs.fa//' > new_clusters.txt

cat new_clusters.txt | \
while read cluster; do \

	mafft --auto ${cluster}_contigs.fa > ${cluster}_contigs.aln.fa; \

	cons -sequence ${cluster}_contigs.aln.fa -outseq ${cluster}_contigs.consensus.fa; \

	sed -i "s/EMBOSS_001/${cluster}/" ${cluster}_contigs.consensus.fa; \

	rep_length=$(grep -v '^>' ${cluster}_contigs.consensus.fa | wc -c); \

	echo ${cluster} ${rep_length}; \

done > new_rep_lengths.tsv

rm *_contigs.aln.fa
rm new_clusters.txt


# Combine the new consensus sequences into one file and blast them against each other
# to find the new consensus sequences made from the same or similar groups of clusters

cat *_contigs.consensus.fa > new_consensus_sequences.fa

rm *_contigs.consensus.fa

makeblastdb -in new_consensus_sequences.fa -dbtype nucl -out new_consensus_sequences

blastn -query new_consensus_sequences.fa -db new_consensus_sequences -outfmt 6 \
> new_consensus_sequences.blast.tsv


# Filter the blast table again

awk '$1 != $2' new_consensus_sequences.blast.tsv | awk '$3 >= 98 && $4 >= 3000 {print $1, $2, $3, $4, $11, $12}' | \
sed '1s/^/Cluster_1 Cluster_2 Identity Alignment_Length E_Value Bit_Score\n/' | column -t -s ' ' \
> new_consensus_sequences.filtered_blast.tbl


# Repeat the alignment percentage step as well

cat new_rep_lengths.tsv | \
while read line; do \

        cluster=$(echo ${line} | awk '{print $1}'); \
        rep_length=$(echo ${line} | awk '{print $2}'); \

        grep "^${cluster} " new_consensus_sequences.filtered_blast.tbl | \
        while read blast_result; do \

                second_cluster=$(echo ${blast_result} | awk '{print $2}'); \
                alignment_length=$(echo ${blast_result} | awk '{print $4}'); \

                percent=$(echo "scale=2; " ${alignment_length} / ${rep_length} " * 100" | bc | sed 's/\.00$//'); \
                echo ${blast_result} ${rep_length} ${percent}; \

        done; \

done | \
sed '1s/^/Cluster_1 Cluster_2 Identity Alignment_Length E_Value Bit_Score Rep_Length Alignment_Percentage\n/' | column -t -s ' ' \
> new_consensus_sequences.alignment_percentages.tbl


# Save all new clusters

grep '^>' new_consensus_sequences.fa | sed 's/>//' | \
while read new_cluster; do \

	similar_cluster=$(grep "${new_cluster}" new_consensus_sequences.alignment_percentages.tbl | awk '$8 >= 98 {print $1, $2}' | sed "s/${new_cluster}//" | head -n 1); \
	if [[ -z ${similar_cluster} ]]; then \
		echo ${new_cluster}; \
	else \
		new_sum=$(echo ${new_cluster} ${similar_cluster} | sed 's/cluster_//g' | sed 's/n//g' | sed 's/ /+/' | bc); \
		echo Combine ${new_cluster} and ${similar_cluster} into cluster_${new_sum}nn; \
	fi; \

done > new_clusters.txt

rm new_consensus_sequences.blast.tsv
rm new_consensus_sequences.filtered_blast.tbl
rm new_consensus_sequences.alignment_percentages.tbl


# For the clusters to be combined, create new consensus sequences
# Put them to the other new cluster representatives

grep 'Combine' new_clusters.txt | \
while read line; do \

	cluster_1=$(echo ${line} | awk '{print $2}'); \
	cluster_2=$(echo ${line} | awk '{print $4}'); \
	latest_cluster=$(echo ${line} | awk '{print $6}'); \

	cat ${cluster_1}_contigs.fa ${cluster_2}_contigs.fa > ${latest_cluster}_contigs.fa; \

	mafft --auto ${latest_cluster}_contigs.fa > ${latest_cluster}_contigs.aln.fa; \

	cons -sequence ${latest_cluster}_contigs.aln.fa -outseq ${latest_cluster}.consensus.fa; \

	sed -i "s/EMBOSS_001/${latest_cluster}/" ${latest_cluster}.consensus.fa; \

done

cat cluster_*nn.consensus.fa > nn_clusters.consensus.fa
rm cluster_*nn.consensus.fa

rm *_contigs.fa
rm *_contigs.aln.fa

cat new_consensus_sequences.fa nn_clusters.consensus.fa | \
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | tr "\t" "\n" \
> latest_consensus_sequences.fasta

rm new_consensus_sequences.fa nn_clusters.consensus.fa


# Create a new file with the filtered clusters
# 1) Copy the previous file
# 2) Remove clusters with very short representatives (less than 3 kb)
# 3) Remove clusters which went into the combined clusters (from both combination steps)
# 4) Add the representatives of the combined clusters

cp all_representatives.correct_direction.fasta filtered_clusters.fasta


awk '$2 < 3000 {print $1}' ../representative_lengths.txt | \
while read cluster; do \

	sed -i "/${cluster}/{N;d;}" filtered_clusters.fasta; \

done

cat groups_of_clusters.txt | tr " " "\n" | sort | uniq | \
while read cluster; do \

	sed -i "/${cluster}/{N;d;}" filtered_clusters.fasta; \

done

awk '{print $NF}' new_clusters.txt | sort | uniq | \
while read new_cluster; do \

	grep -A1 "${new_cluster}" latest_consensus_sequences.fasta >> filtered_clusters.fasta; \

done

rm latest_consensus_sequences.fasta
