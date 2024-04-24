#!/bin/bash


cat sample_names.txt | while read fish; do \

	echo Translate all contigs of haplotype 1 in all frames; \

        transeq -sequence ${fish}.hap1.fasta -clean -frame 6 -outseq ${fish}.hap1.all_frames.fa; \

        echo Translate all contigs of haplotype 2 in all frames; \

        transeq -sequence ${fish}.hap2.fasta -clean -frame 6 -outseq ${fish}.hap2.all_frames.fa; \


        echo Next, search the translated contigs for NLR domains; \

        /home/yschaefer/scripts/bin/hmmsearch \
        --nobias --nonull2 \
        -o ${fish}.hap1.FISNACHT_search.output \
        --domtbl ${fish}.hap1.FISNACHT_search.domtbl \
        /import/yschaefer/pacbio_zebrafish/HMMs/zf_FISNA-NACHT.prot.hmm \
        ${fish}.hap1.all_frames.fa; \

        /home/yschaefer/scripts/bin/hmmsearch \
        --nobias --nonull2 \
        -o ${fish}.hap2.FISNACHT_search.output \
        --domtbl ${fish}.hap2.FISNACHT_search.domtbl \
        /import/yschaefer/pacbio_zebrafish/HMMs/zf_FISNA-NACHT.prot.hmm \
        ${fish}.hap2.all_frames.fa; \


        /home/yschaefer/scripts/bin/hmmsearch \
        --nobias --nonull2 \
        -o ${fish}.hap1.B30.2_search.output \
        --domtbl ${fish}.hap1.B30.2_search.domtbl \
        /import/yschaefer/pacbio_zebrafish/HMMs/zf_B30.2.prot.hmm \
        ${fish}.hap1.all_frames.fa; \

        /home/yschaefer/scripts/bin/hmmsearch \
        --nobias --nonull2 \
        -o ${fish}.hap2.B30.2_search.output \
        --domtbl ${fish}.hap2.B30.2_search.domtbl \
        /import/yschaefer/pacbio_zebrafish/HMMs/zf_B30.2.prot.hmm \
        ${fish}.hap2.all_frames.fa; \

	echo Fish ${fish} was searched for NLR domains; \

done
