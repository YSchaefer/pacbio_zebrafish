#!/bin/bash


transeq -sequence ../WGA/${assembly}.fasta -clean -frame F -outseq ../WGA/${assembly}.fwd_frames.fa; \

/home/yschaefer/scripts/bin/hmmsearch \
--nobias --nonull2 \
-o ../hmmsearch/${assembly}.fwd.FISNACHT_search.output \
--domtblout ../hmmsearch/${assembly}.fwd.FISNACHT_search.domtbl \
../HMMs/zf_FISNA-NACHT.prot.hmm \
../WGA/${assembly}.fwd_frames.fa; \

rm ../WGA/${assembly}.fwd_frames.fa; \

echo The forward FISNACHT exons from assembly ${assembly} were searched for; \

transeq -sequence ../WGA/${assembly}.revcomp.fasta -clean -frame F -outseq ../WGA/${assembly}.rev_frames.fa; \

/home/yschaefer/scripts/bin/hmmsearch \
--nobias --nonull2 \
-o ../hmmsearch/${assembly}.rev.FISNACHT_search.output \
--domtblout ../hmmsearch/${assembly}.rev.FISNACHT_search.domtbl \
../HMMs/zf_FISNA-NACHT.prot.hmm \
../WGA/${assembly}.rev_frames.fa; \

rm ../WGA/${assembly}.rev_frames.fa; \

echo The reverse FISNACHT exons from assembly ${assembly} were searched for; \
