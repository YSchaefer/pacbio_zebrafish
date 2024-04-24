#!/bin/bash
SAMPFILE=$1

## SAMPFILE is a textfile containing a list of sample names, one per row


cat ${SAMPFILE} | while read i; do
        echo "Now processing "$i
        samtools view -F1024 $i.ccs.bam | awk '{print ">"$1"\n"$10}' > $i.fa
        blastn -query barcodes.fasta -subject $i.fa -outfmt 6 | cut -f2 | sort -u > $i.chimeras
        rm $i.fa
        samtools view -F1024 $i.ccs.bam | grep -v -f $i.chimeras | samtools view -Sbh > $i.ccs.nochim.bam
        samtools view $i.ccs.nochim.bam | awk '{print ">"$1"\n"$10}' > $i.ccs.nochim.fa
done
