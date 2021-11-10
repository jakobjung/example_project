#!/bin/bash

# map reads:
~/bin/bbmap/bbmap.sh ref=../data/reference_sequences/genome.fasta\
		     in=raw_reads_R1.fastq.gz \
		     ambig=toss \
		     out=../data/aligned_reads_R1.bam

#get count table:
featureCounts in=../data/aligned_reads_R1.bam \
	      -a ../data/reference_sequences/annotation.gff \
	      out=../analysis/count_table.tab 
