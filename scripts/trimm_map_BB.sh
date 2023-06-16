#!/bin/bash
main(){
    # go to project directory where the reads and reference
    # sequences are stored:
    PROJECT=../data
    echo "Start trimming"
    rename_trim_rna_libs
    echo "Trimming done. Start mapping"
    align_rna_reads_genome
    echo "Finished mapping. Start connecting all tab files"
    featureCounts -T 5 -t CDS,sRNA,tRNA -g Name \
		  -a $PROJECT/reference_sequences/WG_with_sRNAs_Staphylococcus_aureus_101588.gff3\
		  -o $PROJECT/rna_align/counttable.txt \
		  $PROJECT/rna_align/*.bam
    
}

rename_trim_rna_libs(){
    mkdir -pv $PROJECT/libs
    # shellcheck disable=SC2045
    for NAME in $(ls $PROJECT/fastq/*.fq.gz)
    do
        echo "$NAME starts trimming nowwwwwww"
        NEWNAME=${NAME##*/}
        NEWNAME=${NEWNAME%.fastq.gz}_trimmed.fastq.gz
        echo $NEWNAME
       # bbduk trims low quality bases and removes adapters:
        bbduk.sh  in=$NAME \
			     ref=../data/reference_sequences/adapters.fa -Xmx4g t=20\
			     out=$PROJECT/libs/${NEWNAME} ktrim=r k=23 mink=11\
			     hdist=1 qtrim=r trimq=10 ftl=12 
    done
}


align_rna_reads_genome(){
    mkdir -p $PROJECT/rna_align
    DIR=$PROJECT/rna_align
    for i in $(ls $PROJECT/libs/*.fastq.gz)
    do
        NAME=${i##*/}
        NAME=${NAME%_trimmed.fastq.gz}
        echo "Starting mapping for sample: $NAME"
        bbmap.sh in=$i trimreaddescription=t  t=20 \
			     ref=$PROJECT/reference_sequences/Staphylococcus_aureus_101588.fasta \
			     k=12 outm=$DIR/$NAME.sam
        # sort sam file, create BAM file:
        samtools sort -O BAM -@ 40 $DIR/$NAME.sam > $DIR/$NAME.bam
        # remove sam file: (not actually needed)
        rm $DIR/$NAME.sam
        # index BAM file:
        samtools index $DIR/$NAME.bam
        # generate coverage statistics with bamcoverage:
        bamCoverage -b $DIR/$NAME.bam -o $DIR/$NAME.bw -bs 5
	      fastqc $i
	      rm $i
    done
}

main  
