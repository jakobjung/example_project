convert ess. genes to embl:

```

```



This was done to transfer the annotation from k12 to upec:

```
$RATT_HOME/start.ratt.sh ./embl_annotation_k12 ../reference_sequences/ecoli536.fasta reannotated_upec Strain U00096.fasta
```

And convert embl to gff file:

```
seqret sformat embl -sequence CP000247.1.embl -feature -osformat fasta -offormat gff -auto 
```

