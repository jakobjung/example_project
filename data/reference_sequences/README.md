# reference sequences

## sRNA annotation

I blasted the sRNAs from http://srd.genouest.org/ (14-06-2023), and blasted them against my genome: 

```bash
blastn -query sRNA_N315.fasta -subject Staphylococcus_aureus_101588.fasta -max_target_seqs 1 -out blast_srnas.txt -outfmt 6
```

Then I checked for the srna lengths:

```bash
bioawk -c fastx '{ print $name, length($seq) }' < sRNA_N315.fasta > lengths.tsv
```

Then I ran the R-script to merge srna annotations.

```bash
R add_sRNAs_to_gff.R
```

