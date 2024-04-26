
## code

From R code (script 4) we get upecgenes.gff3.  Now we can get the fasta sequences of these genes:

```bash
bedtools getfasta -fi ecoli536.fasta -bed upecgenes.gff3 -name+ -fo genes_upec.fasta -s
```



Show linux bash command to perform blast of all k12 essential genes, and check whether they exist in upec:

```bash
blastn -query egenes_k12.fasta -subject ecoli536.fasta -max_target_seqs 1 -out tot_blast_egenes.txt -outfmt 6  -max_hsps 1 -word_size 8
```



i do a similar thing with proteinortho:

```bash
proteinortho6.pl --project=upec_egenes ./egenes_k12.fasta ./genes_upec.fasta --p=blastn --e=500 --identity=5 --cov=5
```

