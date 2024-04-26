Here, RNAfold is used to calculate mRNA folding at target genes' TIR.

## K12

Includes reference sequences as well as anotation files.

```bash
RNAfold -i tirs_K12.fasta | perl -pe 's/\n/\t/g; s/>//; s/\s+/\t/; s/\(-|\( -/-/; s/\)\t$/\n/' | sed 's/gene-//' >> sec_structure.tsv 
```

Now we can use the data to check sec. structures and MFE.

## UPEC

for upec, we use the generated (in get_essgenes_upec.R) GFF file for all start regions of ess. genes. 



better (new:)

```bash
 RNAfold -i essgenes_ss_upec.fasta > rnafold_output.txt
 perl -pe 's/\n/\t/g; s/>//; s/\s+/\t/; s/\(-|\( -/-/; s/\)\t$/\n/' < rnafold_output.txt  > RNAfold_tab.tsv 
```







```bash
# extract fasta
bedtools getfasta -fi ../reference_sequences/ecoli536.fasta -bed upec_startsites.gff  -name -fo ess_genes_startsites_upec.fasta -s
sed "s/(+)//" ess_genes_startsites_upec.fasta | sed "s/(-)//"  > essgenes_ss_upec.fasta
rm -rf ess_genes_startsites_upec.fasta


# get rnafold:
RNAfold -i essgenes_ss_upec.fasta | perl -pe 'if (/^>/) { s/>//; print "\n$_\t"; } else { chomp; $seq .= $_; } END { print "$seq\n" if $seq; }' | sed 's/gene-//' >> upec_sec_structure.tsv


```

