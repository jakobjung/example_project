# this script goes through all genes (locus tags are stored in ./data/all_genes_all_pnas/upec_essgenes.csv" as tsv file)
# and runs mason on each gene.

# activate conda environment
#source activate /home/jj/miniconda3/envs/MC_new

# walk through each row of the file and run mason on each gene
FILE="../data/all_genes_all_pnas/upec_essgenes.csv"
echo "wut"

# run for loop (only one entry(id) in each row)
while IFS=, read -r gene_id
do
  echo "Running mason on gene: $gene_id"
  mason.sh -f "../data/reference_sequences/ecoli536.fasta" -g \
  "../data/reference_sequences/ecoli536_sRNAs_modified.gff3" -m 2 -i "all_genes/$gene_id" -t "$gene_id" -l 9
  echo "done with gene: $gene_id"
done < $FILE

echo "done with all genes"





