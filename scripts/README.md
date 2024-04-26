# scripts

Here, all scripts needed to produce the results and plots are stored. Below, I go through the scripts in the order 
they should be run to reproduce the results. 

## 1. Run MASON command-line tool
To screen for off-targets and get PNA-specific features, I run the MASON command-line-tool 
([link to paper](https://rnajournal.cshlp.org/content/early/2023/02/07/rna.079263.122.abstract)). For that, you can get
the webserver as described in the GitHub repository [here](https://github.com/BarquistLab/mason_commandline). 
To run it on the PNAs and search for off-targets in our strain, I used the following command to move to the respective
directory and run the script:

```bash
# go to the directory where the mason.sh script is stored
cd ~/Documents/mason_commandline
# run the script
sh mason.shsh  -f ~/Documents/UPEC_K12_essgenes_2023_03/data/reference_sequences/e_coli_K12.fasta \
-g ~/Documents/UPEC_K12_essgenes_2023_03/data/reference_sequences/e_coli_K12.gff3 -m 2 -i ESSSCREEN -\
p ~/Documents/UPEC_K12_essgenes_2023_03/data/reference_sequences/pna_sequences.fasta  
```

The output of this script is a folder called `ESSSCREEN` in the `data` folder. The result of the off-target screen and
the PNA-specific features are then copied to my `./data` folder:
    

```bash
# copy result files to my data directory:
scp -r ./data/ESSSCREEN/outputs/result_table.csv ~/Documents/UPEC_K12_essgenes_2023_03/data/MASON_output/
scp -r ./data/ESSSCREEN/outputs/offtargets_fulltranscripts_sorted.tab \
~/Documents/UPEC_K12_essgenes_2023_03/data/MASON_output/
# go back to the scripts directory
cd ../UPEC_K12_essgenes_2023_03/scripts/
```

## 2. Run `pna_specific_predictors.R`
This script takes the output of the MASON command-line tool and calculates the PNA-specific features. It also
calculates the off-targets and the off-target score. The script is run as follows:

```bash
Rscript pna_specific_predictors_1.R
```

The output of this script is a file called `pna_specific_predictors.csv` in the `data` folder.

## 3. Run `gene_specific_predictors.R`

This script takes the output of the MASON command-line tool and calculates the gene-specific features. The script is
run as follows:

```bash
Rscript gene_specific_predictors_2.R
```

The output of this script is a file called `gene_specific_predictors.csv` in the `data` folder.



