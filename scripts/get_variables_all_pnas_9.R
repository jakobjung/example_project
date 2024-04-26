# this script takes all the designed PNAs and extracts all the important variables: Tm, purine_percentage, sc_bases,
# CAI, MFE_UPEC, gene_vulnerability_crispri, string_interactions, crispri_log2FC_rousset, cas13_log2FC, expression_upec,
# PNA_molecular_weight,
# upec_total_off_targets_0mm, upec_total_off_targets_1mm, upec_total_off_targets_2mm, upec_tir_off_targets_0mm
# upec_tir_off_targets_1mm, upec_tir_off_targets_2mm

# Load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(readxl)

# Load data
pnas_0mm <- read.csv("../mason_commandline/data/all_genes_0mm/final_table.csv", header = TRUE)
pnas_1mm <- read.csv("../mason_commandline/data/all_genes_1mm/final_table.csv", header = TRUE)
pnas_2mm <- read.csv("../mason_commandline/data/all_genes_2mm/final_table.csv", header = TRUE)

# add column names to pnas_0mm: ASO	ASO_seq	target_seq	location	SC_bases	pur_perc	long_pur_stretch
# Tm	OT_tot	OT_TIR
colnames(pnas_0mm) <- c("gene", "PNA", "PNA_sequence", "target_seq", "location", "sc_bases", "purine_percentage",
                        "long_pur_stretch", "Tm", "upec_total_off_targets_0mm", "upec_tir_off_targets_0mm")

# add columns of pnas_1mm and pnas_2mm
all_pnas_all_egenes_upec <- pnas_0mm %>% mutate(upec_total_off_targets_1mm = pnas_1mm$OT_tot,
                                                upec_tir_off_targets_1mm = pnas_1mm$OT_TIR,
                                                upec_total_off_targets_2mm = pnas_2mm$OT_tot,
                                                upec_tir_off_targets_2mm = pnas_2mm$OT_TIR) %>% as_tibble()


sum(is.na(all_pnas_all_egenes_upec))
# add gene names and k12 mappings:
blastn <- read_delim("./data/reference_sequences/blast_egenes.txt", delim = "\t", col_names = FALSE)
proteinortho <- read_delim("./data/reference_sequences/upec_egenes.proteinortho.tsv", delim = "\t")

mappings_egenes <- blastn %>% select(1,2) %>% mutate(K12_locus_tag = gsub("^([^_]+)_.*", "\\1", X1)) %>%
  mutate(K12_genename = gsub("^[^_]+_([^:]+):.*", "\\1", X1)) %>%
  mutate(upec_locus_tag = gsub("^([^:]+):.*", "\\1", X2)) %>% select(-X1, -X2)


all_pnas_all_egenes_upec$K12_locus_tag <- mappings_egenes$K12_locus_tag[match(all_pnas_all_egenes_upec$gene, mappings_egenes$upec_locus_tag)]
all_pnas_all_egenes_upec$K12_genename <- mappings_egenes$K12_genename[match(all_pnas_all_egenes_upec$gene, mappings_egenes$upec_locus_tag)]


# Define a function to calculate the molecular weight of an n-mer of nucleobases connected with a peptide bond
calculate_pna_mw <- function(seq) {
  # Define the molecular weight of each nucleotide (nucleobase only)
  mw_nucleotide <- c(A=135.13, T=125.06, G=152.12, C=111.07)
  # Define the molecular weight of the peptide bond
  mw_peptidebond <- 113.12
  # Split the DNA sequence into individual nucleotides
  nucleotides <- strsplit(seq, "")[[1]]
  # Calculate the molecular weight of the sequence, accounting for the peptide bond
  mw <- sum(mw_nucleotide[nucleotides]) + mw_peptidebond * (length(nucleotides) - 1)
  return(mw)
}

all_pnas_all_egenes_upec$PNA_molecular_weight <- sapply(all_pnas_all_egenes_upec$PNA_sequence, calculate_pna_mw)

# check CAI for essential genes and check whether theres a correlation between cai and expression_upec
CAI_data <- read_delim("./data/CAI_calculation/UPEC_CAI.tsv", delim="\t", col_names = F)
all_pnas_all_egenes_upec$CAI <- CAI_data$X2[match(all_pnas_all_egenes_upec$gene, CAI_data$X1)]


gff_upec <- read.delim("./data/reference_sequences/ecoli536_sRNAs_modified.gff3", skip = 3, header = FALSE)
gff_upec$locus_tag <- gsub(".*locus_tag=([^;]+).*", "\\1", gff_upec$V9)

gff_upec_ess <- gff_upec[gff_upec$locus_tag %in% unique(all_pnas_all_egenes_upec$gene),] %>%
  filter(V3 == "CDS") %>% select(V1, V2, locus_tag, V4, V5, V6, V7, V8, V9)

gff_upec_ess_startregs <- gff_upec_ess %>%
  # change start and end to start-30 and start+30 respectively for V7=+ and end-30 and end+30 for V7=-
    mutate(V4 = ifelse(V7 == "+", V4-30, V5-30)) %>%
    mutate(V5 = ifelse(V7 == "+", V4+60, V5+30))

# save startregs to file for use in bedtools
write.table(gff_upec_ess_startregs, "./data/sec_structure_rnafold/upec_startsites.gff", sep = "\t",
            quote = F, row.names = F, col.names = F)

# pick up the MFE from ./data/sec_structure_rnafold/upec_sec_structure.tsv
sec_structure_upec <- read.table("./data/sec_structure_rnafold/RNAfold_tab_upec.tsv", header = F, sep = "\t")

all_pnas_all_egenes_upec$MFE_UPEC <- sec_structure_upec$V4[match(all_pnas_all_egenes_upec$gene, sec_structure_upec$V1)]

# add gene_vulnerability_crispri
crispr_fitness <- read.table("./data/gene_essentiality/rel_fitness_ecoli_cellsystems.csv", header = T, sep = "\t")
crispr_fitness <- crispr_fitness %>%
  filter(!is.na(relative.fitness..mean.)) %>%
  select(locus_tag,gene, relative.fitness..mean.) %>%
  # summarize locus tag and gene name with relative.fitness..mean.
    group_by(locus_tag, gene) %>%
    summarise(RF_mean_cellsytems = mean(relative.fitness..mean.)) %>%
  mutate(gene_vulnerability_crispri = RF_mean_cellsytems * (-1) + 1)

all_pnas_all_egenes_upec$gene_vulnerability_crispri <- crispr_fitness$gene_vulnerability_crispri[match(all_pnas_all_egenes_upec$K12_locus_tag, crispr_fitness$locus_tag)]
# remove NAs and make them mean
all_pnas_all_egenes_upec$gene_vulnerability_crispri[is.na(all_pnas_all_egenes_upec$gene_vulnerability_crispri)] <- mean(all_pnas_all_egenes_upec$gene_vulnerability_crispri, na.rm = T)

# import crispri data from rousset 2018
crispr_fitness_rousset <- read_delim("./data/gene_essentiality/plos_gen_rousset_2018.csv", col_names = T, delim = ",") %>%
  select(gene, median_coding)

all_pnas_all_egenes_upec$crispri_log2FC_rousset <- crispr_fitness_rousset$median_coding[match(all_pnas_all_egenes_upec$K12_genename, crispr_fitness_rousset$gene)] * (-1)
# remove NAs and make them mean
all_pnas_all_egenes_upec$crispri_log2FC_rousset[is.na(all_pnas_all_egenes_upec$crispri_log2FC_rousset)] <- mean(all_pnas_all_egenes_upec$crispri_log2FC_rousset, na.rm = T)

# import STRING data
string_data <- read_delim("./data/protein_protein_interactions/stringdb_k12.txt", col_names = T, delim = " ")
string_high_confidence_lt <- string_data[string_data$combined_score>=700,] %>% select(protein1, protein2) %>%
  pivot_longer(cols = c(protein1, protein2), names_to = "protein", values_to = "locus_tag") %>%
  mutate(locus_tag = gsub(".*\\.", "", locus_tag)) %>% select(-protein) %>%
  unlist() %>% as.character
# for all lts, get the count of occurences in the string_high_confidence_lt vector
string_high_confidence_lt_count <- table(string_high_confidence_lt)

all_pnas_all_egenes_upec$string_interactions <- string_high_confidence_lt_count[match(all_pnas_all_egenes_upec$K12_locus_tag, names(string_high_confidence_lt_count))]
# make 0 for NAs
all_pnas_all_egenes_upec$string_interactions[is.na(all_pnas_all_egenes_upec$string_interactions)] <- 0



#add cas13 data
cas13_gnames <- read_excel("./data/gene_essentiality/biorxiv_doudna_2023.xlsx", sheet = "Table S12") %>%
  select("oligo_name", "gene")

cas13_data <- read_excel("./data/gene_essentiality/biorxiv_doudna_2023.xlsx", sheet = "Table S16") %>%
  select("oligo_name", "dCas13_LB_S2_normalized_vs_T0_log2FC") %>%
  # add gene name
    left_join(cas13_gnames, by = "oligo_name") %>%
  # get mean of log2FC for each gene
    group_by(gene) %>%
    summarise(mean_log2FC = mean(dCas13_LB_S2_normalized_vs_T0_log2FC))

all_pnas_all_egenes_upec$cas13_log2FC <- cas13_data$mean_log2FC[match(all_pnas_all_egenes_upec$K12_genename, cas13_data$gene)] * (-1)

# read in the expression data
tpm_upec <- read.table("./data/transcriptomic_expression_K12_upec/upec_ctrl_log_tpm_plus1.csv", header = T,
                       sep = ",", row.names = 1)

all_pnas_all_egenes_upec$expression_upec <- rowMeans(tpm_upec[match(all_pnas_all_egenes_upec$gene, rownames(tpm_upec)),])

# save the data as a tsv
write.table(all_pnas_all_egenes_upec, "./data/all_genes_all_pnas/all_pnas_all_egenes_upec.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

data_for_ml <- all_pnas_all_egenes_upec[c('upec_tir_off_targets_2mm',
                     'purine_percentage',
                     'Tm',
                     'upec_tir_off_targets_1mm',
                     'sc_bases',
                     'gene_vulnerability_crispri',
                     'string_interactions',
                     'crispri_log2FC_rousset',
                     'cas13_log2FC',
                     'PNA_molecular_weight',
                     'expression_upec',
                     'MFE_UPEC',
                     'upec_total_off_targets_0mm',
                     'upec_total_off_targets_1mm',
                     'upec_total_off_targets_2mm',
                     'upec_tir_off_targets_0mm')]
# check data for ml whether there are any NAs
sum(is.na(data_for_ml))
# find the rows with NAs
sum(is.na(data_for_ml[complete.cases(data_for_ml),]))





