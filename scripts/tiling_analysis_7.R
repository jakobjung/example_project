# script generates data for the tiling dataset.
#  get libraries:
library(tidyverse)
library(readxl)

# import files:
# pna_info from excel file
pna_info <- read_excel("./data/tiling_test/pnas_giorgia.xlsx") %>% filter(nchar(Name)<15) %>%
  mutate(Number = gsub("-", "_", Number))
# import mic data from csv:
mic_data <- read_csv("./data/tiling_test/tiling_mics.csv") %>% filter(grepl("Jv", PNA))
# import all predictors data from csv:
all_predictors <- read_delim("./data/pnas_predictors_mic_upec.tsv")

# keep data for only 1 gene_name for all duplicate gene_names in in all_predictors
all_predictors <- all_predictors %>% filter(!duplicated(gene_name))

# add info to mic data, use PNA as key in mic_data and Number in pna_info
mic_data_mod <- left_join(mic_data, pna_info, by = c("PNA" = "Number")) %>%
  select(PNA, MIC, Name, "PNA Sequence", concentration)

colnames(mic_data_mod) <- c("pna_name", "MIC", "pna_name_full", "pna_sequence", "concentration")

# save as fasta file with pna_name as header and pna_sequence as sequence
# create fasta file:
library(phylotools)
fasta_df <- data.frame(mic_data_mod[,c(1,4)])
colnames(fasta_df) <- c("seq.name", "seq.text")
dat2fasta(fasta_df, outfile = "./data/tiling_test/tiling_PNAs.fasta")


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

# read mason outputs
restab_mason_0mm <- read_csv("./data/tiling_test/result_table_0mm.csv")
restab_mason_1mm <- read_csv("./data/tiling_test/result_table_1mm.csv")  %>%
  mutate(OT_TIR = OT_TIR-restab_mason_0mm$OT_TIR, OT_tot = OT_tot - restab_mason_0mm$OT_tot)
restab_mason_2mm <- read_csv("./data/tiling_test/result_table_2mm.csv") %>%
  mutate(OT_TIR = OT_TIR-restab_mason_0mm$OT_TIR - restab_mason_1mm$OT_TIR, OT_tot = OT_tot - restab_mason_0mm$OT_tot -
  restab_mason_1mm$OT_tot)


# get target gene names from all_predictors and add gene specific features to mic_data_mod
mic_data <- mic_data_mod %>%
  # extract target gene name
  mutate(gene_name = gsub("KFF-([^_]+)_\\d.*", "\\1", pna_name_full)) %>%
  # get only PNAs that target ATG
    filter(grepl("CAT", pna_sequence)) %>%
  # get only PNAs of length 9
    filter(nchar(pna_sequence) == 9) %>%
  # add purine percentage
    mutate(purine_percentage = ((str_count(pna_sequence, "A") + str_count(pna_sequence, "G"))/9)*100) %>%
  # add gene specific features from all_predictors: expression, MFE, relative_fitness
    left_join(all_predictors, by = "gene_name") %>%
  # select only the columns we need
    select(pna_name.x, pna_name_full, pna_sequence.x, concentration, gene_name, expression_upec, MFE_UPEC,
           crispri_log2FC_rousset, string_interactions,cas13_log2FC,
           gene_vulnerability_crispri, MIC, purine_percentage.x, CAI) %>%
  # remove the .x from the column names
    rename_with(.cols = c(1,3,7,13), ~str_replace(.x, "\\.x", "")) %>%
  # calculate molecular weight
    mutate(PNA_molecular_weight = sapply(pna_sequence, calculate_pna_mw)) %>%
  # add SC_bases from mason output (restab_mason_0mm) by pna_name/ASO
    left_join(restab_mason_0mm, by = c("pna_name" = "ASO")) %>%
  # keep only the columns we need
    select(pna_name, pna_name_full, pna_sequence, concentration, gene_name, expression_upec, MFE_UPEC,
           crispri_log2FC_rousset, string_interactions,cas13_log2FC,
           gene_vulnerability_crispri, MIC, purine_percentage, PNA_molecular_weight, SC_bases, CAI, OT_TIR,OT_tot) %>%
  # rename OT_TIR to tir_off_targets_0mm
    rename(upec_tir_off_targets_0mm = OT_TIR) %>%
  # rename OT_tot to upec_total_off_targets_0mm
    rename(upec_total_off_targets_0mm = OT_tot) %>%
  # add tir offtargets_1mm
    left_join(restab_mason_1mm, by = c("pna_name" = "ASO")) %>%
    # keep only the columns we need
    select(pna_name, pna_name_full, pna_sequence, concentration, gene_name, expression_upec, MFE_UPEC,
           gene_vulnerability_crispri, MIC, purine_percentage, PNA_molecular_weight, CAI, upec_tir_off_targets_0mm, upec_total_off_targets_0mm,
           crispri_log2FC_rousset, string_interactions,cas13_log2FC,
           OT_TIR, OT_tot) %>%
    # rename OT_TIR to upec_tir_off_targets_1mm
    rename(upec_tir_off_targets_1mm = OT_TIR) %>%
    # rename OT_tot to upec_total_off_targets_1mm
    rename(upec_total_off_targets_1mm = OT_tot) %>%
    # add total offtargets_2mm
    left_join(restab_mason_2mm, by = c("pna_name" = "ASO")) %>%
    # keep only the columns we need
    select(pna_name, pna_name_full, pna_sequence, concentration, gene_name, SC_bases, expression_upec, MFE_UPEC,Tm,
           gene_vulnerability_crispri, MIC, purine_percentage, PNA_molecular_weight, CAI, upec_tir_off_targets_0mm, upec_total_off_targets_0mm,
           crispri_log2FC_rousset, string_interactions,cas13_log2FC,
              upec_tir_off_targets_1mm, upec_total_off_targets_1mm, OT_TIR, OT_tot) %>%
    # rename OT_tot to upec_total_off_targets_2mm
    rename(upec_total_off_targets_2mm = OT_tot) %>%
    # rename OT_TIR to upec_tir_off_targets_2mm
    rename(upec_tir_off_targets_2mm = OT_TIR) %>%
  # rename SC_bases to sc_bases
    rename(sc_bases = SC_bases)

# save tiling data as tsv
write_tsv(mic_data, "./data/tiling_test/tiling_data.tsv")

# create barplot of MICs (5,10,20,>20)
mic_data %>%
  # add MIC category
  mutate(MIC_category = ifelse(MIC == 5, "MIC = 5", ifelse(MIC == 10, "MIC = 10", ifelse(MIC == 20, "MIC = 20",
    "MIC > 20")))) %>%
  # count number of PNAs in each MIC category
  count(MIC_category) %>%
  # plot barplot
  ggplot(aes(x = MIC_category, y = n)) +
  geom_bar(stat = "identity") +
  labs(x = "MIC category", y = "Number of PNAs", title = "Number of PNAs in each MIC category")



