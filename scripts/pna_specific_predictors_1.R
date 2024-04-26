##
##   Title:       PNA- essential genes screen - clean table and calculate pna-specific predictors
##
##   Author:      Jakob Jung
##
##   Date:        2023-03-08
##
##   Description: Clean the excel table of the PNAs that has been used for the essential genes screen and calculate
##   pna-specific factors such as melting temperature (Tm), start position of the PNA sequence, off-targets, etc.
##   Also, the script generates several plots to visualize the data.
##

# Load libraries
library(tidyverse)
library(rmelting)
library(stringr)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(eulerr)

# Load data
pnas <- read_excel("./data/PNAscreening_Ecoli_MIC.xlsx")

restab_mason_0mm <- read_csv("./data/MASON_output/result_table_0mm.csv") %>% select(OT_tot, OT_TIR)
restab_mason_1mm <- read_csv("./data/MASON_output/result_table_1mm.csv")  %>%
  mutate(OT_TIR = OT_TIR-restab_mason_0mm$OT_TIR, OT_tot = OT_tot - restab_mason_0mm$OT_tot) %>%
  select(OT_tot, OT_TIR)
restab_mason_2mm <- read_csv("./data/MASON_output/result_table.csv") %>%
  mutate(OT_TIR = OT_TIR-restab_mason_0mm$OT_TIR - restab_mason_1mm$OT_TIR, OT_tot = OT_tot - restab_mason_0mm$OT_tot -
  restab_mason_1mm$OT_tot)


# Clean data
pnas_filtered <- pnas %>%
  select(-c(3,4,6,8,9,10,11,12,13,15)) %>%
  # rename all columns with custom names
  rename(
        "pna_sequence" = "PNA sequence",
        "e_coli_k12_inhibition" = "is_inhibiting(K-12 (2/2))",
        "upec_inhibition" = "is_inhibiting(UPEC (2/2))",
        "pna_name" = "Number Jvpna",
        "gene_name_locus_tag" = "Gene Name/locus tag"
  ) %>%
  # change all - to _ in pna_name
  mutate(pna_name = str_replace_all(pna_name, "-", "_")) %>%
  # for each NA in gene_name_locus_tag, replace with the value of the previous row
  mutate(gene_name_locus_tag = ifelse(is.na(gene_name_locus_tag), lag(gene_name_locus_tag), gene_name_locus_tag)) %>%
  # separate locus tag (in parenthesis and at the end of string) from gene names (no parentheses and at the start)
  mutate(locus_tag = gsub(pattern = "^[a-zA-Z]+ .*\\((b\\d+)\\)", replacement = "\\1", gene_name_locus_tag)) %>%
  # get gene name for each gene in similar way:
  mutate(gene_name = gsub(pattern = "^([a-zA-Z]+).*", replacement = "\\1", gene_name_locus_tag)) %>%
  # split Position column into start_position and end_position separated by one or more spaces
  separate(Position, into = c("start_position", "end_position"), sep = " +") %>%
  # turn start_position and end_position into numeric values and e_coli_k12_inhibition and upec_inhibition into logical
  # values
  mutate(e_coli_k12_inhibition = as.logical(e_coli_k12_inhibition), upec_inhibition = as.logical(upec_inhibition),
         start_position = as.numeric(start_position), end_position = as.numeric(end_position)) %>%
  # subtract one to all positive values of start_position and end_position
  mutate(start_position = ifelse(start_position > 0, start_position - 1, start_position),
           end_position = ifelse(end_position > 0, end_position - 1, end_position)) %>%
  # add total off-targets, TIR off-targets, SC-bases, purine_percentage, longest_purine_stretch, and Tm columns
  mutate(target_seq = restab_mason_2mm$target_seq, total_off_targets_2mm = restab_mason_2mm$OT_tot,
         tir_off_targets_2mm = restab_mason_2mm$OT_TIR, total_off_targets_1mm = restab_mason_1mm$OT_tot,
         tir_off_targets_1mm = restab_mason_1mm$OT_TIR, total_off_targets_0mm = restab_mason_0mm$OT_tot,
         tir_off_targets_0mm = restab_mason_0mm$OT_TIR, sc_bases = restab_mason_2mm$SC_bases,
         purine_percentage = restab_mason_2mm$pur_perc, longest_purine_stretch = restab_mason_2mm$long_pur_stretch,
         Tm = restab_mason_2mm$Tm) %>%
  # create new column called "PNA_GC_content" and calculate GC content of each PNA sequence
  mutate(PNA_GC_content = str_count(pna_sequence, "[GC]") / nchar(pna_sequence)) %>%
  # create new column called "homopolymers" and calculate the Length of longest consecutive nucleotides in each PNA
  # sequence
  mutate(homopolymers = sapply(str_extract_all(pna_sequence, "([ATGC])\\1+"),
                               function (x) ifelse(length(nchar(x)) == 0, 1,  max(nchar(x))))) %>%
  # create column showing whether either e_coli_k12_inhibition or upec_inhibition is TRUE
  mutate(inhibits_either = e_coli_k12_inhibition | upec_inhibition) %>%
  # create new column called "PNA_length" and calculate the length of each PNA sequence
  select(c("pna_name", "gene_name", "locus_tag", "pna_sequence", "target_seq",  "start_position",
           "total_off_targets_2mm", "tir_off_targets_2mm", "total_off_targets_1mm", "tir_off_targets_1mm",
           "total_off_targets_0mm", "tir_off_targets_0mm","sc_bases","purine_percentage", "longest_purine_stretch", "Tm",
           "PNA_GC_content", "homopolymers", "e_coli_k12_inhibition", "upec_inhibition", "inhibits_either")) %>%
  # create new column called "inhibition_of" and assign a value to each PNA based on whether it inhibits K-12, UPEC, or
    # both. Order of levels is important for plotting later
    mutate(inhibition_of = as.factor ( ifelse(e_coli_k12_inhibition & upec_inhibition, "both",
                                  ifelse(upec_inhibition, "UPEC only",
                                         ifelse(e_coli_k12_inhibition, "K12 only", "neither"))))) %>%
  # order the factor levels of inhibition_of
    mutate(inhibition_of = factor(inhibition_of, levels = c("neither","K12 only","UPEC only","both"))) %>%
    # create new column called "OT_TIR_2mm/OT_2mm" and calculate the ratio of TIR off-targets to total off-targets for
    # each PNA sequence
    mutate(OT_TIR_2mm_OT_2mm = tir_off_targets_2mm / total_off_targets_2mm)

# View(pnas_filtered)

# create one-hot encoding of PNAs:upec_inhibition
self_encode <- function(sequence) {
  integer_encoded <- matrix(0, nrow = nchar(sequence), ncol = 4)
  nts <- c('A', 'T', 'C', 'G')
  nt_names <- rep(nts, each = nchar(sequence))
  pos_names <- rep(1:nchar(sequence), 4)
  nt_pos_names <- paste0(nt_names, pos_names)
  for (i in 1:nchar(sequence)) {
    integer_encoded[i, match(substr(sequence, i, i), nts)] <- 1
  }
  sequence_one_hot_encoded <- as.vector(integer_encoded)
  names(sequence_one_hot_encoded) <- nt_pos_names
  return(sequence_one_hot_encoded)
}

# create one-hot encoding of PNA sequences
one_hot_encodings <- as_tibble(t(sapply(pnas_filtered$pna_sequence, self_encode)))

# add one-hot encoding columns to pnas_filtered
pnas_filtered <- cbind(pnas_filtered, one_hot_encodings)

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

# calculate molecular weight of each PNA sequence
pnas_filtered <- pnas_filtered %>%
  mutate(PNA_molecular_weight = sapply(pna_sequence, calculate_pna_mw))

# save data to tsv file
write_tsv(pnas_filtered, "./data/pnas_pna_specific_features.tsv")


## Plots:

# create an euler diagram showing the number of PNAs inhibiting K12 only, UPEC only, both, or neither with the
# euler function from the eulerr package
# get nr of PNAs only inhibiting UPEC
nr_pnas_only_upec <- nrow(pnas_filtered[pnas_filtered$upec_inhibition & !(pnas_filtered$e_coli_k12_inhibition),])
# get nr of PNAs only inhibiting K12
nr_pnas_only_k12 <- nrow(pnas_filtered[!(pnas_filtered$upec_inhibition) & pnas_filtered$e_coli_k12_inhibition,])
# get nr of PNAs inhibiting both
nr_pnas_both <- nrow(pnas_filtered[pnas_filtered$upec_inhibition & pnas_filtered$e_coli_k12_inhibition,])
# get nr of PNAs inhibiting none
nr_pnas_none <- nrow(pnas_filtered[!(pnas_filtered$upec_inhibition) & !(pnas_filtered$e_coli_k12_inhibition),])

# create plot using euler function of eulerr package
fit <- euler(c("UPEC" = nr_pnas_only_upec, "K12" = nr_pnas_only_k12, "UPEC&K12" = nr_pnas_both, "None" = nr_pnas_none))
svg(filename = "./analysis/gene_specific_predictors/euler_plot_normal.svg", width = 8, height = 6)
plot(fit, main = "PNA-mediated growth inhibition at 10uM", quantities = TRUE,
     fills = list(fill = c("#b5cde1", "#66cdaa", "#d2d2d2", "#4682b4"),
                  alpha=1))
dev.off()

svg(filename = "./analysis/gene_specific_predictors/euler_plot_mod.svg", width = 8, height = 6)
plot(fit, main = "PNA-mediated growth inhibition at 10uM", quantities = TRUE,
     fills = list(fill = c("powderblue", "powderblue", "#d2d2d2", "powderblue"),
                  alpha=1))
dev.off()



# generate plot showing start position as x axis and number PNAs inhibiting K12 only, UPEC only, both, or neither as
# y axis in a stacked bar plot (use inhibition_of column). Add a legend with fill colors
# generate a function that creates the plot same as above, but being able to change the x axis
startpos_plot_function <- function(x_axis) {
    startpos_plot <- pnas_filtered %>%
        ggplot(aes(x = !!sym(x_axis), fill = inhibition_of)) +
        geom_bar() +
      # order the legend by the order of the levels in the factor
        scale_fill_manual(values = c("K12 only" = "#66cdaa", "UPEC only" = "#b5cde1", "both" = "#4682b4",
                                     "neither" = "#d2d2d2")) +
        theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13)
    ) +
        theme(axis.text.x = element_text( hjust = 0.7, vjust = 1)) +
        labs(x = x_axis, y = "Number of PNAs", fill = "Inhibition of")
    return(startpos_plot)
}


startpos_plot <- startpos_plot_function("start_position")
startpos_plot
# save plot as svg file
ggsave( filename = "./analysis/pna_specific_predictors/startpos_plot.svg", plot = startpos_plot, width = 10, height = 5)


# generate plots for all possible x axes
plot_gc <- startpos_plot_function("PNA_GC_content")
ggsave( filename = "./analysis/pna_specific_predictors/PNA_GC_content.svg", plot = plot_gc, width = 10, height = 5)
plot_homo <- startpos_plot_function("homopolymers")
plot_homo
ggsave( filename = "./analysis/pna_specific_predictors/PNA_max_homopolymer_length.svg", plot = plot_homo, width = 10,
        height = 5)
plot_longest_purine <- startpos_plot_function("longest_purine_stretch")
plot_longest_purine
ggsave( filename = "./analysis/pna_specific_predictors/PNA_longest_purine_stretch.svg", plot = plot_longest_purine,
        width = 10, height = 5)
plot_purine_perc <- startpos_plot_function("purine_percentage")
plot_purine_perc
ggsave( filename = "./analysis/pna_specific_predictors/PNA_purine_percentage.svg", plot = plot_purine_perc, width = 10,
        height = 5)

plot_sc_bases <- startpos_plot_function("sc_bases")
plot_sc_bases
ggsave( filename = "./analysis/pna_specific_predictors/PNA_sc_bases.svg", plot = plot_sc_bases, width = 10, height = 5)


# get genes that have alternate inhibits_either values for different start positions
genes_with_alternate_inhibits_either <- pnas_filtered %>%
  group_by(gene_name) %>%
  summarise(
    inhibits_either = sum(inhibits_either),
    total = n()
  ) %>%
  filter(inhibits_either != total & inhibits_either != 0) %>%
  select(gene_name)

# plot a plot showing Tm (y axis) and whether the PNA sequence is inhibiting E. coli K12 (x axis) as violon plot with
# showing single values as points.
# write function for continuous variables:
plot_function <- function(x_axis, dataset){
  Tm_plot <- dataset %>%
    ggplot(aes(x = inhibits_either, y = !!sym(x_axis))) +
    geom_violin(aes(fill=inhibits_either),width = 0.8,  alpha=1) +
    geom_boxplot(aes(fill=inhibits_either),width = 0.3,  alpha=1) +
    scale_fill_manual(values = c("TRUE" = "#4682b4", "FALSE" = "#d2d2d2"), limits = c("TRUE", "FALSE"))  +
    geom_jitter(width = 0.1, height = 0, alpha = 0.5, size = 2) +
    # ad p-values using wilcox.test and stat_compare_means and put them on top of the plot, write as "p-value = 0.0001"
    stat_compare_means(method = "wilcox.test",
                       label.y.npc = 1, size = 4,  label.x = 0.55) +
    stat_compare_means(comparisons = list(c("TRUE", "FALSE")), method = "wilcox.test",
                       label.y.npc = 1 , size = 5,  label = "p.signif") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.position = "none"
    ) +
    labs(x = "Growth-inhibition of E. coli K12 and/or UPEC", y = x_axis, fill = NULL, size=12)
  return(Tm_plot)
}

# generate plots for all possible y axes (continuous variables) and save them as svg files:
Tm_plot <- plot_function("Tm", pnas_filtered) +
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
Tm_plot
ggsave( filename = "./analysis/pna_specific_predictors/Tm_plot.svg", plot = Tm_plot, width = 7, height = 5)
Mw_plot <- plot_function("PNA_molecular_weight", pnas_filtered)+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
Mw_plot
ggsave( filename = "./analysis/pna_specific_predictors/Mw_plot.svg", plot = Mw_plot, width = 7, height = 5)
off_targets_1mm_plot <- plot_function("total_off_targets_1mm", pnas_filtered)+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
off_targets_1mm_plot
ggsave( filename = "./analysis/pna_specific_predictors/off_targets_1mm_plot.svg", plot = off_targets_1mm_plot,
        width = 7, height = 5)
off_targets_2mm_plot <- plot_function("total_off_targets_2mm", pnas_filtered)+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
off_targets_2mm_plot
ggsave( filename = "./analysis/pna_specific_predictors/off_targets_2mm_plot.svg", plot = off_targets_2mm_plot,
        width = 7, height = 5)
off_targets_0mm_plot <- plot_function("total_off_targets_0mm", pnas_filtered)+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
off_targets_0mm_plot
ggsave( filename = "./analysis/pna_specific_predictors/off_targets_0mm_plot.svg", plot = off_targets_0mm_plot,
        width = 7, height = 5)
tir_off_targets_0mm_plot <- plot_function("tir_off_targets_0mm", pnas_filtered)+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
tir_off_targets_0mm_plot
ggsave( filename = "./analysis/pna_specific_predictors/tir_off_targets_0mm_plot.svg", plot = tir_off_targets_0mm_plot,
        width = 7, height = 5)
tir_off_targets_1mm_plot <- plot_function("tir_off_targets_1mm", pnas_filtered)+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
tir_off_targets_1mm_plot
ggsave( filename = "./analysis/pna_specific_predictors/tir_off_targets_1mm_plot.svg", plot = tir_off_targets_1mm_plot,
        width = 7, height = 5)
tir_off_targets_2mm_plot <- plot_function("tir_off_targets_2mm", pnas_filtered)+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
tir_off_targets_2mm_plot
ggsave( filename = "./analysis/pna_specific_predictors/tir_off_targets_2mm_plot.svg", plot = tir_off_targets_2mm_plot,
        width = 7, height = 5)
tir_off_targets_2mm_div <- plot_function("OT_TIR_2mm_OT_2mm", pnas_filtered) +
  ylab("TIR off-targets / all off-targets")+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
tir_off_targets_2mm_div
ggsave( filename = "./analysis/pna_specific_predictors/tir_off_targets_2mm_div.svg", plot = tir_off_targets_2mm_div,
        width = 7, height = 5)

# get PNAs that inhibit both, vs PNAs that inhibit none
pnas_only_double <- pnas_filtered %>%
  # filter genes that inhibit either E. coli K12 or UPEC, but not both
  filter((upec_inhibition == "TRUE" & e_coli_k12_inhibition == "TRUE") |
           (upec_inhibition == "FALSE" & e_coli_k12_inhibition == "FALSE"))

# generate all plots above for only these PNAs
Tm_plot_only_double <- plot_function("Tm", pnas_only_double) + labs(x = "Growth-inhibition of E. coli K12 and UPEC or neither",size=12) +
  scale_x_discrete(limits = c("FALSE", "TRUE"), labels = c("FALSE" = "None", "TRUE" = "Both"))
Tm_plot_only_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/Tm_plot_only_double.svg", plot = Tm_plot_only_double, width = 7, height = 5)
Mw_plot_only_double <- plot_function("PNA_molecular_weight", pnas_only_double) + labs(x = "Growth-inhibition of E. coli K12 and UPEC or neither",size=12) +
  scale_x_discrete(limits = c("FALSE", "TRUE"), labels = c("FALSE" = "None", "TRUE" = "both"))
Mw_plot_only_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/Mw_plot_only_double.svg", plot = Mw_plot_only_double, width = 7, height = 5)
off_targets_1mm_plot_only_double <- plot_function("total_off_targets_1mm", pnas_only_double) + labs(x = "Growth-inhibition of E. coli K12 and UPEC or neither",size=12)
off_targets_1mm_plot_only_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/off_targets_1mm_plot_only_double.svg", plot = off_targets_1mm_plot_only_double,
        width = 7, height = 5)
off_targets_2mm_plot_only_double <- plot_function("total_off_targets_2mm", pnas_only_double) + labs(x = "Growth-inhibition of E. coli K12 and UPEC or neither",size=12)
off_targets_2mm_plot_only_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/off_targets_2mm_plot_only_double.svg", plot = off_targets_2mm_plot_only_double,
        width = 7, height = 5)
off_targets_0mm_plot_only_double <- plot_function("total_off_targets_0mm", pnas_only_double) + labs(x = "Growth-inhibition of E. coli K12 and UPEC or neither",size=12)
off_targets_0mm_plot_only_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/off_targets_0mm_plot_only_double.svg", plot = off_targets_0mm_plot_only_double,
        width = 7, height = 5)
tir_off_targets_0mm_plot_only_double <- plot_function("tir_off_targets_0mm", pnas_only_double) + labs(x = "Growth-inhibition of E. coli K12 and UPEC or neither",size=12)
tir_off_targets_0mm_plot_only_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/tir_off_targets_0mm_plot_only_double.svg", plot = tir_off_targets_0mm_plot_only_double,
        width = 7, height = 5)
tir_off_targets_1mm_plot_only_double <- plot_function("tir_off_targets_1mm", pnas_only_double) + labs(x = "Growth-inhibition of E. coli K12 and UPEC or neither",size=12)
tir_off_targets_1mm_plot_only_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/tir_off_targets_1mm_plot_only_double.svg", plot = tir_off_targets_1mm_plot_only_double,
        width = 7, height = 5)
tir_off_targets_2mm_plot_only_double <- plot_function("tir_off_targets_2mm", pnas_only_double) + labs(x = "Growth-inhibition of E. coli K12 and UPEC or neither",size=12)
tir_off_targets_2mm_plot_only_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/tir_off_targets_2mm_plot_only_double.svg", plot = tir_off_targets_2mm_plot_only_double,
        width = 7, height = 5)
tir_off_targets_2mm_div_only_double <- plot_function("OT_TIR_2mm_OT_2mm", pnas_only_double) +
  ylab("TIR off-targets / all off-targets")
tir_off_targets_2mm_div_only_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/tir_off_targets_2mm_div_only_double.svg", plot = tir_off_targets_2mm_div_only_double,
        width = 7, height = 5)






# plot a plot showing Tm (y axis) and whether the PNA sequence is inhibiting E. coli K12 (x axis) as violon plot with
# showing single values as points, but only for genes that have alternate inhibits_either values for different start
# positions. combine points belonging to the same gene with a connecting line
Tm_plot_only_alternating_genes <- pnas_filtered %>% #filter(sc_bases<5) %>%
  filter(gene_name %in% genes_with_alternate_inhibits_either$gene_name) %>%
  ggplot(aes(x = inhibits_either, y = Tm)) +
  geom_violin(aes(fill=inhibits_either),width = 0.8, alpha=0.7, draw_quantiles =  0.5) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "#d2d2d2"), limits = c("TRUE", "FALSE")) +
  geom_jitter(width = 0, height = 0, alpha = 0.5, size = 3) +
  geom_line(aes(group = gene_name), alpha = 0.5) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  labs(x = "Growth-inhibition of E. coli K12 and/or UPEC", y = "Melting temperature (C)", fill = "Inhibition of:  ",
       size=12)

Tm_plot_only_alternating_genes
ggsave( filename = "./analysis/pna_specific_predictors/Tm_plot_only_alternating_genes_paired.svg",
        plot = Tm_plot_only_alternating_genes, width = 7, height=5)

# same with double DF:
# get genes that have alternate inhibits_either values for different start positions
genes_with_alternate_inhibits_either_double <- pnas_only_double %>%
  group_by(gene_name) %>%
  summarise(
    inhibits_either = sum(inhibits_either),
    total = n()
  ) %>%
  filter(inhibits_either != total & inhibits_either != 0) %>%
  select(gene_name)

Tm_plot_only_alternating_genes_double <- pnas_only_double %>% #filter(sc_bases<5) %>%
  filter(gene_name %in% genes_with_alternate_inhibits_either_double$gene_name) %>%
  ggplot(aes(x = inhibits_either, y = Tm)) +
  geom_violin(aes(fill=inhibits_either),width = 0.8, alpha=0.7, draw_quantiles =  0.5) +
   scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "#d2d2d2"), limits = c("TRUE", "FALSE")) +
  geom_jitter(width = 0, height = 0, alpha = 0.5, size = 3) +
  geom_line(aes(group = gene_name), alpha = 0.5) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  labs(x = "Growth-inhibition of E. coli K12 and/or UPEC", y = "Melting temperature (C)", fill = "Inhibition of:  ",
       size=12)
Tm_plot_only_alternating_genes_double
ggsave( filename = "./analysis/pna_specific_predictors/only_double_plots/Tm_plot_only_alternating_genes_paired_double.svg",
        plot = Tm_plot_only_alternating_genes_double, width = 7, height=5)

# now, create a venn euler diagram to show number of PNAs inhibiting E. coli K12, the number of PNAs inhibiting UPEC,
# and the number of PNAs inhibiting both E. coli K12 and UPEC, as well as the number of PNAs inhibiting neither
# E. coli K12 nor UPEC. Use the euler function of the eulrr package for this.

# first, create a data frame with the number of PNAs inhibiting E. coli K12, the number of PNAs inhibiting UPEC,
# and the number of PNAs inhibiting both E. coli K12 and UPEC, as well as the number of PNAs inhibiting neither
# E. coli K12 nor UPEC








# write machine-learning model to predict inhibits_either from Tm, PNA_molecular_weight, total_off_targets_0mm,
# total_off_targets_1mm, total_off_targets_2mm, tir_off_targets_0mm, tir_off_targets_1mm, tir_off_targets_2mm, sc_bases,
# PNA_GC_content, start_position, purine_percentage, longest_purine_stretch, homopolymers.
# before, perform feature selection using random forest and use the top features for the model
# write model to file
library(randomForest)
library("caret")
library(pROC)
library(rpart)
library(rpart.plot)

# select features using random forest for classification
rf_model <- randomForest(factor(inhibits_either) ~ Tm + PNA_molecular_weight + total_off_targets_0mm + total_off_targets_1mm +
                           total_off_targets_2mm + tir_off_targets_0mm + tir_off_targets_1mm + tir_off_targets_2mm +
                           sc_bases + PNA_GC_content + start_position + purine_percentage + longest_purine_stretch +
                           homopolymers,
                         data = pnas_filtered, importance = TRUE, ntree = 1000)

# plot feature importance using ggplot
rf_importance <- data.frame(rf_model$importance)
rf_importance$features <- rownames(rf_importance)
rf_importance$features <- factor(rf_importance$features, levels = rf_importance$features[order(rf_importance$MeanDecreaseAccuracy, decreasing = TRUE)])
rf_importance_plot <- ggplot(rf_importance, aes(x = features, y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill="#571430") +
  theme_classic() +
  # make x axis labels 45 degrees
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 13),
    axis.title = element_text(size=15, )) +
  labs(x = "Feature", y = "Mean decrease in accuracy", title = "Feature importance using random forest")
rf_importance_plot

# plot feature importance using randomforest package
varImpPlot(rf_model, main = "Feature importance using random forest")

# now use the top 10 features for the model and generate multiple models and select the best one using cross-validation.
# write model to file

# select top 10 features
top_10_features <- c("Tm", "PNA_molecular_weight", "total_off_targets_0mm", "total_off_targets_1mm",
                     "total_off_targets_2mm", "tir_off_targets_0mm", "tir_off_targets_1mm", "tir_off_targets_2mm",
                     "sc_bases", "PNA_GC_content")

# now build multiple classification models (svg, random forest, decision tree) and select the best one using
# cross-validation. Write the best model to file.
# first, create a function that takes a model name as input and returns the model
model_function <- function(model_name) {
  if (model_name == "svm") {
    model <- train(factor(inhibits_either) ~ Tm + PNA_molecular_weight + total_off_targets_0mm + total_off_targets_1mm +
                     total_off_targets_2mm + tir_off_targets_0mm + tir_off_targets_1mm + tir_off_targets_2mm +
                     sc_bases + PNA_GC_content,
                   data = pnas_filtered, method = "svmRadial", trControl = trainControl(method = "cv", number = 10))
  } else if (model_name == "rf") {
    model <- train(factor(inhibits_either) ~ Tm + PNA_molecular_weight + total_off_targets_0mm + total_off_targets_1mm +
                     total_off_targets_2mm + tir_off_targets_0mm + tir_off_targets_1mm + tir_off_targets_2mm +
                     sc_bases + PNA_GC_content,
                   data = pnas_filtered, method = "rf", trControl = trainControl(method = "cv", number = 10))
  } else if (model_name == "dt") {
    model <- train(factor(inhibits_either) ~ Tm + PNA_molecular_weight + total_off_targets_0mm + total_off_targets_1mm +
                     total_off_targets_2mm + tir_off_targets_0mm + tir_off_targets_1mm + tir_off_targets_2mm +
                     sc_bases + PNA_GC_content,
                   data = pnas_filtered, method = "rpart", trControl = trainControl(method = "cv", number = 10))
  }
  return(model)
}

# now, create a function that takes a model name as input and returns the model's cross-validated accuracy.
# this function will be used to select the best model:
model_accuracy_function <- function(model_name) {
  if (model_name == "svm") {
    model <- train(factor(inhibits_either) ~ Tm + PNA_molecular_weight + total_off_targets_0mm + total_off_targets_1mm +
                     total_off_targets_2mm + tir_off_targets_0mm + tir_off_targets_1mm + tir_off_targets_2mm +
                     sc_bases + PNA_GC_content,
                   data = pnas_filtered, method = "svmRadial", trControl = trainControl(method = "cv", number = 10))
    accuracy <- model$results$Accuracy
  } else if (model_name == "rf") {
    model <- train(factor(inhibits_either) ~ Tm + PNA_molecular_weight + total_off_targets_0mm + total_off_targets_1mm +
                     total_off_targets_2mm + tir_off_targets_0mm + tir_off_targets_1mm + tir_off_targets_2mm +
                     sc_bases + PNA_GC_content,
                   data = pnas_filtered, method = "rf", trControl = trainControl(method = "cv", number = 10))
    accuracy <- model$results$Accuracy
  } else if (model_name == "dt") {
    model <- train(factor(inhibits_either) ~ Tm + PNA_molecular_weight + total_off_targets_0mm + total_off_targets_1mm +
                     total_off_targets_2mm + tir_off_targets_0mm + tir_off_targets_1mm + tir_off_targets_2mm +
                     sc_bases + PNA_GC_content,
                   data = pnas_filtered, method = "rpart", trControl = trainControl(method = "cv", number = 10))
    accuracy <- model$results$Accuracy
  }
  return(accuracy)
}

# now use the functions above to generate multiple models and select the best one using cross-validation
# first, generate the models
svm_model <- model_function("svm")
rf_model <- model_function("rf")
dt_model <- model_function("dt")

# now, get the cross-validated accuracy for each model
svm_accuracy <- model_accuracy_function("svm")
rf_accuracy <- model_accuracy_function("rf")
dt_accuracy <- model_accuracy_function("dt")



#plot the decision tree model using ggplot
library(rpart.plot)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)

#plot the model
p <- rpart.plot(dt_model$finalModel)
p


