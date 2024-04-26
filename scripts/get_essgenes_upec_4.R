# import gff file e_coli_k12.gff3 and extract the genes that are in our essential genes list.
# This script is used to create the file e_coli_k12_essgenes.gff3

# load libraries
library(data.table)
library(dplyr)
library(readr)
library(Biobase)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(ggbeeswarm)
library(viridis)
library(ggrepel)

# read in the gff file and skip hash lines
gff <- read.delim("./data/reference_sequences/e_coli_K12.gff3", skip = 3, header = FALSE)
gff_upec <- read.delim("./data/reference_sequences/ecoli536_sRNAs_modified.gff3", skip = 3, header = FALSE)

# read in the essential genes table
essgenes <- read.delim("./data/pnas_pna_specific_features.tsv", header = TRUE)
egenes_list <- unique(essgenes$locus_tag)

# extract locus tags from 9th column of gff
gff$locus_tag <- gsub(".*locus_tag=([^;]+).*", "\\1", gff$V9)

# get gene names from 9th column of gff
gff$gene <- gsub(".*gene=([^;]+).*", "\\1", gff$V9)

# combine gene and locus tag columns in one string
gff$gene_locus <- paste( gff$locus_tag,gff$gene, sep = "_")
gff_upec$locus_tag <- gsub(".*locus_tag=([^;]+).*", "\\1", gff_upec$V9)

# extract the essential genes from the gff file
essgenes_gff <- gff[gff$locus_tag %in% egenes_list,] %>% filter(V3 == "gene") %>% select(V1,V2, gene_locus, V4, V5, V6, V7, V8, V9)
genes_gff_upec <- gff_upec %>% filter(V3 == "gene") %>% select(V1,V2, locus_tag, V4, V5, V6, V7, V8, V9)
genes_gff_upec <- gff_upec %>% filter(V3 == "gene") %>% select(V1,V2, locus_tag, V4, V5, V6, V7, V8, V9)

# write the essential genes gff file. do not include colnames or rownames and use tab as separator. use write_de
write_delim(essgenes_gff, "./data/reference_sequences/e_coli_k12_essgenes.gff3", delim = "\t", col_names = FALSE)
write_delim(genes_gff_upec, "./data/reference_sequences/upecgenes.gff3", delim = "\t", col_names = FALSE)

# read in the blastn output file
blastn <- read_delim("./data/reference_sequences/blast_egenes.txt", delim = "\t", col_names = FALSE)
proteinortho <- read_delim("./data/reference_sequences/upec_egenes.proteinortho.tsv", delim = "\t")

mapping_proteinortho <- proteinortho %>% select(4,5) %>% mutate(K12_locus_tag = gsub("^([^_]+)_.*", "\\1", egenes_k12.fasta )) %>%
  mutate(K12_genename = gsub("^[^_]+_([^:]+):.*", "\\1", egenes_k12.fasta)) %>%
  mutate(upec_locus_tag = gsub("^([^:]+):.*", "\\1",  genes_upec.fasta )) %>% select(-egenes_k12.fasta, -genes_upec.fasta)


mappings_egenes <- blastn %>% select(1,2) %>% mutate(K12_locus_tag = gsub("^([^_]+)_.*", "\\1", X1)) %>%
  mutate(K12_genename = gsub("^[^_]+_([^:]+):.*", "\\1", X1)) %>%
  mutate(upec_locus_tag = gsub("^([^:]+):.*", "\\1", X2)) %>% select(-X1, -X2)


# get gff of essential genes from upec
upec_essgenes <- genes_gff_upec %>% filter(locus_tag %in% mappings_egenes$upec_locus_tag) %>% select(V9) %>%
    mutate(locus_tag = gsub(".*locus_tag=(.+).*", "\\1", V9)) %>% select(locus_tag)


upec_essgenes <- mapping_proteinortho %>% filter(upec_locus_tag %in% mappings_egenes$upec_locus_tag) %>% select(upec_locus_tag) %>%
    mutate(locus_tag = gsub(".*locus_tag=(.+).*", "\\1", upec_locus_tag)) %>% select(locus_tag) %>%
  unlist %>% as.character

write_csv(as.data.frame(upec_essgenes), "./data/all_genes_all_pnas/upec_essgenes.csv", col_names = FALSE)

# get ess. genes that are not in upec but in ecoli k12
genes_not_in_upec <- egenes_list[! egenes_list %in%  mapping_proteinortho$K12_locus_tag]
table_unknown_genes <-  mappings_egenes[mappings_egenes$K12_locus_tag %in% genes_not_in_upec,]

# note function of these genes:
# lolE: transmembrane protein; cohE: prophage repressor protein;  racR: Rac prophage predicted DNA-binding
# transcriptional regulator; dicA: Qin prophage predicted regulator; tyrS: tyrosyl-tRNA synthetase;
# plsB: acyltransferase:
table_unknown_genes$product <- c("transmembrane protein", "e14 prophage repressor protein",
                                 "Rac prophage predicted regulator", "Qin prophage predicted regulator",
                                 "tyrosyl-tRNA synthetase", "acyltransferase")

# add note to the table_unknown_genes file. lolE: nonfunction due to frameshift; cohE: Not core essential gene in
# Enterobacteriaceae; racR: Not core essential gene in Enterobacteriaceae; dicA: Not core essential gene in
# Enterobacteriaceae; tyrS: nonfunction due to frameshift; plsB: nonfunction due to frameshift
table_unknown_genes$note <- c("nonfunction due to frameshift in UPEC", "no core essential gene in Enterobacteriaceae",
                              "no core essential gene in Enterobacteriaceae", "no core essential gene in Enterobacteriaceae",
                              "nonfunction due to frameshift in UPEC", "nonfunction due to frameshift in UPEC")

table_unknown_genes <- table_unknown_genes[-c(1,3)]
colnames(table_unknown_genes) <- c("gene name (K12)", "product", "note")

# write the table_unknown_genes to pdf
pdf("./analysis/unconserved_essential_genes_table.pdf", width = 10, height = 10)
ggtexttable(table_unknown_genes, rows = NULL, theme = ttheme("light")) %>%
  table_cell_font(row = 1:nrow(table_unknown_genes)+1, column = 1, face = "italic")
dev.off()


# get on-targets in upec and check whether they are the same for e. coli:
# read in the on-target list:
ot_list_upec_raw <- read_delim("./data/offtargets_fulltranscripts_sorted_upec.csv")

# get the on-targets that are also in the essential genes list
ot_list_upec <- ot_list_upec_raw %>% filter(TIR=="TIR")%>% filter(num_mismatch == 0) %>% filter(locus_tag %in% upec_essgenes) %>%
  filter(!position_from_CDS_start>0)

# modify the mappings:
mappings_egenes <- mappings_egenes[mappings_egenes$upec_locus_tag %in% upec_essgenes,]

# import k12 table:
ot_list_k12 <- read_delim("./data/all_predictors_with_MICs.csv")

# add mapped k12 locus tags to ot_list_upec. use the K12_locus_tag entry from the mappings_egenes table and match with
# the locus_tag entry from the ot_list_upec table, which is called upec_locus_tag in the mappings_egenes table
ot_list_upec <- ot_list_upec %>% left_join(mappings_egenes %>% select(K12_locus_tag, upec_locus_tag),
                                           by = c("locus_tag" = "upec_locus_tag"))

# check which of the upec on-targets are mapping the same position as in k12

# create a vector including on-target locus tags, position, and probe_id from ot_list_upec
unique_pna_vector_upec <- paste0(ot_list_upec$K12_locus_tag, ot_list_upec$position_from_CDS_start, ot_list_upec$ASO_id)

# do same for k12
unique_pna_vector_k12 <- paste0(ot_list_k12$locus_tag, ot_list_k12$start_position-1, ot_list_k12$pna_name)

# get all pnas in

# get number of on-targets that are mapping the same position in k12 and upec
nr_overlapping_pnas_k12_upec <- sum(unique_pna_vector_upec %in% unique_pna_vector_k12)

# get number of pnas targeting no ess. gene in upec
pnas_targeting_no_essgene <- length(unique_pna_vector_k12)-nr_overlapping_pnas_k12_upec

# create a ggplot showing relative number of on-targets mapping the same position in k12 and upec
# create a data frame with the data
data_plot <- data.frame(group = c("Same PNA target gene in UPEC & K12", "PNA does not target the same gene in UPEC"),
                        y = c(nr_overlapping_pnas_k12_upec, pnas_targeting_no_essgene)) %>%
  mutate(prop = y / sum(y) *100) %>%
    mutate(ypos = cumsum(prop) - 0.4 * prop) %>% mutate(group= factor(group, levels = c("PNA does not target the same gene in UPEC",
                                                                                        "Same PNA target gene in UPEC & K12")))

# create the plot
piechart <- ggplot(data_plot, aes(x = "", y = prop, fill=group)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "number of on-targets", x = "") +
  coord_polar("y", start=0) + theme_void() + #theme(legend.position="none") +
  geom_text(aes(y = ypos, label = y), color = c("white"), size=6, vjust=c(7,-6.5), hjust=c(-2.5,1.5)) +
  scale_fill_manual(values = c("#66cdaa" ,"#4682b4"))
piechart

# save piechart to pdf
pdf("./analysis/piechart_on-targets.pdf", width = 10, height = 10)
piechart
dev.off()



# check, how many of the PNAs that target no essential gene in UPEC are effective in K12 or upec

# get vector of pnas with both
targets_both_upec_and_k12 <- gsub(".*(JVpna.*)","\\1",
                                  unique_pna_vector_upec[unique_pna_vector_upec %in% unique_pna_vector_k12])
targets_diff_ess_gene <- unique(ot_list_upec$probe_id) [! unique(ot_list_upec$probe_id) %in% targets_both_upec_and_k12]
targets_no_ess_gene <- ot_list_k12$pna_name[! ot_list_k12$pna_name %in% targets_both_upec_and_k12]

# get the MICs for the pnas targeting no essential gene in upec
data_targeting_no_ess_gene <- ot_list_k12[ot_list_k12$pna_name %in% targets_no_ess_gene,] %>%
  select(pna_name, MIC_K12, MIC_UPEC)

# get table with pnas targeting no essential gene in upec and their MICs in upec and k12
ot_list_no_ess_gene <- ot_list_k12[ot_list_k12$pna_name %in% c(targets_no_ess_gene, targets_diff_ess_gene),]
# add upec locus tag:
ot_list_no_ess_gene <- ot_list_no_ess_gene %>%
  mutate(upec_locus_tag = mappings_egenes$upec_locus_tag[match(locus_tag, mappings_egenes$K12_locus_tag)])

#save it as excel file
write.csv(ot_list_no_ess_gene, "./analysis/pnas_targeting_no_ess_gene.csv")

# make a plot showing the MICs of the pnas targeting no essential gene in upec. make a violin plot as the other one

# adjust df (merge MIC columns and add species column)
data_targeting_no_ess_gene <- data_targeting_no_ess_gene %>%
  pivot_longer(cols = c(MIC_K12, MIC_UPEC), names_to = "species", values_to = "MIC") %>%
  mutate(species = gsub("MIC_", "", species))%>%
  # change MICs at 500 to 20
    mutate(MIC = ifelse(MIC == 500, 20, MIC))

#create violin plot with MIC as y and species as x. add jitter points
plot_mic_untargeted <- ggplot(data_targeting_no_ess_gene, aes(y = log2(MIC), x = species)) +
  # add violin plot
    geom_violin(aes(fill=species), trim = FALSE, alpha = 1,  width=1.1) +
    # add scatterplot
    geom_jitter( size = 3,
                  alpha = 0.2,  width = 0.3, height = 0.1) +
    labs(y = "MIC", x = "") +
    theme_classic() +
    # add title
    ggtitle("MIC UPEC vs K12") +
  # add y axis ticks
    scale_y_continuous(breaks = log2(c(1.25, 2.5, 5,10,20)), labels = c("1.25", "2.5", "5","10",">10"),
                       limits = c(0,5.5)) +
      # add p-values using wilcox.test and stat_compare_means and put them on top of the plot, write as "p-value = 0.0001"
    stat_compare_means(comparisons = list(c("UPEC","K12")), method = "wilcox.test", label = "p.signif",
                       label.y = log2(20+18)) +

    theme(axis.text = element_text( size = 14),
            axis.title = element_text( size = 16),
            plot.title = element_text( size = 18, face = "bold", hjust = 0.5),
            legend.position = "none")+
  scale_fill_manual(values = c("#66cdaa", "#b5cde1"))
  
plot_mic_untargeted


#save plot to pdf
ggsave("./analysis/plot_mic_untargeted.pdf", plot_mic_untargeted, width = 6, height = 5.1)










# filter out pnas targeting no essential gene in upec
data_filterd_UPEC <- ot_list_k12 %>% mutate(ess_target_upec = pna_name %in% targets_both_upec_and_k12) %>%
  # change MICs at 500 to 15:
    mutate(MIC_UPEC = ifelse(MIC_UPEC == 500, 20, MIC_UPEC)) %>%
  # same in K12:
    mutate(MIC_K12 = ifelse(MIC_K12 == 500, 20, MIC_K12)) %>%
  # filter out pnas targeting no essential gene in upec:
    filter(ess_target_upec) %>%
  #add upec locus tag:
  mutate(upec_locus_tag = mappings_egenes$upec_locus_tag[match(locus_tag, mappings_egenes$K12_locus_tag)])

# find out what amount of PNAs have lower MIC in K12 than in UPEC, have lower MIC in UPEC than in K12, and have same
# MIC in UPEC and K12. create a piechart showing the relative amount of pnas in each category

# get the number of pnas with lower MIC in K12 than in UPEC
nr_pnas_lower_mic_k12 <- sum(data_filterd_UPEC$MIC_K12 < data_filterd_UPEC$MIC_UPEC)

# get the number of pnas with lower MIC in UPEC than in K12
nr_pnas_lower_mic_upec <- sum(data_filterd_UPEC$MIC_K12 > data_filterd_UPEC$MIC_UPEC)

# get the number of pnas with same MIC in UPEC and K12
nr_pnas_same_mic <- sum(data_filterd_UPEC$MIC_K12 == data_filterd_UPEC$MIC_UPEC)

# create a data frame with the data
data_plot <- data.frame(group = c("lower MIC in K12 than in UPEC", "lower MIC in UPEC than in K12", "same MIC in UPEC and K12"),
                        y = c(nr_pnas_lower_mic_k12, nr_pnas_lower_mic_upec, nr_pnas_same_mic)) %>%
  mutate(prop = y / sum(y) *100) %>%
    mutate(ypos = cumsum(prop) - 0.5 * prop) %>% mutate(group= factor(group, levels = c("same MIC in UPEC and K12",
                                                                                           "lower MIC in UPEC than in K12",
                                                                                        "lower MIC in K12 than in UPEC")))
# create the plot
piechart <- ggplot(data_plot, aes(x = "", y = prop, fill=group)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y = "number of on-targets", x = "") +
  coord_polar("y", start=0) + theme_void() + #theme(legend.position="none") +
  geom_text(aes(y = ypos, label = y), color = "white", size=8, vjust=c(-5,1,-3), hjust=c(0,-1.5,2)) +
  scale_fill_manual(values = c(  "#4682b4" ,"#b5cde1","#66cdaa"))
piechart

# save piechart to pdf
pdf("./analysis/piechart_mic_upec_k12_comparison.pdf", width = 10, height = 10)
piechart
dev.off()


# save as excel file
library(xlsx)
write.xlsx(data_filterd_UPEC[data_filterd_UPEC$MIC_UPEC>10 & data_filterd_UPEC$MIC_K12<=10,], "./analysis/pnas_only_inhibiting_k12.xlsx")
write.xlsx(data_filterd_UPEC[data_filterd_UPEC$MIC_K12 < data_filterd_UPEC$MIC_UPEC,], "./analysis/pnas_better_inhibiting_k12.xlsx")


# add expression data of upec to the data_filterd_UPEC table:
# read in the expression data
tpm_upec <- read.table("./data/transcriptomic_expression_K12_upec/upec_ctrl_log_tpm_plus1.csv", header = T,
                       sep = ",", row.names = 1)

# add average expression of the genes targeted by the pnas to the data_filterd_UPEC table
data_filterd_UPEC <- data_filterd_UPEC %>%
  # add average expression of the genes targeted by the pnas to the data_filterd_UPEC table
  mutate(expression_upec = rowMeans(tpm_upec[match(upec_locus_tag, rownames(tpm_upec)),])) %>%
  # remove expression and sec_structure columns
    select(-expression, -sec_structure)



# create gff file for only essential genes in UPEC from gff_upec:
gff_upec_ess <- gff_upec[gff_upec$locus_tag %in% unique(data_filterd_UPEC$upec_locus_tag),] %>%
  # get only cds
    filter(V3 == "CDS") %>%
  select(V1, V2, locus_tag, V4, V5, V6, V7, V8, V9)

gff_upec_ess_startregs <- gff_upec_ess %>%
  # change start and end to start-30 and start+30 respectively for V7=+ and end-30 and end+30 for V7=-
    mutate(V4 = ifelse(V7 == "+", V4-30, V5-30)) %>%
    mutate(V5 = ifelse(V7 == "+", V4+60, V5+30))

# save startregs to file for use in bedtools
write.table(gff_upec_ess_startregs, "./data/sec_structure_rnafold/upec_startsites.gff", sep = "\t",
            quote = F, row.names = F, col.names = F)

# pick up the MFE from ./data/sec_structure_rnafold/upec_sec_structure.tsv
sec_structure_upec <- read.table("./data/sec_structure_rnafold/upec_sec_structure.tsv", header = F, sep = "\t")

# add the MFE to the data_filterd_UPEC table
data_filterd_UPEC$MFE_UPEC <- sec_structure_upec$V4[match(data_filterd_UPEC$upec_locus_tag, sec_structure_upec$V1)]

#View(data_filterd_UPEC)

library(eulerr)
# make eler of nr of genes that have MIC>=10 in upec, k12, both, and none
nr_pnas_only_k12 <- sum(ot_list_k12$MIC_UPEC>10 & ot_list_k12$MIC_K12<=10)
nr_pnas_only_upec <- sum(ot_list_k12$MIC_UPEC<=10 & ot_list_k12$MIC_K12>10)
nr_pnas_none <- sum(ot_list_k12$MIC_UPEC>10 & ot_list_k12$MIC_K12>10)
nr_pnas_both <- sum(ot_list_k12$MIC_UPEC<=10 & ot_list_k12$MIC_K12<=10)



# create a euler plot
fit <- euler(c("UPEC" = nr_pnas_only_upec, "K12" = nr_pnas_only_k12, "UPEC&K12" = nr_pnas_both, "None" = nr_pnas_none))
svg(filename = "./analysis/pna_specific_predictors/euler_plot_normal_updated.svg", width = 8, height = 6)
plot(fit, main = "PNA-mediated growth inhibition at 10uM", quantities = TRUE,
     fills = list(fill = c("#b5cde1", "#66cdaa", "#d2d2d2", "#4682b4"),
                  alpha=1))
dev.off()

# create plot with startpos_plots
startpos_plot_function <- function(x_axis) {
    startpos_plot <- data_filterd_UPEC %>%
      # sort inhibition_of factor levels
        mutate(inhibition_of = factor(inhibition_of, levels = c("neither", "K12 only","UPEC only","both"    ))) %>%
        ggplot(aes(x = !!sym(x_axis), fill = inhibition_of)) +
        geom_bar() + ggtitle("PNA-mediated growth inhibition at 10uM") +
      # order the legend by the order of the levels in the factor
        scale_fill_manual(values = c("K12 only" = "#66cdaa", "UPEC only" = "#b5cde1", "both" = "#4682b4",
                                     "neither" = "#d2d2d2")) +
        theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 13),
    # make title bold and bigger and center it
            plot.title = element_text(size = 20, hjust = 0.5))+
        theme(axis.text.x = element_text( hjust = 0.7, vjust = 1)) +
        labs(x = x_axis, y = "Number of PNAs", fill = "Inhibition of")
    return(startpos_plot)
}


startpos_plot <- startpos_plot_function("start_position")
startpos_plot
# save plot as svg file
ggsave( filename = "./analysis/pna_specific_predictors/startpos_plot.svg", plot = startpos_plot, width = 10, height = 5)


# make plot with all MICs and show how many PNAs have that mic in upec.

# create a data frame with the data
data_plot_upec <- data_filterd_UPEC %>%
  group_by(MIC_UPEC) %>%
  summarise(nr_pnas = n())



# do same plot for K12
# create a data frame with the data
data_plot_k12 <- data_filterd_UPEC %>%
  group_by(MIC_K12) %>%
  summarise(nr_pnas = n())



# create df with both upec and k12 (bind)
data_plot_both <- data_plot_upec %>%
  mutate(species = "UPEC", MIC = MIC_UPEC) %>%
  select(species, MIC, nr_pnas) %>%
  bind_rows(data_plot_k12 %>%
              mutate(species = "K12", MIC = MIC_K12) %>%
                select(species, MIC, nr_pnas)) %>%
  # add an entry for mic 1.25 in k12 with 0 pnas
    add_row(species = "K12", MIC = 1.25, nr_pnas = 0)

# create the plot. make x axis log2 scale abd add fill for species (species next to each other, not stacked!!)
plot_mic_both <- ggplot(data_plot_both, aes(x = log2(MIC), y = nr_pnas, fill = species)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  scale_x_continuous(breaks = log2(c(1.25,2.5,5,10,20)), labels = c("1.25","2.5","5","10",">10"))+
  labs(y = "number of PNAs", x = "MIC") +
  theme_classic() +
  theme(axis.text = element_text( size = 12),
        axis.title = element_text( size = 15)) +
  #add numbers to the bars so that they are next to each other
    geom_text(aes(label = nr_pnas), position = position_dodge(width = 0.9), vjust = -0.3, size = 5) +
  # add coloring as in the other plots
    scale_fill_manual(values = c(  "#66cdaa", "#b5cde1"))
plot_mic_both

# save plot to svg
svg("./analysis/plot_mic_both_nr_pnas.svg", width = 6, height = 5.1)
plot_mic_both
dev.off()


# import data on gene essentiality for e. coli
ess_genes_ecoli <- read.table("./data/gene_essentiality/mbio_paper_2021.csv", header = T, sep = "\t")

# import other data on crispri fitness
crispr_fitness <- read.table("./data/gene_essentiality/rel_fitness_ecoli_cellsystems.csv", header = T, sep = "\t")
crispr_fitness <- crispr_fitness %>%
  filter(!is.na(relative.fitness..mean.)) %>%
  select(locus_tag,gene, relative.fitness..mean.) %>%
  # summarize locus tag and gene name with relative.fitness..mean.
    group_by(locus_tag, gene) %>%
    summarise(RF_mean_cellsytems = mean(relative.fitness..mean.))

# import mrna decay data
mrna_decay <- read.table("./data/mRNA_decay/esquerre_2015_nar.csv", header = T, sep = "\t")
mrna_decay <- mrna_decay %>%
  select("GeneId", GeneSym, mRNA.half.lifea.in.min.at.µ.0.40.h.1) %>%
  # rename column names
    rename(locus_tag = GeneId, gene = GeneSym, mRNA_half_life = mRNA.half.lifea.in.min.at.µ.0.40.h.1) %>%
  as_tibble()

# import STRING data
string_data <- read_delim("./data/protein_protein_interactions/stringdb_k12.txt", col_names = T, delim = " ")
string_high_confidence_lt <- string_data[string_data$combined_score>=700,] %>% select(protein1, protein2) %>%
  pivot_longer(cols = c(protein1, protein2), names_to = "protein", values_to = "locus_tag") %>%
  mutate(locus_tag = gsub(".*\\.", "", locus_tag)) %>% select(-protein) %>%
  unlist() %>% as.character
# for all lts, get the count of occurences in the string_high_confidence_lt vector
string_high_confidence_lt_count <- table(string_high_confidence_lt)


# import crispri data from rousset 2018
crispr_fitness_rousset <- read_delim("./data/gene_essentiality/plos_gen_rousset_2018.csv", col_names = T, delim = ",") %>%
  select(gene, median_coding)

#import crispri_a_baum (excel)
library(readxl)
crispr_fitness_a_baum <- read_excel("./data/gene_essentiality/a_baum_crispri_2024_mbio.xlsx")%>%
  filter(grepl("None_0_T2 - None_0_T0", condition)) %>%
    select("unique name", "median LFC")
colnames(crispr_fitness_a_baum) <- c("gene", "median_log2FC_a_baum")

# add protein half life excel
protein_half_life <- read_excel("./data/protein_half_life/msystms_nagar_etal_2021.xlsx") %>%
  select("Gene names", half_life) %>%
  mutate(half_life = as.numeric(half_life)) %>%
  mutate(gene_name = `Gene names`) %>%
  select(gene_name, half_life)%>%
  # filter inf half life
  filter(half_life != Inf)


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

# add RF.mean to the data_filterd_UPEC table by gene name
data_filterd_UPEC <- data_filterd_UPEC %>%
  mutate(RF_mean_u = ess_genes_ecoli$RF.mean_u[match(gene_name, ess_genes_ecoli$Gene)]) %>%
  mutate(RF_mean_i = ess_genes_ecoli$RF.mean_i[match(gene_name, ess_genes_ecoli$Gene)]) %>%
  # make mean from NA or NaN
    mutate(RF_mean_u = ifelse(is.na(RF_mean_u), mean(na.omit(ess_genes_ecoli$RF.mean_u)), RF_mean_u)) %>%
    mutate(RF_mean_i = ifelse(is.na(RF_mean_i), mean(na.omit(ess_genes_ecoli$RF.mean_i)), RF_mean_i)) %>%
      # add RF_mean_cellsytems to the data_filterd_UPEC table by gene name
      mutate(RF_mean_cellsytems = crispr_fitness$RF_mean_cellsytems[match(gene_name, crispr_fitness$gene)]) %>%
      # make mean from NA or NaN
      mutate(RF_mean_cellsytems = ifelse(is.na(RF_mean_cellsytems), mean(na.omit(crispr_fitness$RF_mean_cellsytems)),
                                         RF_mean_cellsytems)) %>%
  # # rename RF_mean_cellsystems to gene_vulnerability_crispri and multiply it by -1
       rename(gene_vulnerability_crispri = RF_mean_cellsytems) %>%
          mutate(gene_vulnerability_crispri = gene_vulnerability_crispri * (-1) +1) %>%
  # round vulnerability to 2 digits after comma
          mutate(gene_vulnerability_crispri = round(gene_vulnerability_crispri, 2)) %>%
  # round expression to 2 digits after comma
          mutate(expression_upec = round(expression_upec, 2)) %>%
        # add mRNA_half_life to the data_filterd_UPEC table by gene name
          mutate(mRNA_half_life = mrna_decay$mRNA_half_life[match(gene_name, mrna_decay$gene)]) %>%
            # make mean from NA or NaN
  mutate(mRNA_half_life = ifelse(is.na(mRNA_half_life), median(na.omit(mrna_decay$mRNA_half_life[match(gene_name, mrna_decay$gene)])), mRNA_half_life)) %>%
  # add string high confidence interactions to the data_filterd_UPEC table by gene name
    mutate(string_interactions = string_high_confidence_lt_count[locus_tag]) %>%
    # make mean from NA or NaN
    mutate(string_interactions = ifelse(is.na(string_interactions), median(na.omit(as.numeric(string_high_confidence_lt_count[locus_tag]))), string_interactions)) %>%
  # add rousset 2018 data to the data_filterd_UPEC table by gene name
    mutate(crispri_log2FC_rousset = crispr_fitness_rousset$median_coding[match(gene_name, crispr_fitness_rousset$gene)]) %>%
    # make mean from NA or NaN
    mutate(crispri_log2FC_rousset = ifelse(is.na(crispri_log2FC_rousset), -median(na.omit(crispr_fitness_rousset$median_coding[match(gene_name, crispr_fitness_rousset$gene)])), -crispri_log2FC_rousset)) %>%
    # add a_baum 2024 data to the data_filterd_UPEC table by gene name
    mutate(crispri_log2FC_a_baum = crispr_fitness_a_baum$median_log2FC_a_baum[match(gene_name, crispr_fitness_a_baum$gene)]) %>%
    # make median from NA or NaN
    mutate(crispri_log2FC_a_baum = ifelse(is.na(crispri_log2FC_a_baum), -median(na.omit(crispr_fitness_a_baum$median_log2FC_a_baum[match(gene_name, crispr_fitness_a_baum$gene)])), -crispri_log2FC_a_baum)) %>%
    # add protein half life to the data_filterd_UPEC table by gene name
    mutate(protein_half_life_elife = protein_half_life$half_life[match(gene_name, protein_half_life$gene_name)]) %>%
    # make mean from NA or NaN
    mutate(protein_half_life_elife = ifelse(is.na(protein_half_life_elife), log2(median(na.omit(protein_half_life$half_life[match(gene_name, protein_half_life$gene_name)]))), log2(protein_half_life_elife))) %>%
    # add cas13 data to the data_filterd_UPEC table by gene name
    mutate(cas13_log2FC = -cas13_data$mean_log2FC[match(gene_name, cas13_data$gene)])









# Now I want to create some plots which show  as a violin plot the distribution of different factors as y and the
# MIC as x. I first have to make MIC a factor and call 20 >10.
data_filterd_UPEC_plot <- data_filterd_UPEC %>%
  # change 20 to >10
    mutate(MIC_UPEC = factor(ifelse(MIC_UPEC == 20, ">10", MIC_UPEC), levels = c("1.25","2.5","5","10",">10")))

# create violin plots for the different factors. make a function for this
violin_plot <- function(factor, ylab, title){
  #create comparisons for stat_compare_means
    comparisons <- list(c("1.25","2.5"), c("1.25","5"), c("1.25","10"), c("1.25",">10"),
                        c("2.5","5"), c("2.5","10"), c("2.5",">10"),
                        c("5","10"), c("5",">10"),
                        c("10",">10"))
  # get the p-values of the comparisons. use factor
    p_values <- list()
    for (i in 1:length(comparisons)){
      p_values[[i]] <- wilcox.test(factor[data_filterd_UPEC_plot$MIC_UPEC == comparisons[[i]][1]],
                                   factor[data_filterd_UPEC_plot$MIC_UPEC == comparisons[[i]][2]])$p.value
    }
  pvals_logic <- sapply(p_values,function(x) x <0.05)


  plot <- ggplot(data_filterd_UPEC_plot, aes(x = MIC_UPEC, y = factor)) +
    # add violin plot
    geom_violin(color = "black", alpha=0.5, width = 1, fill="steelblue") +
    # add scatterplot
      geom_jitter( size = 2, alpha = 0.5, color = "black", width = 0.2, height = 0) +
    # add boxplot
      geom_boxplot(width=0.2, fill = "white", color = "black", alpha=0.7, outlier.shape = NA) +
    # add p-values using wilcox.test and stat_compare_means and put them on top of the plot, write as "p-value = 0.0001"
        stat_compare_means(comparisons = comparisons[pvals_logic], method = "wilcox.test", label = "p.signif") +
    labs(y = ylab, x = "MIC") +
    #scale_fill_manual(values= viridis(8, direction = -1, option = "inferno")[2:6]) +
    theme_classic() +
    # add title
    ggtitle(title) +
    theme(axis.text = element_text( size = 14),
          axis.title = element_text( size = 16),
            plot.title = element_text( size = 18, face = "bold", hjust = 0.5),
    legend.position = "none")
  return(plot)
}

# create the plots
plot_mRNA_half_life <- violin_plot(data_filterd_UPEC_plot$mRNA_half_life, "mRNA half life", "mRNA half life")
plot_mRNA_half_life
plot_gene_vulnerability_crispri <- violin_plot(data_filterd_UPEC_plot$gene_vulnerability_crispri, "gene vulnerability (CRISPRi)", "gene vulnerability (CRISPRi)")
plot_gene_vulnerability_crispri
plot_Tm <- violin_plot(data_filterd_UPEC_plot$Tm, "Tm", "PNA/mRNA melting temperature (°C)")
plot_Tm
plot_ot_tot_2mm <- violin_plot(data_filterd_UPEC_plot$upec_total_off_targets_2mm, "total off targets 2mm", "total off targets 2mm")
plot_ot_tot_2mm
plot_ot_tot_1mm <- violin_plot(data_filterd_UPEC_plot$upec_total_off_targets_1mm, "total off targets 1mm", "total off targets 1mm")
plot_ot_tot_1mm
plot_ot_tot_0mm <- violin_plot(data_filterd_UPEC_plot$upec_total_off_targets_0mm, "total off targets 0mm", "total off targets 0mm")
plot_ot_tot_0mm
plot_tir_off_targets_2mm <- violin_plot(data_filterd_UPEC_plot$upec_tir_off_targets_2mm, "tir off targets 2mm", "TIR off targets 2mm")
plot_tir_off_targets_2mm
plot_tir_off_targets_1mm <- violin_plot(data_filterd_UPEC_plot$upec_tir_off_targets_1mm, "tir off targets 1mm", "TIR off targets 1mm")
plot_tir_off_targets_1mm
plot_tir_off_targets_0mm <- violin_plot(data_filterd_UPEC_plot$upec_tir_off_targets_0mm, "tir off targets 0mm", "TIR off targets 0mm")
plot_tir_off_targets_0mm
plot_sc_bases <- violin_plot(data_filterd_UPEC_plot$sc_bases, "sc bases", "self complementary bases")
plot_sc_bases
plot_mfe <- violin_plot(data_filterd_UPEC_plot$MFE_UPEC, "mfe", "minimum free energy in TIR")
plot_mfe
plot_expression_upec <- violin_plot(data_filterd_UPEC_plot$expression_upec, "expression upec", "gene expression in UPEC")
plot_expression_upec
plot_gene_length <- violin_plot(data_filterd_UPEC_plot$gene_length, "gene length", "gene length")
plot_gene_length
plot_purine_percentage <- violin_plot(data_filterd_UPEC_plot$purine_percentage, "purine percentage", "purine %")
plot_purine_percentage
plot_molecular_weight <- violin_plot(data_filterd_UPEC_plot$PNA_molecular_weight, "PNA molecular weight", "PNA molecular weight")
plot_molecular_weight
plot_longest_purine_stretch <- violin_plot(data_filterd_UPEC_plot$longest_purine_stretch, "longest purine stretch", "longest purine stretch")
plot_longest_purine_stretch

plot_string <- violin_plot(data_filterd_UPEC_plot$string_interactions, "# STRING interactions", "STRING interactions (high confidence)")
plot_string
plot_rousset_rf <- violin_plot(data_filterd_UPEC_plot$crispri_log2FC_rousset, "CRISPRi log2FC", "CRISPRi depletion Rousset 2018")
plot_rousset_rf
plot_a_baum_rf <- violin_plot(data_filterd_UPEC_plot$crispri_log2FC_a_baum, "CRISPRi log2FC", "CRISPRi depletion A. Baum 2024")
plot_a_baum_rf
plot_protein_half_life <- violin_plot(data_filterd_UPEC_plot$protein_half_life_elife, "protein half life", "protein half life")
plot_protein_half_life
plot_cas13_log2FC <- violin_plot(data_filterd_UPEC_plot$cas13_log2FC, "Cas13d -log2FC", "Cas13d -log2FC")
plot_cas13_log2FC
ggsave( "./analysis/cas13_log2FC.svg",plot_cas13_log2FC, width = 6, height = 6)


# save plots in one big pdf using cowplot
library(cowplot)
# save plots in one big pdf using cowplot
pdf("./analysis/pnas_predictors_mic_upec.pdf", width = 18, height = 28)
plot_grid(plot_Tm, plot_sc_bases,   plot_expression_upec,
          plot_ot_tot_2mm, plot_ot_tot_1mm, plot_ot_tot_0mm,
          plot_tir_off_targets_2mm, plot_tir_off_targets_1mm, plot_tir_off_targets_0mm,
          plot_gene_vulnerability_crispri, plot_mRNA_half_life, plot_mfe,
          plot_purine_percentage, plot_molecular_weight, plot_longest_purine_stretch,
            plot_string, plot_rousset_rf, plot_a_baum_rf,
          ncol = 3, nrow = 6, scale = 0.92 , labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J","K","L","M",
                                                        "N", "O"),
      label_size = 20)
dev.off()

# save plots in one big svg using cowplot
svg("./analysis/pnas_predictors_mic_upec.svg", width = 18, height = 26)
plot_grid(plot_Tm, plot_sc_bases,   plot_expression_upec,
          plot_ot_tot_2mm, plot_ot_tot_1mm, plot_ot_tot_0mm,
          plot_tir_off_targets_2mm, plot_tir_off_targets_1mm, plot_tir_off_targets_0mm,
          plot_gene_vulnerability_crispri, plot_mRNA_half_life, plot_mfe,
          plot_purine_percentage, plot_molecular_weight, plot_longest_purine_stretch,
          ncol = 3, nrow = 5, scale = 0.92 , labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J","K","L","M",
                                                        "N", "O"),
      label_size = 20)
dev.off()

# create one only for Tm, sc_bases and expression_upec
svg("./analysis/pnas_predictors_mic_upec_2.svg", width = 18, height = 6)
plot_grid(plot_Tm, plot_sc_bases, plot_expression_upec,
          ncol = 3, nrow = 1, scale = 0.92 , labels = c("A", "B", "C"),
          label_size = 20)
dev.off()

# create plots for additional factors:
plot_upstream_genes_operon <- violin_plot(data_filterd_UPEC_plot$upstream_genes_operon, "upstream genes operon", "upstream genes operon")
plot_upstream_genes_operon
plot_essential_genes_upstream <- violin_plot(data_filterd_UPEC_plot$essential_genes_upstream, "essential genes upstream", "essential genes upstream")
plot_essential_genes_upstream
plot_essential_genes_downstream <- violin_plot(data_filterd_UPEC_plot$essential_genes_downstream, "essential genes downstream", "essential genes downstream")
plot_essential_genes_downstream
plot_downstream_genes_operon <- violin_plot(data_filterd_UPEC_plot$downstream_genes_operon, "downstream genes operon", "downstream genes operon")
plot_downstream_genes_operon
plotnr_operon <- violin_plot(data_filterd_UPEC_plot$nr_genes_operon, "nr genes operon", "nr genes operon")
plotnr_operon

# get overlap plot
plot_overlap <- data_filterd_UPEC_plot %>% ggplot(aes(x = overlapping_genes, y = MIC_UPEC)) +
  # point plot with jittering
    geom_point(aes(color = overlapping_genes), position = position_jitter(width = 0.4, height = 0.2))


plot_overlap


# save the opron plots in one pdf, 2 rows and 2 columns
pdf("./analysis/pnas_predictors_mic_upec_operon.pdf", width = 12, height = 15)
plot_grid(plot_upstream_genes_operon, plot_essential_genes_upstream,
          plot_downstream_genes_operon, plot_essential_genes_downstream,
          plotnr_operon, plot_overlap,
          ncol = 2, nrow = 3, scale = 0.92 , labels = c("A", "B", "C", "D", "E", "F"),
          label_size = 20)
dev.off()




# write marked_violin as a function where the user can specify the factor to plot, the title, the top or bottom genes,
# the number of top and bottom genes and the color of the dots
marked_violin_function <- function(fact, title, top_bottom, n, color, laby){
  if (top_bottom == "top"){
    top_vuln_genes <- data_filterd_UPEC_plot %>%
      arrange(desc(.data[[fact]])) %>%
      head(n)
    wo_top_vuln_genes <- data_filterd_UPEC_plot %>%
      arrange(desc(.data[[fact]])) %>%
      # remove the top 5 genes
      slice(-1:-n)
  } else if (top_bottom == "bottom"){
    top_vuln_genes <- data_filterd_UPEC_plot %>%
      arrange(.data[[fact]]) %>%
      head(n)
    wo_top_vuln_genes <- data_filterd_UPEC_plot %>%
      arrange(.data[[fact]]) %>%
      # remove the top 5 genes
      slice(-1:-n)
  }
  #create comparisons for stat_compare_means
  comparisons <- list(c("1.25","2.5"), c("1.25","5"), c("1.25","10"), c("1.25",">10"),
                      c("2.5","5"), c("2.5","10"), c("2.5",">10"),
                      c("5","10"), c("5",">10"),
                      c("10",">10"))
  # get the p-values of the comparisons. use factor
  p_values <- list()
  for (i in 1:length(comparisons)){
    p_values[[i]] <- wilcox.test(data_filterd_UPEC_plot[[fact]][data_filterd_UPEC_plot$MIC_UPEC == comparisons[[i]][1]],
                                 data_filterd_UPEC_plot[[fact]][data_filterd_UPEC_plot$MIC_UPEC == comparisons[[i]][2]])$p.value
  }
  pvals_logic <- sapply(p_values,function(x) x <0.05)

  marked_violin <- ggplot(data_filterd_UPEC_plot, aes(x = MIC_UPEC, y = .data[[fact]])) +
    # add violin plot
    geom_violin(color = "black", alpha=0.2, width = 1, fill="steelblue") +
    # add scatterplot
    geom_jitter( data = wo_top_vuln_genes, size = 2, alpha = 0.5, color = "black", width = 0.2, height = 0) +
    # add boxplot
    geom_boxplot(width=0.2, fill = "white", color = "black", alpha=0.7, outlier.shape = NA) +
    # add p-values using wilcox.test and stat_compare_means and put them on top of the plot, write as "p-value = 0.0001"
    # make threshold of 0.05 for p-values
    stat_compare_means(comparisons = comparisons[pvals_logic], method = "wilcox.test", label = "p.signif") +
    labs(y = laby, x = "MIC") +
    theme_classic() +
    # add title
    ggtitle(title) +
    theme(axis.text = element_text( size = 14),
          axis.title = element_text( size = 16),
          plot.title = element_text( size = 18, face = "bold", hjust = 0.5),
          legend.position = "none") +
    # add the top 5 genes as red dots
    geom_point(data = top_vuln_genes, aes(x = MIC_UPEC, y = .data[[fact]]), color = color, alpha=0.7, size = 2,
               position = position_jitter(width = 0.2, height = 0, seed = 1)) +
    # add the gene names as text
    geom_label_repel(data = top_vuln_genes, aes(x = MIC_UPEC, y = .data[[fact]], label = gene_name),
                    color = color, size = 3,fontface = "bold.italic", alpha=0.7, max.overlaps = 1000,
                    position = position_jitter(width = 0.2, height = 0, seed = 1), box.padding = 0.5)

    return(marked_violin)
}

# create the plot for the top 10 genes with the highest gene vulnerability
plot_top_vuln_genes <- marked_violin_function("gene_vulnerability_crispri", "gene vuln. CRISPRi - Hawkins 2020",
                                                                           "top", 10, "darkred",
                                              "gene vulnerability (CRISPRi)")
# same for crispr cas13 log2FC
plot_top_cas13_genes <- marked_violin_function("cas13_log2FC", "Cas13d log2FC - Adler 2023",
                                              "top", 10, "darkred",
                                              "Cas13d -log2FC")
# same for plot_rousset_rf rousset et al. 2018
plot_top_rousset_genes <- marked_violin_function("crispri_log2FC_rousset", "CRISPRi log2FC - Rousset 2018",
                                              "top", 10, "darkred",
                                              "CRISPRi log2FC")
# same for a baum et al. 2024
plot_top_a_baum_genes <- marked_violin_function("crispri_log2FC_a_baum", "CRISPRi log2FC - Ward 2024",
                                              "top", 10, "darkred",
                                              "CRISPRi log2FC")

#  same for the string interactions
plot_top_string_genes <- marked_violin_function("string_interactions", "# STRING interactions",
                                              "top", 10, "darkred",
                                              "# STRING interactions")

# same for expression (popella)
plot_top_expression_genes <- marked_violin_function("expression_upec", "mRNA level - Popella 2022",
                                              "top", 10, "darkred",
                                              "gene expression (TPM)")

# same for half life
plot_top_mRNA_half_life_genes <- marked_violin_function("mRNA_half_life", "mRNA half life - Esquerré 2015",
                                              "top", 10, "darkred",
                                              "mRNA half life (min)")

# same for mfe
plot_top_mfe_genes <- marked_violin_function("MFE_UPEC", "MFE in TIR",
                                              "bottom", 10, "darkred",
                                              "MFE (delta E)")

# put all in 1 plot using cowplot and save as pdf
pdf("./analysis/pnas_predictors_mic_upec_top_genes_top10.pdf", width = 10, height = 20)
plot_grid(plot_top_vuln_genes, plot_top_cas13_genes, plot_top_rousset_genes,
          plot_top_a_baum_genes, plot_top_string_genes, plot_top_expression_genes,
          plot_top_mRNA_half_life_genes, plot_top_mfe_genes,
          ncol = 2, nrow = 4, scale = 0.92 , labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          label_size = 20)
dev.off()



# create a simple violin plot which shows distribution of MICs for upec and k12 as x axis and mic as y axis
# create a dataframe with the MICs of upec and k12

df_new_plot <- data.frame(MIC = c(data_filterd_UPEC$MIC_UPEC, data_filterd_UPEC$MIC_K12),
                          strain = factor(c(rep("UPEC", length(data_filterd_UPEC$MIC_UPEC)),
                                     rep("K12", length(data_filterd_UPEC$MIC_K12))), levels = c("UPEC", "K12")))
# create the plot
plot <- ggplot(df_new_plot, aes(y = log2(MIC), x = strain)) +
  # add violin plot
    geom_violin(aes(fill=strain), trim = FALSE, alpha = 1,  width=1.1) +
    # add scatterplot
    geom_jitter( size = 3,
                  alpha = 0.2,  width = 0.3, height = 0.1) +
    labs(y = "MIC", x = "") +
    theme_classic() +
    # add title
    ggtitle("MIC UPEC vs K12") +
  # add y axis ticks
    scale_y_continuous(breaks = log2(c(1.25, 2.5, 5,10,20)), labels = c("1.25", "2.5", "5","10",">10"),
                       limits = c(0,5.5)) +
      # add p-values using wilcox.test and stat_compare_means and put them on top of the plot, write as "p-value = 0.0001"
    stat_compare_means(comparisons = list(c("UPEC","K12")), method = "wilcox.test", label = "p.signif",
                       label.y = log2(20+18)) +

    theme(axis.text = element_text( size = 14),
            axis.title = element_text( size = 16),
            plot.title = element_text( size = 18, face = "bold", hjust = 0.5),
            legend.position = "none")+
  scale_fill_manual(values = c( "#b5cde1", "#66cdaa"))
plot

# save plot
ggsave("./analysis/pnas_mic_upec_vs_k12.pdf", plot, width = 6, height = 6)
# as svg
ggsave("./analysis/pnas_mic_upec_vs_k12.svg", plot, width = 6, height = 6)


# check for all PNA sequence bast its abundance across whole dataset and make a heatmap from it (columns are position,
# 2 rows, one of PNAs with MIC < 5 and one of PNAs with MIC > 5, and one of all PNAs). The heatmap should be colored
# by the abundance of the PNA sequence in the dataset (the more often the PNA sequence is found in the dataset, the
# darker the color should be).

# create a dataframe with the PNA sequences and their abundance per position. column 1 is position in PNA, column 2 is
# A abundance, column 3 is C abundance, column 4 is G abundance, column 5 is T abundance. Use all PNAs first
df_PNA_abundance <- data.frame(A = rep(0,9), T = rep(0,9), G = rep(0,9), C = rep(0,9))

# loop over all PNAs and count the abundance of each base per position. save as vector of length 9 for each base
for (i in 1:nrow(data_filterd_UPEC_plot)){
  for (j in 1:9){
    if (substr(data_filterd_UPEC_plot$pna_sequence[i], j, j) == "A"){
      df_PNA_abundance$A[j] <- df_PNA_abundance$A[j] + 1
    } else if (substr(data_filterd_UPEC_plot$pna_sequence[i], j, j) == "C"){
      df_PNA_abundance$C[j] <- df_PNA_abundance$C[j] + 1
    } else if (substr(data_filterd_UPEC_plot$pna_sequence[i], j, j) == "G"){
      df_PNA_abundance$G[j] <- df_PNA_abundance$G[j] + 1
    } else if (substr(data_filterd_UPEC_plot$pna_sequence[i], j, j) == "T"){
      df_PNA_abundance$T[j] <- df_PNA_abundance$T[j] + 1
    }
  }
}
# transpose df and add column names as position
df_PNA_abundance <- data.frame(t(df_PNA_abundance))
colnames(df_PNA_abundance) <- c(1:9)

# make relative abundance per position by dividing by sum of pnas (512)
df_PNA_abundance <- df_PNA_abundance/512



# create a complexheatmap object from the dataframe. use viridis inferno color scheme
library(ComplexHeatmap)
library(viridis)

# create heatmap
heatmap_PNA_abundance <- Heatmap(as.matrix(df_PNA_abundance), name = "rel. ab.", col = inferno(100),
                                 cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE,
                                 column_title = "PNA base position",  column_title_gp = gpar(fontsize = 14),
                                 column_title_side = "bottom", row_title_side = "right",
                                 column_names_rot = 0, row_names_side = "left",
                                 row_title_gp = gpar(fontsize = 14))
heatmap_PNA_abundance

# save heatmap as pdf
pdf("./analysis/pnas_abundance_heatmap.pdf", width = 4, height = 2)
draw(heatmap_PNA_abundance)
dev.off()


# crate same for PNAs with MIC <= 5
# create a dataframe with the PNA sequences and their abundance per position. column 1 is position in PNA, column 2 is
# A abundance, column 3 is C abundance, column 4 is G abundance, column 5 is T abundance. Use all PNAs first
df_PNA_abundance_5 <- data.frame(A = rep(0,9), T = rep(0,9), G = rep(0,9), C = rep(0,9))

# loop over all PNAs and count the abundance of each base per position. save as vector of length 9 for each base
for (i in 1:nrow(data_filterd_UPEC)){
  if (data_filterd_UPEC$MIC_UPEC[i] <= 5){
    for (j in 1:9){
      if (substr(data_filterd_UPEC$pna_sequence[i], j, j) == "A"){
        df_PNA_abundance_5$A[j] <- df_PNA_abundance_5$A[j] + 1
      } else if (substr(data_filterd_UPEC$pna_sequence[i], j, j) == "C"){
        df_PNA_abundance_5$C[j] <- df_PNA_abundance_5$C[j] + 1
      } else if (substr(data_filterd_UPEC$pna_sequence[i], j, j) == "G"){
        df_PNA_abundance_5$G[j] <- df_PNA_abundance_5$G[j] + 1
      } else if (substr(data_filterd_UPEC$pna_sequence[i], j, j) == "T"){
        df_PNA_abundance_5$T[j] <- df_PNA_abundance_5$T[j] + 1
      }
    }
  }
}

# transpose df and add column names as position
df_PNA_abundance_5 <- data.frame(t(df_PNA_abundance_5))
colnames(df_PNA_abundance_5) <- c(1:9)

# make relative abundance per position by dividing by sum of pnas (211)
df_PNA_abundance_5 <- df_PNA_abundance_5/211

# create heatmap
heatmap_PNA_abundance_5 <- Heatmap(as.matrix(df_PNA_abundance_5), name = "rel. ab.", col = inferno(100),
                                 cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE,
                                 column_title = "PNA base position",  column_title_gp = gpar(fontsize = 14),
                                 column_title_side = "bottom", row_title_side = "right",
                                 column_names_rot = 0, row_names_side = "left",
                                 row_title_gp = gpar(fontsize = 14))
heatmap_PNA_abundance_5

# save heatmap as pdf
pdf("./analysis/pnas_abundance_heatmap_5.pdf", width = 4, height = 2)
draw(heatmap_PNA_abundance_5)
dev.off()


# get contrast of PNAs with MIC <= 5 vs all PNAs
df_contrast_5 <- df_PNA_abundance_5 - df_PNA_abundance

# create color scheme for contrast heatmap. darkorange is positive, steelblue is negative, white is 0
library(RColorBrewer)

# create color scheme
my_palette <- colorRampPalette(c("steelblue", "white", "darkorange"))(n = 100)
# create heatmap
heatmap_contrast_5 <- Heatmap(as.matrix(df_contrast_5), name = "rel. ab.", col = my_palette,
                              cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE,
                              column_title = "PNA base position",  column_title_gp = gpar(fontsize = 14),
                              column_title_side = "bottom", row_title_side = "right",
                              column_names_rot = 0, row_names_side = "left",
                              row_title_gp = gpar(fontsize = 14))
heatmap_contrast_5

# save heatmap as pdf
pdf("./analysis/pnas_abundance_contrast_5.pdf", width = 4, height = 2)
draw(heatmap_contrast_5)
dev.off()


# create same for PNAs with MIC > 5
# create a dataframe with the PNA sequences and their abundance per position. column 1 is position in PNA, column 2 is
# A abundance, column 3 is C abundance, column 4 is G abundance, column 5 is T abundance. Use all PNAs first
df_PNA_abundance_10 <- data.frame(A = rep(0,9), T = rep(0,9), G = rep(0,9), C = rep(0,9))

# loop over all PNAs and count the abundance of each base per position. save as vector of length 9 for each base
for (i in 1:nrow(data_filterd_UPEC)){
  if (data_filterd_UPEC$MIC_UPEC[i] > 5){
    for (j in 1:9){
      if (substr(data_filterd_UPEC$pna_sequence[i], j, j) == "A"){
        df_PNA_abundance_10$A[j] <- df_PNA_abundance_10$A[j] + 1
      } else if (substr(data_filterd_UPEC$pna_sequence[i], j, j) == "C"){
        df_PNA_abundance_10$C[j] <- df_PNA_abundance_10$C[j] + 1
      } else if (substr(data_filterd_UPEC$pna_sequence[i], j, j) == "G"){
        df_PNA_abundance_10$G[j] <- df_PNA_abundance_10$G[j] + 1
      } else if (substr(data_filterd_UPEC$pna_sequence[i], j, j) == "T"){
        df_PNA_abundance_10$T[j] <- df_PNA_abundance_10$T[j] + 1
      }
    }
  }
}

# transpose df and add column names as position
df_PNA_abundance_10 <- data.frame(t(df_PNA_abundance_10))
colnames(df_PNA_abundance_10) <- c(1:9)

# make relative abundance per position by dividing by sum of pnas (301)
df_PNA_abundance_10 <- df_PNA_abundance_10/301

# create heatmap
heatmap_PNA_abundance_10 <- Heatmap(as.matrix(df_PNA_abundance_10), name = "rel. ab.", col = inferno(100),
                                 cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE,
                                 column_title = "PNA base position",  column_title_gp = gpar(fontsize = 14),
                                 column_title_side = "bottom", row_title_side = "right",
                                 column_names_rot = 0, row_names_side = "left",
                                 row_title_gp = gpar(fontsize = 14))
heatmap_PNA_abundance_10

# save heatmap as pdf
pdf("./analysis/pnas_abundance_heatmap_10.pdf", width = 4, height = 2)
draw(heatmap_PNA_abundance_10)
dev.off()

# make contrast heatmap btw PNAs with MIC <= 5 and MIC > 5
df_contrast_10 <- df_PNA_abundance_10 - df_PNA_abundance_5

# create heatmap
heatmap_contrast_10 <- Heatmap(as.matrix(df_contrast_10), name = "rel. ab.", col = my_palette,
                              cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE,
                              column_title = "PNA base position",  column_title_gp = gpar(fontsize = 14),
                              column_title_side = "bottom", row_title_side = "right",
                              column_names_rot = 0, row_names_side = "left",
                              row_title_gp = gpar(fontsize = 14))
heatmap_contrast_10

# save heatmap as pdf
pdf("./analysis/pnas_abundance_contrast_10.pdf", width = 4, height = 2)
draw(heatmap_contrast_10)
dev.off()





# check CAI for essential genes and check whether theres a correlation between cai and expression_upec
CAI_data <- read_delim("./data/CAI_calculation/UPEC_CAI.tsv", delim="\t", col_names = F)

# add CAI to data_filterd_UPEC by upec_locus_tag and X1
data_filterd_UPEC <- data_filterd_UPEC %>%
  mutate(CAI = CAI_data$X2[match(upec_locus_tag, CAI_data$X1)])

# plot correlation between CAI and expression_upec. add a correlation R value and p-value
plot_CAI <- data_filterd_UPEC %>%
  select(CAI, expression_upec, upec_locus_tag) %>%
  # keep only one entry per gene
    distinct() %>%
  ggplot(aes(x = CAI, y = expression_upec)) +
  geom_point() + xlim(c(0.63, 0.86))+
  labs(y = "expression upec", x = "Codon Adaption Index (CAI)") +
  theme_classic() +
  theme(axis.text = element_text( size = 12),
        axis.title = element_text( size = 15)) +
  # add correlation R value and p-value
  annotate("text", x = 0.64, y = 3.2,
           label = paste("R = ",
                         round(cor(data_filterd_UPEC$CAI, data_filterd_UPEC$expression_upec),2),
                         "\n                  p-value = ",
                         round(cor.test(data_filterd_UPEC$CAI, data_filterd_UPEC$expression_upec)$p.value,15)),
           size = 5)
plot_CAI

# save plot to pdf
pdf("./data/CAI_calculation/plot_CAI_expression_upec.pdf", width = 6, height = 5.1)
plot_CAI
dev.off()

x
# constrain the data frame to important features, delete columns "gene_length", "e_coli_k12_inhibition", "upec_inhibition",
# "inhibits_either", "ess_target_upec" ,"inhibition_of", "RF_mean_u", "RF_mean_i"
data_filterd_UPEC_final <- data_filterd_UPEC %>%
  select(-e_coli_k12_inhibition, -upec_inhibition, -inhibits_either, -ess_target_upec, -inhibition_of,
         -RF_mean_u, -RF_mean_i, -total_off_targets_0mm, -total_off_targets_1mm, -total_off_targets_2mm, -gene_length,
                -tir_off_targets_0mm, -tir_off_targets_1mm, -tir_off_targets_2mm)

# write data_filterd_UPEC to file
write.table(data_filterd_UPEC_final, "./data/pnas_predictors_mic_upec.tsv", sep = "\t", quote = F, row.names = F)

# make data for corr matrix
cor_data <- data_filterd_UPEC_final %>%
  # remove upec locus tag, gene name, pna sequence, target seq, all one-hot encoded bases
    select(-upec_locus_tag, -gene_name, -pna_sequence, -target_seq, -pna_name,  -A1, -A2, -A3, -A4, -A5, -A6, -A7, -A8, -A9,
             -C1, -C2, -C3, -C4, -C5, -C6, -C7, -C8, -C9, -G1, -G2, -G3, -G4, -G5, -G6, -G7, -G8, -G9, -T1, -T2, -T3,
             -T4, -T5, -T6, -T7, -T8, -T9)



# now create a correlation plot for all features. Use ComplexHeatmap for this. First create a correlation matrix
cor_matrix <- cor(cor_data[,-1])

# create a color scheme for the correlation plot. darkred is  positive (limit 1), steelblue is negative (limit -1),
# white is 0
library(circlize)
my_palette <- colorRamp2(c(-1,0,1),c("steelblue", "white", "darkred"))

# create the correlation plot
cor_plot <- Heatmap(cor_matrix, name = "correlation", col = my_palette,
                    cluster_rows = T, cluster_columns = T, show_column_names = TRUE,
                    column_title_gp = gpar(fontsize = 14),
                    column_names_rot = 90, row_names_side = "right",
                    show_heatmap_legend = F,
                    height = unit(20, "cm"), width = unit(20, "cm"),
                    row_title_gp = gpar(fontsize = 14))
Legend <- Legend(col_fun = my_palette, title = "Pearson correlation", title_position = "leftcenter-rot",
                 legend_height = unit(10, "cm"), grid_width = unit(1, "cm"),
                 # increase legend text size
                    labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 20))


cor_plot

svg("./analysis/correlation_plot_features.svg", width = 15, height = 15)
draw(cor_plot, heatmap_legend_side = "right", annotation_legend_side = "left", annotation_legend_list = list(Legend))
dev.off()






