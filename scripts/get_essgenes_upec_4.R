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


# read in the gff file and skip hash lines
gff <- read.delim("./data/reference_sequences/e_coli_K12.gff3", skip = 3, header = FALSE)
gff_upec <- read.delim("./data/reference_sequences/ecoli536_sRNAs.gff3", skip = 3, header = FALSE)

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

# write the essential genes gff file. do not include colnames or rownames and use tab as separator. use write_de
write_delim(essgenes_gff, "./data/reference_sequences/e_coli_k12_essgenes.gff3", delim = "\t", col_names = FALSE)
write_delim(genes_gff_upec, "./data/reference_sequences/upecgenes.gff3", delim = "\t", col_names = FALSE)

# read in the blastn utput file
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
ot_list_upec_raw <- read_delim("./data/upec_offtargets_fulltranscripts_sorted.csv")

# get the on-targets that are also in the essential genes list
ot_list_upec <- ot_list_upec_raw %>% filter(TIR=="TIR")%>% filter(num_mismatch == 0) %>% filter(locus_tag %in% upec_essgenes) %>%
  filter(!trans_coord>0)

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
unique_pna_vector_upec <- paste0(ot_list_upec$K12_locus_tag, ot_list_upec$trans_coord, ot_list_upec$probe_id)

# do same for k12
unique_pna_vector_k12 <- paste0(ot_list_k12$locus_tag, ot_list_k12$start_position-1, ot_list_k12$pna_name)

# get number of on-targets that are mapping the same position in k12 and upec
nr_overlapping_pnas_k12_upec <- sum(unique_pna_vector_upec %in% unique_pna_vector_k12)
nr_pnas_targeting_essgene_in_upec_but_different <-length(unique(ot_list_upec$probe_id))-nr_overlapping_pnas_k12_upec
pnas_targeting_no_essgene <- length(unique_pna_vector_k12)-nr_overlapping_pnas_k12_upec-nr_pnas_targeting_essgene_in_upec_but_different

# create a ggplot showing relative number of on-targets mapping the same position in k12 and upec
# create a data frame with the data
data_plot <- data.frame(group = c("same PNA target in UPEC as in K12", "targets different ess. gene in UPEC", "targets no essential UPEC gene"),
                        y = c(nr_overlapping_pnas_k12_upec, nr_pnas_targeting_essgene_in_upec_but_different, pnas_targeting_no_essgene)) %>%
  mutate(prop = y / sum(y) *100) %>%
    mutate(ypos = cumsum(prop) - 0.5 * prop) %>% mutate(group= factor(group, levels = c("targets no essential UPEC gene",
                                                                                           "targets different ess. gene in UPEC",
                                                                                        "same PNA target in UPEC as in K12")))

# create the plot
piechart <- ggplot(data_plot, aes(x = "", y = prop, fill=group)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y = "number of on-targets", x = "") +
  coord_polar("y", start=0) + theme_void() + #theme(legend.position="none") +
  geom_text(aes(y = ypos, label = y), color = "white", size=6, vjust=c(7,-5,-6.5), hjust=c(1,3.8,2)) +
  scale_fill_manual(values = c("#66cdaa" ,"#b5cde1" ,"#4682b4"))
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
targets_no_ess_gene <- ot_list_k12$pna_name[! ot_list_k12$pna_name %in% c(targets_both_upec_and_k12,targets_diff_ess_gene)]

# get the MICs for the pnas targeting no essential gene in upec
data_targeting_no_ess_gene <- ot_list_k12[ot_list_k12$pna_name %in% targets_no_ess_gene,] %>%
  select(pna_name, MIC_K12, MIC_UPEC)

# get table with pnas targeting no essential gene in upec and their MICs in upec and k12
ot_list_no_ess_gene <- ot_list_k12[ot_list_k12$pna_name %in% targets_no_ess_gene,]
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
  scale_fill_manual(values = c( "#b5cde1", "#4682b4" ,"#66cdaa"))
piechart

# save piechart to pdf
pdf("./analysis/piechart_mic_upec_k12_comparison.pdf", width = 10, height = 10)
piechart
dev.off()

# add expression data of upec to the data_filterd_UPEC table:
# read in the expression data
tpm_upec <- read.table("./data/transcriptomic_expression_K12_upec/upec_ctrl_log_tpm_plus1.csv", header = T,
                       sep = ",", row.names = 1)

# add average expression of the genes targeted by the pnas to the data_filterd_UPEC table
data_filterd_UPEC <- data_filterd_UPEC %>%
  # add average expression of the genes targeted by the pnas to the data_filterd_UPEC table
  mutate(expression_upec = rowMeans(tpm_upec[match(upec_locus_tag, rownames(tpm_upec)),])) %>%
  # remove expression and sec_structure columns
    select(-expression, -sec_structure, -locus_tag)



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

View(data_filterd_UPEC)




# make plot with all MICs and show how many PNAs have that mic in upec.

# create a data frame with the data
data_plot_upec <- data_filterd_UPEC %>%
  group_by(MIC_UPEC) %>%
  summarise(nr_pnas = n())

# create the plot.
plot_mic_upec <- ggplot(data_plot_upec, aes(x = log2(MIC_UPEC), y = nr_pnas)) +
  geom_bar(stat = "identity", fill="steelblue", color="black") +
  scale_x_continuous(breaks = log2(c(1.25,2.5,5,10,20)), labels = c("1.25","2.5","5","10",">10"))+
  labs(y = "number of PNAs", x = "MIC") +
  theme_classic() +
  theme(axis.text = element_text( size = 12),
        axis.title = element_text( size = 15)) +
  #add numbers to the bars
    geom_text(aes(label = nr_pnas), vjust = -0.3, size = 5)
plot_mic_upec

# save plot to pdf
pdf("./analysis/plot_mic_upec_nr_pnas.pdf", width = 6, height = 5.1)
plot_mic_upec
dev.off()

# do same plot for K12
# create a data frame with the data
data_plot_k12 <- data_filterd_UPEC %>%
  group_by(MIC_K12) %>%
  summarise(nr_pnas = n())

# create the plot. make x axis log2 scale
plot_mic_k12 <- ggplot(data_plot_k12, aes(x = log2(MIC_K12), y = nr_pnas)) +
  geom_bar(stat = "identity", fill="steelblue", color="black") +
  scale_x_continuous(breaks = log2(c(1.25,2.5,5,10,20)), labels = c("1.25","2.5","5","10",">10"))+
  labs(y = "number of PNAs", x = "MIC") +
  theme_classic() +
  theme(axis.text = element_text( size = 12),
        axis.title = element_text( size = 15)) +
  #add numbers to the bars
  geom_text(aes(label = nr_pnas), vjust = -0.3, size = 5)
plot_mic_k12

# save plot to pdf
pdf("./analysis/plot_mic_k12_nr_pnas.pdf", width = 6, height = 5.1)
plot_mic_k12
dev.off()

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
        # add mRNA_half_life to the data_filterd_UPEC table by gene name
          mutate(mRNA_half_life = mrna_decay$mRNA_half_life[match(gene_name, mrna_decay$gene)]) %>%
            # make mean from NA or NaN
  mutate(mRNA_half_life = ifelse(is.na(mRNA_half_life), median(na.omit(mrna_decay$mRNA_half_life)), mRNA_half_life))






# write data_filterd_UPEC to file
write.table(data_filterd_UPEC, "./data/pnas_predictors_mic_upec.tsv", sep = "\t", quote = F, row.names = F)

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


  plot <- ggplot(data_filterd_UPEC_plot, aes(x = MIC_UPEC, y = factor, fill=MIC_UPEC)) +
    # add violin plot
    geom_violin(color = "black", alpha=0.7, width = 1) +
    # add scatterplot
      geom_jitter( size = 2, alpha = 0.5, color = "black", width = 0.2, height = 0) +
    # add boxplot
      geom_boxplot(width=0.2, fill = "white", color = "black", alpha=0.7) +
    # add p-values using wilcox.test and stat_compare_means and put them on top of the plot, write as "p-value = 0.0001"
        stat_compare_means(comparisons = comparisons[pvals_logic], method = "wilcox.test", label = "p.signif") +
    labs(y = ylab, x = "MIC") +
    scale_fill_manual(values= viridis(8, direction = -1, option = "inferno")[2:6]) +
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
plot_RF_mean_cellsytems <- violin_plot(data_filterd_UPEC_plot$RF_mean_cellsytems, "relative fitness", "relative fitness")
plot_RF_mean_cellsytems
plot_Tm <- violin_plot(data_filterd_UPEC_plot$Tm, "Tm", "PNA/mRNA melting temperature (°C)")
plot_Tm
plot_ot_tot_2mm <- violin_plot(data_filterd_UPEC_plot$total_off_targets_2mm, "total off targets 2mm", "total off targets 2mm")
plot_ot_tot_2mm
plot_sc_bases <- violin_plot(data_filterd_UPEC_plot$sc_bases, "sc bases", "self complementary bases")
plot_sc_bases
plot_mfe <- violin_plot(data_filterd_UPEC_plot$MFE_UPEC, "mfe", "minimum free energy in TIR")
plot_mfe
plot_expression_upec <- violin_plot(data_filterd_UPEC_plot$expression_upec, "expression upec", "gene expression in UPEC")
plot_expression_upec
plot_tir_off_targets_1mm <- violin_plot(data_filterd_UPEC_plot$tir_off_targets_1mm, "tir off targets 1mm", "TIR off targets 1mm")
plot_tir_off_targets_1mm
plot_gene_length <- violin_plot(data_filterd_UPEC_plot$gene_length, "gene length", "gene length")
plot_gene_length
plot_purine_percentage <- violin_plot(data_filterd_UPEC_plot$purine_percentage, "purine percentage", "purine %")
plot_purine_percentage
plot_molecular_weight <- violin_plot(data_filterd_UPEC_plot$PNA_molecular_weight, "PNA molecular weight", "PNA molecular weight")
plot_molecular_weight
plot_longest_purine_stretch <- violin_plot(data_filterd_UPEC_plot$longest_purine_stretch, "longest purine stretch", "longest purine stretch")
plot_longest_purine_stretch



# save plots in one big pdf using cowplot
library(cowplot)
# save plots in one big pdf using cowplot
pdf("./analysis/pnas_predictors_mic_upec.pdf", width = 18, height = 20)
plot_grid(plot_Tm, plot_sc_bases, plot_RF_mean_cellsytems,  plot_expression_upec,plot_mRNA_half_life,
          plot_tir_off_targets_1mm, plot_ot_tot_2mm, plot_mfe, plot_gene_length, plot_purine_percentage,
          plot_molecular_weight,plot_longest_purine_stretch,
          ncol = 3, nrow = 4, scale = 0.92 , labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J","K","L"),
      label_size = 20)
dev.off()

# create one only for Tm, sc_bases and expression_upec
svg("./analysis/pnas_predictors_mic_upec_2.svg", width = 18, height = 6)
plot_grid(plot_Tm, plot_sc_bases, plot_expression_upec,
          ncol = 3, nrow = 1, scale = 0.92 , labels = c("A", "B", "C"),
          label_size = 20)
dev.off()

# create a simple violin plot which shows distribution of MICs for upec and k12 as x axis and mic as y axis
# create a dataframe with the MICs of upec and k12
df_new_plot <- data.frame(MIC = c(data_filterd_UPEC$MIC_UPEC, data_filterd_UPEC$MIC_K12),
                          strain = c(rep("UPEC", length(data_filterd_UPEC$MIC_UPEC)),
                                     rep("K12", length(data_filterd_UPEC$MIC_K12))))
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
  scale_fill_manual(values = c("#66cdaa", "#b5cde1"))
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


