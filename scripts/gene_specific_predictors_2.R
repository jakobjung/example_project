##
##   Title:       PNA- essential genes screen - calculate gene-specific predictors
##
##   Author:      Jakob Jung
##
##   Date:        2022-03-14
##
##   Description: Calculate gene-specific features such as kegg pathways, gene expression, operon membership, etc.
##

# load libraries
source("./scripts/pna_specific_predictors.R")
library(Biostrings)
library(stringr)

# load data
pnas_filtered <- read_delim("data/pnas_pna_specific_features.tsv")


# get kegg pathways for each gene name
library(KEGGREST)

# get link and list to get kegg info:
link_kegg <- keggLink("pathway", "eco")
list_kegg <- keggList("pathway", "eco")
kegg_pw_ids <- names(list_kegg)

# remove eco: from link_kegg names:
names(link_kegg) <- gsub(pattern = "eco:", replacement = "", names(link_kegg))

# remove kegg links not found in our dataset:
link_kegg <- link_kegg[names(link_kegg) %in% pnas_filtered$locus_tag] #remove genes not in data

# get kegg ids with genes associated with them:
idx_kegg <- sapply(kegg_pw_ids, function(x){
  x <- unique(names(link_kegg[grepl(x, link_kegg)])) # choose all genes, except duplucates
})

# for all entries of pnas_filtered$locus_tag, get all pathways associated with them, so for each ocus tag go through
# idx_kegg list and if the locus tag is found, add the pathway id to the list.
kegg_pathways <- sapply(unique(pnas_filtered$locus_tag), function(x){
  pws <- c()
  for (i in seq_along(idx_kegg)){
    if (x %in% idx_kegg[[i]]){
      pws <- c(pws, kegg_pw_ids[i])
    }
  }
  unique(pws)
})
kegg_pathways$b1094 <- c(kegg_pathways$b1094, "eco00071")

# remove empty entries in idx_kegg
idx_kegg <- idx_kegg[!sapply(idx_kegg, function(x)  length(x)==0)]

# find all pathways that "b1094" is part of
pw_acpp <- kegg_pathways[which(names(kegg_pathways) == "b1094")][[1]]
# get their kegg terms
list_kegg[pw_acpp]

# compress data frame to one locus tag
gene_specific_DF <- pnas_filtered %>%
  select(locus_tag, gene_name, e_coli_k12_inhibition, upec_inhibition, inhibits_either) %>%
  # remove duplicates of locustags and keep TRUE if one of e. coli or upec is inhibited, same for gene_inhibition.
  group_by(locus_tag) %>%
  summarise(
    gene_name = unique(gene_name),
    e_coli_k12_inhibition = any(e_coli_k12_inhibition),
    upec_inhibition = any(upec_inhibition),
    inhibits_either = any(inhibits_either)
  )

# create df for kegg_pws
df_kegg <- as_tibble(t(sapply(seq_along(idx_kegg), function (i) {
  # get kegg pathway id
  pw_id <- names(list_kegg)[i]
  # get kegg pathway name
  pw_name <- list_kegg[i]
  # get number of genes in pathway
  tot_genes <- length(idx_kegg[[i]])
  # get number of inhibited genes in pathway
  tot_inhibited_genes <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$inhibits_either]])
  # get number of inhibited genes for UPEC
  inhibited_genes_upec <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$upec_inhibition]])
  # same for K12
  inhibited_genes_k12 <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$e_coli_k12_inhibition]])
  c(pw_id, pw_name, tot_genes, tot_inhibited_genes, inhibited_genes_upec, inhibited_genes_k12)
})))

# change colnames of df_kegg
colnames(df_kegg) <- c("pw_id", "pw_name", "tot_genes", "tot_inhibited_genes", "inhibited_genes_upec", "inhibited_genes_k12")

# make gene numbers integers
df_kegg$tot_genes <- as.integer(df_kegg$tot_genes)
df_kegg$tot_inhibited_genes <- as.integer(df_kegg$tot_inhibited_genes)
df_kegg$inhibited_genes_upec <- as.integer(df_kegg$inhibited_genes_upec)
df_kegg$inhibited_genes_k12 <- as.integer(df_kegg$inhibited_genes_k12)

# get rid of - Escherichia coli K-12 MG1655 in pw_name
df_kegg$pw_name <- gsub(pattern = " - Escherichia coli K-12 MG1655", replacement = "", df_kegg$pw_name)

# create stacked barplot with the pathway name on the x axis and the number of total and inhibited genes for each
# bacteria on the y axis. Use only pathways with more than 10 genes.
# The plot is split into two parts, one for the total number of genes and one for the inhibited genes.
# The total number of genes is shown in grey and the inhibited genes in red.
# The plot is saved as svg file. Order plot by total number of genes.

df_kegg$pw_name <- factor(df_kegg$pw_name, levels = df_kegg$pw_name[order(df_kegg$tot_genes)])

# create stacked barplot
kegg_plot <- df_kegg %>%
  filter(tot_genes > 5) %>%
  ggplot(aes(x = !!sym(x_axis), fill = inhibition_of)) +
        geom_bar() +
      # order the legend by the order of the levels in the factor
        scale_fill_manual(values = c("K12 only" = "#66cdaa", "UPEC only" = "#b5cde1", "both" = "#4682b4",
                                     "neither" = "#d2d2d2"))
  theme_classic() +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold")
  ) +
  labs(x = "KEGG pathway", y = "Number of genes", fill = "Inhibition of:  ", size=12)+
  scale_fill_manual(values = c("UPEC" = "steelblue", "E. coli K12" = "lightblue", "total genes" = "grey"),
                      limits = c("UPEC", "E. coli K12", "total genes"))+
  # add alpha to e_coli_k12_inhibition
  scale_alpha_manual(values = c("E. coli K12" = 0.7, "UPEC" = 1),  guide = "none") +
  coord_flip()

kegg_plot

# save plot as svg file
svg(filename = "./analysis/gene_specific_predictors/KEGG_plot.svg", width = 10, height = 10, pointsize = 12)
print(kegg_plot)
dev.off()

# get number of pathways per gene:
gene_specific_DF$nr_pathways <- sapply(gene_specific_DF$locus_tag, function(x){
  sum(sapply(idx_kegg, function(y) x %in% y))
})


# import transcriptome data and give average value per gene:
transcriptome <- read_xlsx("./data/transcriptomic_expression_K12/palsson_2019_13483_MOESM4_ESM.xlsx",
                           sheet = "Expression Data")[1:3] %>%
  mutate(avg_expression = rowMeans(select(., "control__wt_glc__1", "control__wt_glc__2"), na.rm = TRUE)) %>%
  select(`log-TPM`, avg_expression) %>% add_row(`log-TPM`= "b1085", avg_expression=5.70) %>%
  add_row(`log-TPM`= "b2891", avg_expression=9.53)

# add avg expression to gene_specific_DF
gene_specific_DF$gene_expression <- unlist(transcriptome[transcriptome$`log-TPM` %in% gene_specific_DF$locus_tag ,
                                                  "avg_expression"])




# add nr of pathways to each row of pnas_filtered, using the locus_tag as key
pnas_filtered$nr_pathways <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,6]))
pnas_filtered$expression <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,7]))

pnas_only_double <- pnas_filtered %>%
  # filter genes that inhibit either E. coli K12 or UPEC, but not both
  filter((upec_inhibition == "TRUE" & e_coli_k12_inhibition == "TRUE") |
           (upec_inhibition == "FALSE" & e_coli_k12_inhibition == "FALSE"))

# make same plot for pnas_only_double df
boxplot_expression_onlydouble <- plot_function("expression", pnas_only_double) + ylab("gene expression (in log TPM)")
ggsave( filename = "./analysis/gene_specific_predictors/only_double_plots/expression_plot_onlydouble.svg", plot = boxplot_expression_onlydouble, width = 7, height = 5)

boxplot_expression <- plot_function("expression", pnas_filtered) + ylab("gene expression (in log TPM)")+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
ggsave( filename = "./analysis/gene_specific_predictors/expression_plot.svg", plot = boxplot_expression, width = 7, height = 5)


nr_pw_plot <- startpos_plot_function("nr_pathways")
expr_plot <- plot_function("expression", pnas_filtered)+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))

# save plot as svg file using ggsave
ggsave( filename = "./analysis/gene_specific_predictors/pathway_nr_plot.svg", plot = nr_pw_plot, width = 8, height = 6)


# add operon information to gene_specific_DF
operon_data <- read_delim("./data/operon_annotation/OperonSet.txt", delim = "\t", skip = 38)

# go through each gene_name and check if it is in the operon_data in "genes_in_operon" column. if yes, get number of
# genes in operon by counting the number of commas in the string and add 1. if not, add 0.
gene_specific_DF$nr_genes_in_operon <- unlist(sapply(gene_specific_DF$gene_name, function(x){
  l = 0
  ngenes <- unlist(sapply(operon_data$genes_in_operon, function(y){
    if(grepl(x, y)){
      # if yes, count number of commas and add 1
      return(length(strsplit(y, ",")[[1]]) + 1)
    }
  }))
  if(is.double(ngenes)){
    l <- ngenes
  } else {
    l <- 2
  }
  l
}))


pnas_filtered$nr_genes_operon <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,8])) -1

# create plot
nr_genes_in_operon_plot <- startpos_plot_function("nr_genes_operon")
nr_genes_in_operon_plot

# save plot as svg file using ggsave
ggsave( filename = "./analysis/gene_specific_predictors/nr_genes_in_operon_plot.svg", plot = nr_genes_in_operon_plot, width = 8)

# get the number of downstream genes for each operon, i.e. the nimber of commas in the "genes_in_operon" column that
# are followed by the matching gene_name. if no match, add 0.
gene_specific_DF$nr_downstream_genes_operon <- unlist(sapply(gene_specific_DF$gene_name, function(x){
  l <- 0
  ngenes <- unlist(sapply(operon_data$genes_in_operon, function(y){
    if(grepl(x, y)){
      operon <- strsplit(y, ",")[[1]]
      # if yes, count number of commas and add 1
      return(length(operon) - which(operon==x))
    }
  }))
  if(is.integer(ngenes)){
    l <- ngenes
  }
  l
}))

pnas_filtered$downstream_genes_operon <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,9]))
nr_ds_in_operon_plot <- startpos_plot_function("downstream_genes_operon")
nr_ds_in_operon_plot

# save plot as svg file using ggsave
ggsave( filename = "./analysis/gene_specific_predictors/nr_ds_in_operon_plot.svg", plot = nr_ds_in_operon_plot, width = 8)

# check whether there's an essential gene downstream of the genes in the operon.
essgenes <- gene_specific_DF$gene_name

gene_specific_DF$essential_genes_downstream <- unlist(sapply(gene_specific_DF$gene_name, function(x){
  l <- 0
  ngenes <- unlist(sapply(operon_data$genes_in_operon, function(y){
    if(grepl(x, y)){
      operon <- strsplit(y, ",")[[1]]
      return(sum(sapply(operon[which(operon==x):length(operon)], function (x) x %in% essgenes)))
    }
  }))
  if(is.integer(ngenes)){
    l <- ngenes-1
  }
  l
}))

# add essential downstream genes to pnas_filtered
pnas_filtered$essential_genes_downstream <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,10]))

# create plot
ess_ds_in_operon_plot <- startpos_plot_function("essential_genes_downstream")
ess_ds_in_operon_plot

# save plot as svg file using ggsave
ggsave( filename = "./analysis/gene_specific_predictors/ess_ds_in_operon_plot.svg", plot = ess_ds_in_operon_plot, width = 8)


# now I want to get the GC-content and the gene length for all genes.
# I will use the fasta file and GFF-file for this
# read in fasta file
fasta <- readDNAStringSet("./data/reference_sequences/e_coli_K12.fasta")

# read in GFF file
gff <- read_delim("./data/reference_sequences/e_coli_K12.gff3", delim = "\t", skip = 3,
                  col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"))

# add locus_tag colmn to gff
gff_egenes <- gff %>% mutate(locus_tag = gsub(".+;locus_tag=(b\\d+).*", "\\1", attributes)) %>%
  # keep only one of 2 entries for each gene/CDS if it has same start, end position and locus_tag
    group_by(locus_tag) %>% filter(row_number()==1) %>% ungroup() %>%
    # keep only genes and CDS
    filter(feature %in% c("gene", "CDS")) %>%
    # keep only locus tags that are in gene_specific_DF$locus_tag
    filter(locus_tag %in% gene_specific_DF$locus_tag) %>%
    # add length column to gff
    mutate(length = end - start + 1)

# add length column to gene_specific_DF
gene_specific_DF$length <- unlist(sapply(gene_specific_DF$locus_tag, function(x) gff_egenes[gff_egenes$locus_tag==x,11]))

# add gene length column to pnas_filtered
pnas_filtered$gene_length <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gff_egenes$locus_tag==x,11]))

# add GC-content column to gff using fasta file and gff file start end end to extract fasta record for each gene and get
# GC-content
gff_egenes$GC_content <- unlist(sapply(gff_egenes$locus_tag, function(x){
  # get start and end position for gene
  start <- gff_egenes[gff_egenes$locus_tag==x,4][[1]]
  end <- gff_egenes[gff_egenes$locus_tag==x,5][[1]]
  # extract fasta record for gene
  fasta_record <- fasta[[1]][start:end]
  # get GC-content
  GC_content <- (str_count(fasta_record, "C") + str_count(fasta_record, "G")) / length(fasta_record)
  GC_content
}))

# add GC-content column to gene_specific_DF
gene_specific_DF$gene_GC_content <- unlist(sapply(gene_specific_DF$locus_tag, function(x) gff_egenes[gff_egenes$locus_tag==x,12]))

# create plot
gene_length_plot <- plot_function("gene_length", pnas_filtered)+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
plot_function("gene_length", pnas_filtered)

# save plot as svg file using ggsave
ggsave( filename = "./analysis/gene_specific_predictors/gene_length_plot.svg", plot = gene_length_plot, width = 7, height = 5)

pnas_only_double <- pnas_filtered %>%
  # filter genes that inhibit either E. coli K12 or UPEC, but not both
  filter((upec_inhibition == "TRUE" & e_coli_k12_inhibition == "TRUE") |
           (upec_inhibition == "FALSE" & e_coli_k12_inhibition == "FALSE"))
# do same plt for  pnas_only_double
gene_length_plot <- plot_function("gene_length", pnas_only_double)
plot_function("gene_length", pnas_only_double)
ggsave( filename = "./analysis/gene_specific_predictors/only_double_plots/gene_length_onlydouble.svg", plot = gene_length_plot, width = 7, height = 5)

# add GC-content column to pnas_filtered
pnas_filtered$gene_GC_content <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,12]))

# for each gene in gff, extract the sequence from -30 to +15 nt relative to the start codon. Take into account the strand.
# I will use the fasta file and GFF-file for this

# add sequence column to gff
gff_egenes$tir_sequence <- unlist(sapply(gff_egenes$locus_tag, function(x){
  # get start and end position for gene
  start <- gff_egenes[gff_egenes$locus_tag==x,4][[1]]
  end <- gff_egenes[gff_egenes$locus_tag==x,5][[1]]
  # get strand
  strand <- gff_egenes[gff_egenes$locus_tag==x,7][[1]]
  # get sequence
  if(strand == "+"){
    sequence <- fasta[[1]][(start-30):(start+15)]
  } else {
    sequence <- complement(rev(fasta[[1]][(end-15):(end+30)]))
  }
  sequence
}))

# save gff_egenes$sequence to fasta file with locus tag as header:
writeXStringSet(DNAStringSet(gff_egenes$tir_sequence), "./data/sec_structure_rnafold/tirs_K12.fasta", format = "fasta")

# after running RNAfold, we download the ./data/sec_structure_rnafold/sec_structure.tsv file and add it to
# gene_specific_DF
sec_structure <- read_delim("./data/sec_structure_rnafold/sec_structure.tsv", delim = "\t",
                            col_names = c("locus_tag","seq", "binding", "delta_G"))

# add sec_structure to gene_specific_DF
gene_specific_DF$sec_structure <- sec_structure$delta_G

# add sec_structure to pnas_filtered
pnas_filtered$sec_structure <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,13]))

# create plot
sec_structure_plot <- plot_function("sec_structure", pnas_filtered) + ylab("free energy of mRNA in TIR (ΔG)")+
    scale_fill_manual(values = c("TRUE" = "powderblue", "FALSE" = "#d2d2d2"))
sec_structure_plot
# save plot as svg file using ggsave
ggsave( filename = "./analysis/gene_specific_predictors/sec_structure_plot.svg", plot = sec_structure_plot, width = 7, height=5)

# get PNAs that inhibit both, vs PNAs that inhibit none
pnas_only_double <- pnas_filtered %>%
  # filter genes that inhibit either E. coli K12 or UPEC, but not both
  filter((upec_inhibition == "TRUE" & e_coli_k12_inhibition == "TRUE") |
           (upec_inhibition == "FALSE" & e_coli_k12_inhibition == "FALSE"))

# do same with pnas_only_double data frame
sec_structure_plot_double <- plot_function("sec_structure", pnas_only_double) + ylab("free energy of mRNA in TIR (ΔG)")
sec_structure_plot_double
# save plot as svg file using ggsave
ggsave( filename = "./analysis/gene_specific_predictors/only_double_plots/sec_structure_plot_double.svg", plot = sec_structure_plot_double, width = 7, height=5)

# save gene_specific_DF to file
write_csv(pnas_filtered, "./data/all_predictors.csv")


# Additional analysis of ess. genes and checking wchich are also in P.

# get gene name from gff_egenes
gff_egenes$gene_name <- unlist(sapply(gff_egenes$locus_tag, function(x){
  # get gene name
  gene_name <-  gsub( ".*Name=([^; ]+).*",  "\\1", gff_egenes[gff_egenes$locus_tag==x,9][[1]])
  gene_name
}))

# import gff from P. laumondii
gff_p_laumondii <- read_delim("../paramita/data/P_laumandii.gff3", delim = "\t", skip = 6,
                  col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"))
gff_p_laumondii$gene_name <- gsub( ".*Name=([^; ]+).*",  "\\1", gff_p_laumondii$attributes)

# import eggnogg annotation:
eggnogg <- read_delim("../paramita/data/p_laumanii.annotations.tsv", delim = "\t")
names_lt_eggnogg <- eggnogg %>% mutate(gene_name=gsub("(.+);.+", "\\1", query)) %>%
  # get lt:
  mutate(locus_tag=gsub(".+;(.+)", "\\1", query)) %>%
  select(gene_name, locus_tag, Preferred_name)


genes_in_both <- gff_egenes$gene_name[(gff_egenes$gene_name %in% names_lt_eggnogg$gene_name)|(gff_egenes$gene_name %in% names_lt_eggnogg$Preferred_name) ]
df_gboth <- tibble(gene_name = genes_in_both, locus_tag_ecoli=names(genes_in_both),
                   locus_tag_plaumondii=names_lt_eggnogg$locus_tag[names_lt_eggnogg$gene_name %in% genes_in_both | names_lt_eggnogg$Preferred_name %in% genes_in_both])

lt_lauman_both <- names_lt_eggnogg[(names_lt_eggnogg$gene_name %in% gff_egenes$gene_name)|(names_lt_eggnogg$Preferred_name %in% gff_egenes$gene_name),]
lt_lauman_both$gname <- ifelse(nchar(lt_lauman_both$gene_name) > 4, lt_lauman_both$Preferred_name, lt_lauman_both$gene_name)

mason_res_laumanii <- read_csv("../mason_commandline/data/paramita_screen/outputs/offtargets_startregions_sorted.csv")

mason_res_laumanii_ess <- mason_res_laumanii[mason_res_laumanii$locus_tag %in% lt_lauman_both$locus_tag,] %>% filter(trans_coord<1)

# get locus tags and number of rows in which they appear. and create a plot visualizing the distribution
lt_counts <- table(mason_res_laumanii_ess$locus_tag)
lt_counts <- data.frame(lt_counts)
names(lt_counts) <- c("locus_tag", "count")
lt_counts <- lt_counts[order(lt_counts$count, decreasing = TRUE),]

# plot
lt_counts_plot <- ggplot(lt_counts, aes(x = count)) +
  geom_histogram(binwidth = 1) +
  labs(x = "Number of off-targets", y = "Number of genes") +
  theme_bw() +
  theme(axis.text.x = element_text( hjust = 1, vjust = 0.5))
lt_counts_plot

x <- lt_lauman_both$gname
names(x) <- lt_lauman_both$locus_tag

laumonii_pnas$gene_name <- x[laumonii_pnas$locus_tag]

laumonii_pnas <- mason_res_laumanii_ess %>% select(probe_id, locus_tag, gene_name, trans_coord, strand)
View(laumonii_pnas)

laumonii_pnas$trans_coord <- laumonii_pnas$trans_coord + 1

laumonii_pnas$PNA_Sequence <- pnas_filtered$pna_sequence[match(laumonii_pnas$PNA_ID, pnas_filtered$pna_name)]


laumonii_pnas$PNA_ID <- gsub("_", "-", laumonii_pnas$PNA_ID)

names(laumonii_pnas) <- c("PNA_ID", "locus_tag_P.laumonii", "gene_name", "binding_start_from_CDS_start", "strand", "PNA_Sequence")
write.xlsx(laumonii_pnas, "../paramita/data/laumonii_pnas.xlsx")