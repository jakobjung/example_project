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
# specify working directory
setwd("~/Documents/UPEC_K12_essgenes_2023_03")
library(Biostrings)
library(stringr)
library(tidyverse)
library(readxl)
library(phylotools)
library(KEGGREST)

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
idx_kegg[["eco00061"]] <- c(idx_kegg[["eco00061"]], "b1094")

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
  pw_id <- names(idx_kegg)[i]
  # get kegg pathway name
  pw_name <- list_kegg[[pw_id]]
  # get number of genes in pathway
  tot_genes <- length(idx_kegg[[i]])
  #get number of inhibited genes in both UPEC and K12
  both_inhibited_genes <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$upec_inhibition & gene_specific_DF$e_coli_k12_inhibition]])
  # get number of inhibited genes for UPEC but not K12
  inhibited_genes_upec <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$upec_inhibition & !gene_specific_DF$e_coli_k12_inhibition]])
  # same for K12
  inhibited_genes_k12 <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$e_coli_k12_inhibition & !gene_specific_DF$upec_inhibition]])
  c(pw_id, pw_name, tot_genes, both_inhibited_genes, inhibited_genes_upec, inhibited_genes_k12)
})))

# change colnames of df_kegg
colnames(df_kegg) <- c("pw_id", "pw_name", "tot_genes", "both", "UPEC only", "K12 only")

# make gene numbers integers
df_kegg$tot_genes <- as.integer(df_kegg$tot_genes)
df_kegg$both <- as.integer(df_kegg$both)
df_kegg$"UPEC only" <- as.integer(df_kegg$"UPEC only")
df_kegg$"K12 only" <- as.integer(df_kegg$"K12 only")

# get rid of - Escherichia coli K-12 MG1655 in pw_name
df_kegg$pw_name <- gsub(pattern = " - Escherichia coli K-12 MG1655", replacement = "", df_kegg$pw_name)

# create stacked barplot with the pathway name on the x axis and the number of total and inhibited genes for each
# bacteria on the y axis. Use only pathways with more than 10 genes.
# The plot is split into two parts, one for the total number of genes and one for the inhibited genes.
# The total number of genes is shown in grey and the inhibited genes in red.
# The plot is saved as svg file. Order plot by total number of genes.

df_kegg$pw_name <- factor(df_kegg$pw_name, levels = df_kegg$pw_name[order(df_kegg$tot_genes)])

# create stacked barplot
kegg_plot <- ggplot(df_kegg[df_kegg$tot_genes > 4,], aes(x = pw_name)) +
  geom_bar(aes(y = tot_genes), stat = "identity", fill = "lightgrey") +
  geom_bar(aes(y = `K12 only`+ `UPEC only` + both), stat = "identity", fill = "#66cdaa") +
  geom_bar(aes(y = `UPEC only` + both), stat = "identity", fill = "#b5cde1") +
  geom_bar(aes(y = both), stat = "identity", fill = "#4682b4") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
      axis.title = element_text(size = 13)) +
  labs(x = "KEGG pathway", y = "Number of genes") +
  scale_fill_manual(values = c("grey", "red")) +
  guides(fill = guide_legend(title = "Inhibited genes", nrow = 2, byrow = TRUE))

kegg_plot

# save plot as svg file
svg(filename = "./analysis/gene_specific_predictors/KEGG_plot.svg", width = 10, height = 9, pointsize = 12)
print(kegg_plot)
dev.off()


# get number of pathways per gene:
gene_specific_DF$nr_pathways <- sapply(gene_specific_DF$locus_tag, function(x){
  sum(sapply(idx_kegg, function(y) x %in% y))
})


# import transcriptome data and give average value per gene:
transcriptome <- read_xlsx("./data/transcriptomic_expression_K12_upec/palsson_2019_13483_MOESM4_ESM.xlsx",
                           sheet = "Expression Data")[1:3] %>%
  mutate(avg_expression = rowMeans(select(., "control__wt_glc__1", "control__wt_glc__2"), na.rm = TRUE)) %>%
  select(`log-TPM`, avg_expression) %>% add_row(`log-TPM`= "b1085", avg_expression=5.70) %>%
  add_row(`log-TPM`= "b2891", avg_expression=9.53)

# add avg expression to gene_specific_DF
gene_specific_DF$gene_expression <- unlist(transcriptome[transcriptome$`log-TPM` %in% gene_specific_DF$locus_tag ,
                                                  "avg_expression"])




# add nr of pathways to each row of pnas_filtered, using the locus_tag as key
pnas_filtered$nr_pathways <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$nr_pathways))
pnas_filtered$expression <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$gene_expression))

pnas_only_double <- pnas_filtered %>%
  # filter genes that inhibit either E. coli K12 or UPEC, but not both
  filter((upec_inhibition == "TRUE" & e_coli_k12_inhibition == "TRUE") |
           (upec_inhibition == "FALSE" & e_coli_k12_inhibition == "FALSE"))


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


pnas_filtered$nr_genes_operon <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$nr_genes_in_operon)) -1

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

pnas_filtered$downstream_genes_operon <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$nr_downstream_genes_operon))

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
pnas_filtered$essential_genes_downstream <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$essential_genes_downstream))

# do the same thing for upstream genes and ustream essential genes
gene_specific_DF$nr_upstream_genes_operon <- unlist(sapply(gene_specific_DF$gene_name, function(x){
  l <- 0
  ngenes <- unlist(sapply(operon_data$genes_in_operon, function(y){
    if(grepl(x, y)){
      operon <- strsplit(y, ",")[[1]]
      # if yes, count number of commas and add 1
      upstream_genes_nr <- which(operon==x)-1
      return(as.integer(upstream_genes_nr))
    }
  }))
  if(is.integer(ngenes)){
    l <- ngenes
  }
  l
}))

pnas_filtered$upstream_genes_operon <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$nr_upstream_genes_operon))

# check whether there's an essential gene upstream of the genes in the operon.
essgenes <- gene_specific_DF$gene_name

gene_specific_DF$essential_genes_upstream <- unlist(sapply(gene_specific_DF$gene_name, function(x){
  l <- 0
  ngenes <- unlist(sapply(operon_data$genes_in_operon, function(y){
    if(grepl(x, y)){
      operon <- strsplit(y, ",")[[1]]
      return(sum(sapply(operon[1:which(operon==x)], function (x) x %in% essgenes)))
    }
  }))
  if(is.integer(ngenes)){
    l <- ngenes-1
  }
  l
}))

# add essential upstream genes to pnas_filtered
pnas_filtered$essential_genes_upstream <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$essential_genes_upstream))

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
gene_specific_DF$length <- unlist(sapply(gene_specific_DF$locus_tag, function(x) gff_egenes[gff_egenes$locus_tag==x,]$length))

# add gene length column to pnas_filtered
pnas_filtered$gene_length <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gff_egenes$locus_tag==x,]$length))

# add GC-content column to gff using fasta file and gff file start end end to extract fasta record for each gene and get
# GC-content
gff_egenes$GC_content <- unlist(sapply(gff_egenes$locus_tag, function(x){
  # get start and end position for gene
  start <- gff_egenes[gff_egenes$locus_tag==x,4][[1]]
  end <- gff_egenes[gff_egenes$locus_tag==x,5][[1]]
  # extract fasta record for gene
  fasta_record <- as.character(fasta[[1]][start:end])
  # get GC-content
  GC_content <- (str_count(fasta_record, "C") + str_count(fasta_record, "G")) / nchar(fasta_record)
  GC_content
}))

# add GC-content column to gene_specific_DF
gene_specific_DF$gene_GC_content <- unlist(sapply(gene_specific_DF$locus_tag, function(x) gff_egenes[gff_egenes$locus_tag==x,12]))


# add GC-content column to pnas_filtered
pnas_filtered$gene_GC_content <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$gene_GC_content))

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
    sequence <- reverseComplement(fasta[[1]][(end-15):(end+30)])
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
pnas_filtered$sec_structure <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$sec_structure))

# import genes that have an upstrream gene overlapping
overlapping_genes <- read_delim("./data/overlapping_genes_lt.tab", delim = "\t", col_names = c("gene", "locus_tag"))
# get a yes/no column for each gene in gene_specific_DF
gene_specific_DF$overlapping_genes <- unlist(sapply(gene_specific_DF$locus_tag, function(x) x %in% overlapping_genes$locus_tag))
# add overlapping_genes column to pnas_filtered
pnas_filtered$overlapping_genes <- unlist(sapply(pnas_filtered$locus_tag,function(x) gene_specific_DF[gene_specific_DF$locus_tag==x,]$overlapping_genes))


# save gene_specific_DF to file
write_csv(pnas_filtered, "./data/all_predictors.csv")



gene_specific_DF$lowest_mic <- unlist(sapply(gene_specific_DF$locus_tag, function(x) min(pnas_filtered[pnas_filtered$locus_tag==x,]$MIC)))
