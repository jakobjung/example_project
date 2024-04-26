##
##   Title:       Create MIC output
##
##   Author:      Jakob Jung
##
##   Date:        2023-11-10
##
##   Description: This script creates a kegg plot for the genes that are targeted by the PNAs.

# load libraries
library(KEGGREST)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(stringr)
library(readr)
library(tidyr)

# load data
pnas_filtered <- read_tsv("./data/pnas_predictors_mic_upec.tsv")

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
  select(locus_tag, gene_name, MIC_K12, MIC_UPEC) %>%
  # remove duplicates of locustags and keep TRUE if one of e. coli or upec is inhibited, same for gene_inhibition.
  group_by(locus_tag) %>%
  # get smallest MIC value for each locus tag
  summarise(
    gene_name = unique(gene_name),
    MIC_K12 = min(MIC_K12),
    MIC_UPEC = min(MIC_UPEC)
  )

# make function for kegg plot:
make_kegg_plot <- function (threshold){
  # create df for kegg_pws
  df_kegg <- as_tibble(t(sapply(seq_along(idx_kegg), function (i) {
    # get kegg pathway id
    pw_id <- names(idx_kegg)[i]
    # get kegg pathway name
    pw_name <- list_kegg[[pw_id]]
    # get number of genes in pathway
    tot_genes <- length(idx_kegg[[i]])
    #get number of inhibited genes in both UPEC and K12
    both_inhibited_genes <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_UPEC <= threshold & gene_specific_DF$MIC_K12 <= threshold]])
    # get number of inhibited genes for UPEC but not K12
    inhibited_genes_upec <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_UPEC <= threshold & gene_specific_DF$MIC_K12 > threshold]])
    # same for K12
    inhibited_genes_k12 <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_UPEC > threshold & gene_specific_DF$MIC_K12 <= threshold]])
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
    geom_bar(aes(y = tot_genes, fill = "neither"), stat = "identity") +
    geom_bar(aes(y = `K12 only`+ `UPEC only` + both, fill = "only K12"), stat = "identity") +
    geom_bar(aes(y = `UPEC only` + both, fill = "only UPEC"), stat = "identity") +
    geom_bar(aes(y = both, fill = "UPEC & K12"), stat = "identity") +
    scale_fill_manual(values = c("only K12" = "#66cdaa", "only UPEC" = "#b5cde1", "UPEC & K12" = "#4682b4",
                                  "neither" = "lightgrey")) +
    coord_flip() +
    theme_classic() +
    theme(axis.text.y = ggtext::element_markdown(size = 12, colour = rev(ifelse(df_kegg[df_kegg$tot_genes > 4,]$pw_id %in% pw_acpp,
                                                                "blue", "black"))),
          axis.text.x = element_text(size = 12),
          axis.title = element_text(size = 13),
    #put legend in lower right part of figure
            legend.position = c(0.85, 0.2),
          ) +
    labs(x = "KEGG pathway", y = "Number of genes", fill = paste0("# genes with \nMIC < = ", threshold, " in..."))

  kegg_plot
}


k_plot_10 <- make_kegg_plot(10)
k_plot_10

# save plot as svg file
svg(filename = "./analysis/gene_specific_predictors/KEGG_plot_th_10.svg", width = 10, height = 9, pointsize = 12)
print(k_plot_10)
dev.off()

k_plot_5 <- make_kegg_plot(5)

# save plot as svg file
svg(filename = "./analysis/gene_specific_predictors/KEGG_plot_th_5.svg", width = 10, height = 9, pointsize = 12)
print(k_plot_5)
dev.off()

k_plot_2_5 <- make_kegg_plot(2.5)




df_kegg <- as_tibble(t(sapply(seq_along(idx_kegg), function (i) {
    # get kegg pathway id
    pw_id <- names(idx_kegg)[i]
    # get kegg pathway name
    pw_name <- list_kegg[[pw_id]]
    #get number of genes wit MIC=1.25 in k12
    genes_1_25_K12 <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_K12 == 1.25]])
    #get number of genes wit MIC=1.25 in upec
    genes_1_25_UPEC <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_UPEC == 1.25]])
    #get number of genes wit MIC=2.5 in k12
    genes_2_5_K12 <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_K12 == 2.5]])
    #get number of genes wit MIC=2.5 in upec
    genes_2_5_UPEC <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_UPEC == 2.5]])
    #get number of genes wit MIC=5 in k12
    genes_5_K12 <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_K12 == 5]])
    #get number of genes wit MIC=5 in upec
    genes_5_UPEC <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_UPEC == 5]])
    #get number of genes wit MIC=10 in k12
    genes_10_K12 <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_K12 == 10]])
    #get number of genes wit MIC=10 in upec
    genes_10_UPEC <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_UPEC == 10]])
    #get number of genes wit MIC=20 in k12
    genes_20_K12 <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_K12 == 20]])
    #get number of genes wit MIC=20 in upec
    genes_20_UPEC <- length(idx_kegg[[i]][idx_kegg[[i]] %in% gene_specific_DF$locus_tag[gene_specific_DF$MIC_UPEC == 20]])
    #put out all the numbers
    c(pw_id, pw_name, genes_1_25_K12, genes_1_25_UPEC, genes_2_5_K12, genes_2_5_UPEC, genes_5_K12, genes_5_UPEC,
      genes_10_K12, genes_10_UPEC, genes_20_K12, genes_20_UPEC)
})))

# change colnames of df_kegg
colnames(df_kegg) <- c("pw_id", "pw_name", "genes_1_25_K12", "genes_1_25_UPEC", "genes_2_5_K12", "genes_2_5_UPEC",
                       "genes_5_K12", "genes_5_UPEC", "genes_10_K12", "genes_10_UPEC", "genes_20_K12", "genes_20_UPEC")

df_kegg_plot <- df_kegg %>% pivot_longer(cols = c("genes_1_25_K12", "genes_1_25_UPEC", "genes_2_5_K12", "genes_2_5_UPEC",
                                             "genes_5_K12", "genes_5_UPEC", "genes_10_K12", "genes_10_UPEC",
                                             "genes_20_K12", "genes_20_UPEC"), names_to = "MIC_str", values_to = "number") %>%
    mutate(strain = ifelse(str_detect(MIC_str, "K12"), "K12", "UPEC")) %>%
    mutate(MIC = as.numeric(gsub("[^_]+_(\\d+)_.*", "\\1", MIC_str))) %>%
    mutate(MIC = ifelse(MIC == 1, 1.25, ifelse(MIC == 2, 2.5, MIC))) %>%
    mutate(pw_name = str_remove(pw_name, " - Escherichia coli K-12 MG1655")) %>%
    select(-MIC_str) %>%
    mutate(number = as.integer(number))



df_kegg_hm_UPEC <- df_kegg[,c(4,6,8,10,12)] %>%
  # make it as integers
    mutate_all(as.integer) %>%
  # divide all values by the respective rowsum
    mutate_all(function(x) x/rowSums(.[,])) %>%
  # make 0 out of NA
    replace(is.na(.), 0)
colnames(df_kegg_hm_UPEC) <- c("1.25", "2.5", "5", "10", ">10")
df_kegg_hm_UPEC <- as.matrix(df_kegg_hm_UPEC)
rownames(df_kegg_hm_UPEC) <- gsub(" - Escherichia coli K-12 MG1655", "", df_kegg$pw_name)

# do same for K12
df_kegg_hm_K12 <- df_kegg[,c(3,5,7,9,11)] %>%
  # make it as integers
  mutate_all(as.integer) %>%
  # divide all values by the respective rowsum
  mutate_all(function(x) x/rowSums(.[,])) %>%
  # make 0 out of NA
  replace(is.na(.), 0)
colnames(df_kegg_hm_K12) <- c("1.25", "2.5", "5", "10", ">10")
df_kegg_hm_K12 <- as.matrix(df_kegg_hm_K12)
rownames(df_kegg_hm_K12) <- gsub(" - Escherichia coli K-12 MG1655", "", df_kegg$pw_name)

# get total nr of genes per pw (rowsums of df_kegg[,c(4,6,8,10,12)])
tot_genes <- rowSums(sapply(df_kegg[,c(4,6,8,10,12)], as.integer))

# get rid of pathways with <5 genes
df_kegg_hm_UPEC <- df_kegg_hm_UPEC[tot_genes > 4,]
df_kegg_hm_K12 <- df_kegg_hm_K12[tot_genes > 4,]
tot_genes <- tot_genes[tot_genes > 4]

# order both matrices by tot_genes (reverse order)
df_kegg_hm_UPEC <- df_kegg_hm_UPEC[order(tot_genes, decreasing = T),]
df_kegg_hm_K12 <- df_kegg_hm_K12[order(tot_genes, decreasing = T),]
tot_genes <- tot_genes[order(tot_genes, decreasing = T)]

# create a heatmap of the matrix:
library(ComplexHeatmap)
library(viridis)
library(circlize)

ht4 <- Heatmap(df_kegg_hm_UPEC, name = "% of genes",
               col = viridis(1000, option = "viridis"),
               cluster_rows = F, cluster_columns = F, show_heatmap_legend = F,
               border = TRUE,
               # put row names on the left
                row_names_side = "left",
               #  create a coumn title
                column_title = "UPEC",
               # increase size of column title
                column_title_gp = gpar(fontsize = 20),
               # make row names bigger and color them blue if they are in the acpp pathway
                row_names_gp = gpar(fontsize = 12, col = ifelse(rownames(df_kegg_hm_UPEC) %in% gsub(" - Escherichia coli K-12 MG1655", "", list_kegg[pw_acpp]), "blue", "black") ),
               # make column and row names bigger
                column_names_gp = gpar(fontsize = 15),
               height = unit(20, "cm"), width = unit(5, "cm"),
               column_names_rot = 45)

ht5 <- Heatmap(df_kegg_hm_K12, name = "% of genes",
               col = viridis(1000, option = "viridis"),
               cluster_rows = F, cluster_columns = F, show_heatmap_legend = F,
               border = TRUE,
               # put row names on the left
               row_names_side = "left",
               #  create a coumn title
                column_title = "K12",
               # increase size of column title
                column_title_gp = gpar(fontsize = 20),
               # make column and row names bigger
               row_names_gp = gpar(fontsize = 12),
               column_names_gp = gpar(fontsize = 15),
               height = unit(20, "cm"), width = unit(5, "cm"),
               column_names_rot = 45,
               right_annotation = rowAnnotation("# genes in\npathway" = anno_barplot(tot_genes), width = unit(3, "cm"),
                                                annotation_name_gp = gpar(fontsize = 20), gp = gpar(fontsize = 20)))

# combine both heatmaps
ht_list <- ht4 + ht5

col_fun = colorRamp2(seq(0,1, length=100), viridis(100))
# create a legend for the heatmap
lgd = Legend(col_fun = col_fun, title = "% of genes with respective MIC", #direction = "horizontal",
             title_gp = gpar(fontsize = 20), labels = c("0%", "50%","100%"), legend_height = unit(20, "cm"),
                legend_width = unit(3, "cm"),
                grid_width = unit(1, "cm"),
             labels_gp = gpar(fontsize = 20),
             at = c(0,0.5, 1), border = "black",
             title_position = "leftcenter-rot")

# save heatmap as svg file
svg(filename = "./analysis/gene_specific_predictors/KEGG_heatmap_UPEC.svg", width = 15, height = 10, pointsize = 12)
draw(ht_list)
draw(lgd, x = unit(31, "cm"), y = unit(3.4, "cm"), just = c("left", "bottom"))
dev.off()











