



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
library(circlize)

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

  # order df_kegg by tot_genes and then by both
  df_kegg$pw_name <- factor(df_kegg$pw_name, levels = df_kegg$pw_name[order(df_kegg$tot_genes, df_kegg$`UPEC only`)])

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
    theme(axis.text.y = ggtext::element_markdown(size = 12, colour = rev(ifelse(df_kegg[df_kegg$tot_genes > 4,]$pw_name %in% gsub(" - Escherichia coli K-12 MG1655", "", list_kegg[pw_acpp]), "blue", "black"))),
                                                 #colour = rev(ifelse(df_kegg[df_kegg$tot_genes > 4,]$pw_id %in% pw_acpp, "blue", "black"))),
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

# pnas_filtered: take only the ones which have 2 entries per gene_name
pnas_filtered_only_2 <- pnas_filtered %>%
  group_by(gene_name) %>%
  filter(n() == 2)%>% select(gene_name) %>% unique %>% unlist()


# create per gene the lower of MICs in upec
lowest_mic_gene <- pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  group_by(gene_name) %>%
  summarise(lowest_mic_upec = min(MIC_UPEC, na.rm = T)) %>%
  select(gene_name, lowest_mic_upec)

avg_mic_gene <- pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  group_by(gene_name) %>%
  summarise(avg_mic_upec = mean(MIC_UPEC, na.rm = T)) %>%
  select(gene_name, avg_mic_upec)

# get difference between 2 ppnas per gene
diff_mic_gene <- pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  mutate(MIC_UPEC = log2(MIC_UPEC)) %>%
  group_by(gene_name) %>%
  summarise(diff_mic_upec = diff(range(MIC_UPEC, na.rm = T))) %>%
    select(gene_name, diff_mic_upec)

pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  # add a row for average MIC of upec for each gene
  left_join(lowest_mic_gene, by = "gene_name") %>%
  # add a row for difference in MIC of upec for each gene
  left_join(diff_mic_gene, by = "gene_name") %>% filter(lowest_mic_upec<5) %>% select(gene_name, diff_mic_upec, MIC_UPEC)

# now create a plot from pnas_filtered, in which y axis is genes and x axis is MIC of UPEC. So there should be 2
# connedted dots per gene.
mic_plot <- pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  # add a row for average MIC of upec for each gene
    left_join(lowest_mic_gene, by = "gene_name") %>%
  # add avg mic
    left_join(avg_mic_gene, by = "gene_name") %>%
  # add a row for difference in MIC of upec for each gene
    left_join(diff_mic_gene, by = "gene_name") %>%
  mutate(MIC_UPEC = log2(MIC_UPEC)) %>%
  # make difference a factor with increasing levels
  mutate(diff_mic_upec = factor(diff_mic_upec, levels = c(0, 1, 2, 3, 4))) %>%
  # make gene_name a factor and sort levels by lowest MIC and then by avg MIC
    mutate(gene_name = factor(gene_name, levels = lowest_mic_gene$gene_name[order(avg_mic_gene$avg_mic_upec, lowest_mic_gene$lowest_mic_upec,

                                                                                  decreasing = T)])) %>%
 #   mutate(gene_name = factor(gene_name, levels = gene_name[order(MIC_UPEC)])) %>%
  ggplot(aes(x = MIC_UPEC, y = gene_name, color = diff_mic_upec)) +
  # make points, but jitter if there are multiple points at the same x value
    geom_point(size = 0.5) +
    # make points connected by gene_name
    geom_line() +
  # make x axis ticks at log2(1.25, 2.5, 5,10,20) and manually add labels
    scale_x_continuous(breaks = log2(c(1.25, 2.5, 5, 10, 20)), labels = c("1.25", "2.5", "5", "10", ">10")) +
  # add theme
    theme_classic() +
  ylab("Target genes") +
  # add label for legend
    labs(color = "Difference in MIC") +
  # increase size of x axis labels
    theme(axis.text.x = element_text(size = 12),
          axis.title = element_text(size = 15),
    # remove row names
            axis.text.y = element_blank(),
          # add legend to inside of plot upper right
            legend.position = c(0.85, 0.9),
    legend.text = element_text(size = 10)
    )  +
  # change legend ticks to 0, 1, 2, 3, 4
  scale_color_viridis_d(labels = c("0", "2", "4", "8", ">8"),
                        breaks = c("0", "1", "2", "3", "4"))

mic_plot

# save plot as svg file
svg(filename = "./analysis/gene_specific_predictors/MIC_plot.svg", width = 5, height = 10, pointsize = 12)
print(mic_plot + guides(color = guide_legend(override.aes = list(size=3))))
dev.off()


# get difference between 2 ppnas per gene
diff_mic_gene_k12 <- pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  mutate(MIC_K12 = log2(MIC_K12)) %>%
  group_by(gene_name) %>%
  summarise(diff_mic_k12 = diff(range(MIC_K12, na.rm = T))) %>%
  select(gene_name, diff_mic_k12)

lowest_mic_gene_k12 <- pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  group_by(gene_name) %>%
  summarise(lowest_mic_k12 = min(MIC_K12, na.rm = T)) %>%
  select(gene_name, lowest_mic_k12)

avg_mic_gene_k12 <- pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  group_by(gene_name) %>%
  summarise(avg_mic_k12 = mean(MIC_K12, na.rm = T)) %>%
  select(gene_name, avg_mic_k12)

# now create a plot from pnas_filtered, in which y axis is genes and x axis is MIC of K12. So there should be 2
# connedted dots per gene.
mic_plot_k12 <- pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  # add a row for average MIC of upec for each gene
  left_join(lowest_mic_gene, by = "gene_name") %>%
  # add a row for difference in MIC of upec for each gene
  left_join(diff_mic_gene_k12, by = "gene_name") %>%
  mutate(MIC_K12 = log2(MIC_K12)) %>%
  # make difference a factor with increasing levels
  mutate(diff_mic_k12 = factor(diff_mic_k12, levels = c(0, 1, 2, 3, 4))) %>%
  # make gene_name a factor and sort levels by lowest MIC
  mutate(gene_name = factor(gene_name, levels = lowest_mic_gene_k12$gene_name[order(avg_mic_gene$avg_mic_upec,
                                                                                    lowest_mic_gene$lowest_mic_upec,
                                                                                    decreasing = T)])) %>%
  #   mutate(gene_name = factor(gene_name, levels = gene_name[order(MIC_UPEC)])) %>%
  ggplot(aes(x = MIC_K12, y = gene_name, color = diff_mic_k12)) +
  # make points, but jitter if there are multiple points at the same x value
  geom_point(size = 0.5) +
  # make points connected by gene_name
  geom_line() +
  # make x axis ticks at log2(1.25, 2.5, 5,10,20) and manually add labels
  scale_x_continuous(breaks = log2(c(1.25, 2.5, 5, 10, 20)), labels = c("1.25", "2.5", "5", "10", ">10"),
                        limits = log2(c(1.25, 20))) +
  # add theme
  theme_classic() +
  ylab("Target genes") +
  # increase size of x axis labels
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = "none",
        # remove row names
        axis.text.y = element_blank()
    )  +
    # change legend to viridis, keep color scale similar to previous plot
    scale_color_manual(values = c(viridis(5)[1], viridis(5)[2], viridis(5)[3], viridis(5)[4], viridis(5)[5]),
                       breaks = c("0", "1", "2", "3", "4"),
                       labels = c("0", "1", "2", "3", "4"))

mic_plot_k12

# save plot as svg file
svg(filename = "./analysis/gene_specific_predictors/MIC_plot_k12.svg", width = 5, height = 10, pointsize = 12)
print(mic_plot_k12)
dev.off()

# save both mic plots next to each other as svg file using cowplot
library(cowplot)
svg(filename = "./analysis/gene_specific_predictors/MIC_plot_combined_same_order.svg", width = 10, height = 10, pointsize = 12)
plot_grid(mic_plot+ guides(color = guide_legend(override.aes = list(size=3))), mic_plot_k12, ncol = 2)
dev.off()


# do same plots but



mic_plot_k12 <- pnas_filtered %>%
  filter(gene_name %in% pnas_filtered_only_2) %>%
  # add a row for average MIC of upec for each gene
  left_join(lowest_mic_gene, by = "gene_name") %>%
  # add a row for difference in MIC of upec for each gene
  left_join(diff_mic_gene_k12, by = "gene_name") %>%
  mutate(MIC_K12 = log2(MIC_K12)) %>%
  # make difference a factor with increasing levels
  mutate(diff_mic_k12 = factor(diff_mic_k12, levels = c(0, 1, 2, 3, 4))) %>%
  # make gene_name a factor and sort levels by lowest MIC
  mutate(gene_name = factor(gene_name, levels = lowest_mic_gene_k12$gene_name[order(avg_mic_gene_k12$avg_mic_k12,
                                                                                    lowest_mic_gene_k12$lowest_mic_k12,
                                                                                    decreasing = T)])) %>%
  #   mutate(gene_name = factor(gene_name, levels = gene_name[order(MIC_UPEC)])) %>%
  ggplot(aes(x = MIC_K12, y = gene_name, color = diff_mic_k12)) +
  # make points, but jitter if there are multiple points at the same x value
  geom_point(size = 0.5) +
  # make points connected by gene_name
  geom_line() +
  # make x axis ticks at log2(1.25, 2.5, 5,10,20) and manually add labels
  scale_x_continuous(breaks = log2(c(1.25, 2.5, 5, 10, 20)), labels = c("1.25", "2.5", "5", "10", ">10"),
                        limits = log2(c(1.25, 20))) +
  # add theme
  theme_classic() +
  ylab("Target genes") +
  # increase size of x axis labels
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = "none",
        # remove row names
        axis.text.y = element_blank()
    )  +
    # change legend to viridis, keep color scale similar to previous plot
    scale_color_manual(values = c(viridis(5)[1], viridis(5)[2], viridis(5)[3], viridis(5)[4], viridis(5)[5]),
                       breaks = c("0", "1", "2", "3", "4"),
                       labels = c("0", "1", "2", "3", "4"))


library(cowplot)
svg(filename = "./analysis/gene_specific_predictors/MIC_plot_combined_both_sorted.svg", width = 10, height = 10, pointsize = 12)
plot_grid(mic_plot+ guides(color = guide_legend(override.aes = list(size=3))), mic_plot_k12, ncol = 2)
dev.off()






# OK now create a circlize plot which shows for each essential gene the location, the MICs of UPEC and the MICs of K12.

# we first need to generate a df with one gene per row and the MICs of UPEC and K12 in the columns (each have 2 values)
df_gene_wise <- pnas_filtered %>%
  select(gene_name, locus_tag, MIC_UPEC, MIC_K12) %>%
  # make new column PNA. if there are duplicate entries for gene_name, enumerate them (1 and 2)
    group_by(gene_name) %>%
    mutate(PNA = row_number()) %>%
    ungroup() %>%
  # compress df, so that each gene has only one row from MIC_K12 make MIC_K12_1 and MIC_K12_2 if there are 2 entries
  # and the same for MIC_UPEC
    pivot_wider(names_from = PNA, values_from = c(MIC_UPEC, MIC_K12))
  # add nr 1 or 2 to val


# OK now import the ess. gene gff file and extract the start and end position of each gene
# import gff file
gff <- read_delim("./data/reference_sequences/e_coli_K12.gff3", delim = "\t", skip = 3,
                  col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"))

# get only ess. genes
gff_egenes <- gff %>% mutate(locus_tag = gsub(".+;locus_tag=(b\\d+).*", "\\1", attributes)) %>%
  # keep only one of 2 entries for each gene/CDS if it has same start, end position and locus_tag
    group_by(locus_tag) %>% filter(row_number()==1) %>% ungroup() %>%
    # keep only genes and CDS
    filter(feature %in% c("gene", "CDS")) %>%
    # keep only locus tags that are in gene_specific_DF$locus_tag
    filter(locus_tag %in% df_gene_wise$locus_tag) %>%
    # add length column to gff
    mutate(length = end - start + 1)

# add start and end position to df_gene_wise
df_gene_wise <- df_gene_wise %>%
  left_join(gff_egenes %>% select(locus_tag, start, end, length), by = "locus_tag")


# make df with 1 gene and start at 1 and end at 4599600
df_plot <- data.frame(gene_name = "E. coli genome", start = 1, end = 4650000)

gene_locations_upec <- df_gene_wise %>% mutate(gene = "E. coli genome") %>%
    mutate(MIC_UPEC_1 = log2(MIC_UPEC_1), MIC_UPEC_2 = log2(MIC_UPEC_2)) %>%
  select(gene, start, end, MIC_UPEC_1, MIC_UPEC_2) %>% as.data.frame()

# same for K12
gene_locations_k12 <- df_gene_wise %>% mutate(gene = "E. coli genome") %>%
  mutate(MIC_K12_1 = log2(MIC_K12_1), MIC_K12_2 = log2(MIC_K12_2)) %>%
  select(gene, start, end, MIC_K12_1, MIC_K12_2) %>% as.data.frame()



svg("./analysis/circle_plot_mics.svg")
circos.clear()
# initialize circos plot but remove numbers on the outside of the plot
circos.genomicInitialize(df_plot)
text(0, 0, "E. coli genome", cex = 1, italics = TRUE)

# now agg ess. genes to the circle

# circos.initialize(df_plot, xlim = c(0, 1))
# add the essential genes to the plot
#circos.genomicHeatmap(gene_locations_upec, side = "outside")
circos.genomicTrack(gene_locations_upec, panel.fun = function(region, value, ...) {
  #add rect, blue!
    circos.genomicRect(region, value, border = alpha("black", 0.3), ...)
})

circos.genomicTrack(gene_locations_upec, panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, col = alpha("black", 0.5), bg = alpha("#b5cde1", 0.5), cex = 1,  pch = 21)
})
# same for K12
circos.genomicTrack(gene_locations_k12, panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, value, col = alpha("black", 0.5), bg = alpha("#66cdaa", 0.5), cex = 1,  pch = 21)
})
dev.off()
circos.clear()




