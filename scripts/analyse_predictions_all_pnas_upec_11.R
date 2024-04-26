# analyse the predictions that we have made with the upec data.
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# read in the data
all_pnas_all_egenes_upec <- read.table("./data/all_genes_all_pnas/all_pnas_all_egenes_upec_predicted.tsv", header = T, sep = "\t")

# get table of tried pnas
tried_pnas <- read.table("./data/pnas_predictors_mic_upec.tsv", header = T, sep = "\t")

all_pnas_all_egenes_upec$tested <- all_pnas_all_egenes_upec$PNA_sequence %in% tried_pnas$pna_sequence


topgenes <- all_pnas_all_egenes_upec %>%
  # filter genes with highest avg predicted mic
    group_by(K12_genename) %>%
    summarise(avg_predicted_MIC = mean(predicted_MIC)) %>%
    arrange(desc(avg_predicted_MIC)) %>%
    tail(30) %>% select(K12_genename) %>% unlist() %>% as.vector()

#same with min
topgenes_min <- all_pnas_all_egenes_upec %>%
  # filter genes with highest avg predicted mic
    group_by(K12_genename) %>%
    summarise(avg_predicted_MIC = min(predicted_MIC)) %>%
    arrange(desc(avg_predicted_MIC)) %>%
    tail(30) %>% select(K12_genename) %>% unlist() %>% as.vector()

tgenes <- unique(c(topgenes, topgenes_min))

# make a point plot of the predicted mic
plottop <- all_pnas_all_egenes_upec %>% filter(K12_genename %in% tgenes) %>%
  ggplot(aes(x = K12_genename, y = predicted_MIC, color = tested, size= tested)) +
  geom_point(alpha = 0.5) +
  labs(x = "gene", y = "predicted MIC") +
  theme_minimal()+
  # make x labels vertical
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
   scale_y_continuous(breaks = log2(c(1.25, 2.5, 5,10,20)), labels = c("1.25", "2.5", "5","10",">10"), limits = c(0.3, 5))+
  scale_color_manual(values = c("TRUE" = "darkorange", "FALSE" = "purple")) +
  scale_size_manual(values = c("TRUE" = 2, "FALSE" = 1))+
  # make x labels italic
    theme(axis.text.x = element_text(face = "italic"))
  # make y labes at log2(1.25, 2.5, 5, 10, 20)

plotall <-  all_pnas_all_egenes_upec %>%
  ggplot(aes(x = K12_genename, y = predicted_MIC, color = tested, size= tested)) +
  geom_point(alpha = 0.5) +
  labs(x = "gene", y = "predicted MIC") +
  theme_minimal()+
  # make x labels vertical
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
   scale_y_continuous(breaks = log2(c(1.25, 2.5, 5,10,20)), labels = c("1.25", "2.5", "5","10",">10"), limits = c(0.3, 5))+
  scale_color_manual(values = c("TRUE" = "darkorange", "FALSE" = "purple")) +
  scale_size_manual(values = c("TRUE" = 2, "FALSE" = 1))+
  # make x labels italic
    theme(axis.text.x = element_text(face = "italic"))

plottop

svg("./analysis/allgenes_allpnas/top_genes_predicted_mic.svg", height = 4, width = 10)
plottop
dev.off()

svg("./analysis/allgenes_allpnas/all_genes_predicted_mic.svg", height = 4, width = 40)
plotall
dev.off()

plottop_mason <- all_pnas_all_egenes_upec %>% filter(K12_genename %in% tgenes) %>%
  ggplot(aes(x = K12_genename, y = predicted_MIC_mason, color = tested, size= tested)) +
  geom_point(alpha = 0.5) +
  labs(x = "gene", y = "predicted MIC") +
  theme_minimal()+
  # make x labels vertical
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
   scale_y_continuous(breaks = log2(c(1.25, 2.5, 5,10,20)), labels = c("1.25", "2.5", "5","10",">10"), limits = c(0.3, 5))+
  scale_color_manual(values = c("TRUE" = "darkorange", "FALSE" = "purple")) +
  scale_size_manual(values = c("TRUE" = 2, "FALSE" = 1)) +
  # make x labels italic
    theme(axis.text.x = element_text(face = "italic"))
  # make y labes at log2(1.25, 2.5, 5, 10, 20)





# check correlation between predicted mic of mason and other alg:
ggplot(all_pnas_all_egenes_upec, aes(x = predicted_MIC, y = predicted_MIC_mason)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "predicted MIC", y = "mason-adjusted algorithm predicted MIC") +
  theme_minimal()


# get data from mason and nonmason
diagnostics_mason <- read_csv("./data/mason_modeling/df_results_all_metrics_optimized_models.csv")
diagnostics_nonmason <- read_csv("./data/df_results_all_metrics_optimized_models.csv")

# combine these 2 dataframes (add column to indicate mason or nonmason)
diagnostics_mason$algorithm <- "MASON-adjusted"
diagnostics_nonmason$algorithm <- "Raw algorithm"

diagnostics <- rbind(diagnostics_nonmason, diagnostics_mason)

# make algorithm a factor and make mason second
diagnostics$algorithm <- factor(diagnostics$algorithm, levels = c("Raw algorithm", "MASON-adjusted"))

# create a boxplot of within-1-tier accuracy of mason vs non mason algorithms
# make this a function inputting only the metric and the diagnostics dataframe and ylim
plot_metric <- function(metric, diagnostics, ylim, laby){
  diagnostics %>% filter(scoring_metric == metric) %>%
    ggplot(aes(x = algorithm, y = Median_Scores)) +
    geom_boxplot(fill = "steelblue", alpha=0.5) +
    geom_jitter(width = 0.05, alpha = 0.2, color = "black") +
    #geom_jitter(width = 0.01, alpha = 0.5) +
    labs(x = "algorithm", y = metric) +
    theme_bw() + ylim(ylim) +
    ylab(laby) +
    # add stars if theres a significant difference (non-parametric test- wilcoxon)
    stat_compare_means(comparisons = list(c("Raw algorithm", "MASON-adjusted")), method = "wilcox.test",
                       label = "p.signif", size = 3, label.y = ylim[2] - 0.1 * diff(ylim))
}
w1a_plot <- plot_metric("within-1-tier accuracy", diagnostics, c(0.65, 0.85), "within-1-tier accuracy (higher = better)")

ggsave("./analysis/mason_modeling/within-1-tier-accuracy.svg", w1a_plot ,  height = 3.5, width = 3.5)

# use r2
r2_plot <- plot_metric("R2 score", diagnostics, c(0.05, 0.4), "R2 score (higher = better)")
ggsave("./analysis/mason_modeling/r2.svg", r2_plot, height = 3.5, width = 3.5)


# use "Root mean squared error"
rmse_plot <- plot_metric("Root mean squared error", diagnostics, c(0.75, 1), "Root mean squared error (lower = better)")
ggsave("./analysis/mason_modeling/rmse.svg", rmse_plot, height = 3.5, width = 3.5)


# get the shap value data from mason
shap_mason_model <- read_csv("./data/mason_modeling/shap_values_rf_optimized_model.csv")
data_mason_model <- read_csv("./data/mason_modeling/X_selected.csv")

# add _shap to all column names in shap_mason_model
colnames(shap_mason_model) <- paste0(colnames(shap_mason_model), "_shap")

# cbind the shap values to the data
data_mason_model_shap <- as.tibble(cbind(data_mason_model, shap_mason_model))

library(viridis)
# make a scatter plot showing Tm, TM_shap
plot_TM <- data_mason_model_shap %>%
  ggplot(aes(x = Tm, y = -Tm_shap, fill=sc_bases)) +
  geom_point(alpha=0.8, shape=21, color="darkgrey", size=2, stroke=0.2) +
  labs(x = "Tm", y = "Tm SHAP value") +
  theme_bw() +
  # add horizontal line at 0
    geom_hline(yintercept = 0, linetype = "dashed", alpha=0.5) +
  ylab("SHAP value") +
    xlab("predicted PNA-mRNA Tm (°C)")+
  # # add color scale   # adust legend ticks of color legend
  #   scale_color_viridis(discrete=F, breaks = c( 0.3, 0.6),
  #                       labels = c( "0.3", "0.6"), name = "% SC \nbases")+
  # add color scale (0-0.5 and 0.5-1). just use 2 colors
    scale_fill_viridis(name = "self-comp.\n  bases:", breaks= c(0, 1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9),
                        labels = c("0", "","2","","4","","6","","8"), direction = 1, option = "inferno")+
  ylim(-0.95, 0.95)

plot_TM

# save the plot
ggsave("./analysis/mason_modeling/Tm_shap.svg", plot_TM, height = 3.5, width = 5)


# make a plot for sc_bases interacting with Tm
p_sc <- data_mason_model_shap %>%
  ggplot(aes(x = jitter(sc_bases, 0.5), y = -sc_bases_shap)) +
  geom_point(size=2,alpha=0.4, color="black", shape=21, fill="steelblue", stroke=0.3) +
  labs(x = "% SC bases", y = "SC bases SHAP value") +
  theme_bw() +
  # add horizontal line at 0
    geom_hline(yintercept = 0, linetype = "dashed", alpha=0.5) +
  ylab("SHAP value") +
    xlab("self-complementary bases")+
  # add color scale   # adust legend ticks of color legend
    scale_color_viridis(name = "CAI")+
  ylim(-0.95, 0.95) +
  # add x axis ticks and labels manually:
  scale_x_continuous(breaks= c(0,  2/9, 4/9,  6/9,  8/9),
                     labels = c("0", "2","4","6","8"))

p_sc

ggsave("./analysis/mason_modeling/sc_bases_shap.svg", p_sc, height = 3.5, width = 3.5)

# do same for non-mason model
shap_nonmason_model <- read_csv("./data/shap_values_rf_optimized_model.csv")

