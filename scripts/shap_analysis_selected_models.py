"""
Title: shap_analysis_selected_models.py
Description: This script is used to perform SHAP analysis on the selected models.
Author: Jakob Jung
Date: 2023-11-21
"""

# Importing libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import shap
import pickle
from matplotlib.colors import LinearSegmentedColormap

#import data
pna_data = pd.read_csv("../data/pnas_predictors_mic_upec.tsv", delimiter="\t", header=0)
pna_data.head()


# step 2: split the data into a training and test set
# split the data into a training and test set
X = pna_data.drop(["pna_name", "gene_name", "pna_sequence", "target_seq", "upec_locus_tag", "MIC_UPEC", "MIC_K12",
                   "locus_tag"],
                  axis=1)
selected_features_mason = ['upec_tir_off_targets_2mm',
                     'purine_percentage', 'Tm',
                     'upec_tir_off_targets_1mm',
                     'PNA_molecular_weight',
                     'sc_bases',
                     'C9',
                     'CAI',
                     'MFE_UPEC',
                     'upec_total_off_targets_0mm',
                     'upec_total_off_targets_1mm',
                     'upec_total_off_targets_2mm',
                     'upec_tir_off_targets_0mm']

selected_features_normal = ['upec_tir_off_targets_2mm',
                     'purine_percentage', 'Tm',
                     'upec_tir_off_targets_1mm',
                     'sc_bases',
                     'gene_vulnerability_crispri',
                #     'PNA_molecular_weight',
                     'expression_upec',
                     'CAI',
                     'MFE_UPEC',
                     'upec_total_off_targets_0mm',
                     'upec_total_off_targets_1mm',
                     'upec_total_off_targets_2mm',
                     'upec_tir_off_targets_0mm']

X_normal = X[selected_features_normal]
X_mason = X[selected_features_mason]
y = np.log2(pna_data["MIC_UPEC"])

# Importing pickle file mason model
with open('../data/models_trained/rf_optimized_model_mason.sav', 'rb') as f:
    rf_mason = pickle.load(f)

# Importing pickle file mason model
with open('../data/models_trained/rf_optimized_model.sav', 'rb') as f:
    rf_normal = pickle.load(f)


def plot_shap_values(model, X, y, file_name):
    """
    Function to plot the shap values of a model
    :param model: model to use
    :param X: X data
    :param y: y data
    :param file_name: name of the file to save the plot as
    :return: plot of the shap values
    """
    # use shap to explain the predictions and the feature importance of the model
    model_shap = shap.TreeExplainer(model).shap_values(X)
    # plot the shap values
    f = plt.figure(figsize=(4, 6))
    # plot the shap values
    shap.summary_plot(model_shap, X, plot_type="bar", show=False, max_display=9, color="#4682B4")  # , max_display=10,

    # save the plot as svg
    plt.savefig("../analysis/" + file_name + "_summary_bar.svg", dpi=1200)

    # plot shap local values
    # Define the colors for the colormap
    colors = ['#4682B4', 'lightgray', '#B11226']  # Orange, Gray, Steel Blue
    # Create a colormap using LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

    f = plt.figure(figsize=(4, 6))
    shap.summary_plot(
        -model_shap, X, plot_type="dot", auto_size_plot=False, show=False, max_display=9  # , max_display=10
    )
    # Change the colormap of the artists
    for fc in plt.gcf().get_children():
        for fcc in fc.get_children():
            if hasattr(fcc, "set_cmap"):
                fcc.set_cmap(cmap)
    # save the plot as svg
    plt.savefig("../analysis/shap_plots/" + file_name + "_summary_dot.svg", dpi=1200)

    return plt


# plot the shap values
plot_shap_values(rf_normal, X_normal, y, "shap_values_rf_optimized_model")
plot_shap_values(rf_mason, X_mason, y, "shap_values_rf_optimized_model_mason")


