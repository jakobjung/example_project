# script that uses the trained model in pickle format (./data/models/rf_optimized_model.sav)
# to predict the PNAS class of all UPEC genomes in the dataset

import pandas as pd
import numpy as np
import pickle

from matplotlib.colors import LinearSegmentedColormap
from sklearn.ensemble import RandomForestRegressor
import shap
import matplotlib.pyplot as plt


# Load the trained model
model_rf = pickle.load(open('../data/models_trained/rf_optimized_model.sav', 'rb'))
model_rf_mason = pickle.load(open('../data/models_trained/rf_optimized_model_mason.sav', 'rb'))

# load the data (tsv) to predict
data_upec = pd.read_csv('../data/all_genes_all_pnas/all_pnas_all_egenes_upec.tsv', sep='\t')
data_upec.columns

# get feature order from the model
selected_features = ['upec_tir_off_targets_2mm',
                     'purine_percentage',
                     'Tm',
                     'upec_tir_off_targets_1mm',
                     'sc_bases',
                     #'gene_vulnerability_crispri',
                     'string_interactions',
                     #'crispri_log2FC_rousset',
                     'cas13_log2FC',
                     'PNA_molecular_weight',
                     'expression_upec',
                     'MFE_UPEC',
                     'upec_total_off_targets_0mm',
                     'upec_total_off_targets_1mm',
                     'upec_total_off_targets_2mm',
                     'upec_tir_off_targets_0mm']

selected_features_mason = ['Tm', 'sc_bases', 'PNA_molecular_weight', 'CAI', 'upec_total_off_targets_1mm', 'MFE_UPEC',
                           'purine_percentage']


upec_data_for_pred = data_upec[selected_features]
upec_data_for_pred_mason = data_upec[selected_features_mason]

# make sc bases relative to the length of the pna
upec_data_for_pred_mason["sc_bases"] = upec_data_for_pred_mason["sc_bases"] / 9

# predict the PNA MIC
y = model_rf.predict(upec_data_for_pred)
y_mason = model_rf_mason.predict(upec_data_for_pred_mason)

data_upec["predicted_MIC"] = y
data_upec["predicted_MIC_mason"] = y_mason

# save the data
data_upec.to_csv('../data/all_genes_all_pnas/all_pnas_all_egenes_upec_predicted.tsv', sep='\t', index=False)


# make a shap plot for the predictions of PNAs targeting the glmS gene
explainer = shap.TreeExplainer(model_rf)
shap_values = explainer.shap_values(upec_data_for_pred)
explanation = explainer(upec_data_for_pred)


shap.summary_plot(shap_values, upec_data_for_pred, plot_type="bar", show=False)
shap.summary_plot(shap_values, upec_data_for_pred, plot_type="dot", show=False)

shap.force_plot(explainer.expected_value, shap_values[0,:], upec_data_for_pred.iloc[0,:])

# find glmS gene from data_upec
glmS = data_upec[data_upec["K12_genename"] == "glmS"]
glmS

shap.save_html("../analysis/allgenes_allpnas/shapforce_glmS_best.html", shap.force_plot(explainer.expected_value,
                                                                              shap_values[2089,:],
                                                                              upec_data_for_pred.iloc[2089,:],
                                                                              show=False, plot_cmap='DrDb',
                                                                              text_rotation=0))


shap.save_html("../analysis/allgenes_allpnas/shapforce_glms_worst.html", shap.force_plot(explainer.expected_value,
                                                                              shap_values[2085,:],
                                                                              upec_data_for_pred.iloc[2085,:],
                                                                              show=False, plot_cmap='DrDb',
                                                                              text_rotation=0))

shap.save_html("../analysis/allgenes_allpnas/shapforce_glms_mediocre.html", shap.force_plot(explainer.expected_value,
                                                                              shap_values[2091,:],
                                                                              upec_data_for_pred.iloc[2091,:],
                                                                              show=False, plot_cmap='DrDb',

                                                                                            # Define the colors for the colormap
                                                                                            colors=['#4682B4',
                                                                                                    'lightgray',
                                                                                                    '#B11226'])) # Orange, Gray, Steel Blue

# Create a colormap using LinearSegmentedColormap
colors = ['#4682B4', 'lightgray', '#B11226']  # Orange, Gray, Steel Blue
cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

plt.figure(figsize=(6, 6))
shap.dependence_plot("Tm", shap_values, upec_data_for_pred, interaction_index="sc_bases",
                     cmap=cmap)

# save the plot as svg
plt.savefig("../analysis/allgenes_allpnas/b.svg", dpi=1200)

# close the plot
plt.close()


plt.figure(figsize=(6, 6))
shap.dependence_plot("Tm", -shap_values_mason, upec_data_for_pred_mason, interaction_index="sc_bases",
                     cmap=cmap)

# save the plot as svg
plt.savefig("../analysis/allgenes_allpnas/b.svg", dpi=1200)

# close the plot
plt.close()





# DO SAME FOR MASON MODEL
explainer_mason = shap.TreeExplainer(model_rf_mason)
shap_values_mason = explainer_mason.shap_values(upec_data_for_pred_mason)
explanation_mason = explainer_mason(upec_data_for_pred_mason)

shap.summary_plot(shap_values_mason, upec_data_for_pred_mason, plot_type="bar", show=False)
shap.summary_plot(shap_values_mason, upec_data_for_pred_mason, plot_type="dot", show=False)

shap.force_plot(explainer_mason.expected_value, shap_values_mason[0,:], upec_data_for_pred_mason.iloc[0,:])

# find glmS gene from data_upec
glmS = data_upec[data_upec["K12_genename"] == "glmS"]
glmS

shap.save_html("../analysis/allgenes_allpnas/shapforce_glmS_best_mason.html", shap.force_plot(explainer_mason.expected_value,
                                                                                shap_values_mason[2089,:],
                                                                                upec_data_for_pred_mason.iloc[2089,:],
                                                                                show=False, plot_cmap='DrDb',
                                                                                text_rotation=0))

shap.save_html("../analysis/allgenes_allpnas/shapforce_glms_worst_mason.html", shap.force_plot(explainer_mason.expected_value,
                                                                                shap_values_mason[2086,:],
                                                                                upec_data_for_pred_mason.iloc[2086,:],
                                                                                show=False, plot_cmap='DrDb',
                                                                                text_rotation=0))

shap.save_html("../analysis/allgenes_allpnas/shapforce_glms_mediocre_mason.html", shap.force_plot(explainer_mason.expected_value,
                                                                                shap_values_mason[2091,:],
                                                                                upec_data_for_pred_mason.iloc[2091,:],
                                                                                show=False, plot_cmap='DrDb',
                                                                                text_rotation=0))





