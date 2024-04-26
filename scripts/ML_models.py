"""
Title:      Test machine learning models on the PNA dataset
Author:     Jakob Jung
Date:       22-03-2023
Details:    I read in the data with the PNAs and all curated features. I then split the data into a training and test
            set. I then train a machine learning model on the training set and test it on the test set.
            I then save the model and the test set predictions. I use scikit-learn to train the models (random forests,
            SVMs, decision trees, etc.). I will make use of the shap packages to explain the predictions
            and the feature importance of the tree-based models.
"""

# Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import shap
import pickle
import os
import sys
import time
import datetime
import argparse
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
import xgboost
import matplotlib.pyplot as pl
from lightgbm import LGBMClassifier
from tpot import TPOTClassifier
#import autosklearn.classification



# Import data from the data folder
pna_data = pd.read_csv("../data/pnas_predictors_mic_upec.tsv",  delimiter="\t", header=0)

# view the data frame
pna_data.head()

# I will use cross-validation to train the models. I will use 10-fold cross-validation and 3 repeats.
cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=10, random_state=2)

# I will use the following models: random forest, SVM, decision tree, and logistic regression.
# I will also use auto-sklearn to do automated machine learning. I will use the shap package to explain the predictions
# and the feature importance of the tree-based models. I will also do automated machine learning with auto-sklearn.

# I will start with defining the features and the target variable.
# I will use all features except the PNA sequence and the target variable.
X = pna_data.drop(["pna_name", "gene_name", "pna_sequence", "target_seq", "e_coli_k12_inhibition", "upec_inhibition",
                   "inhibits_either", "upec_locus_tag", "MIC_UPEC", "MIC_K12", "ess_target_upec" ,"inhibition_of",
                   "OT_TIR_2mm_OT_2mm"], axis=1)
y = pna_data["MIC_UPEC"] < 6


r = 2
# split the data into a training and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=r, stratify=y)

# Use tpot to do automated machine learning
tpot = TPOTClassifier(generations=5, population_size=20, random_state=r, verbosity=2, n_jobs=-1)
tpot.fit(X_train, y_train)
print(tpot.score(X_test, y_test))
tpot.export('./tpot_pna_pipeline.py')


# train/test a random forest model using stratified 10-fold cross-validation using RepeatedStratifiedKFold
model = RandomForestClassifier()
scores_rf = cross_val_score(model, X, y, scoring='accuracy', cv=cv, n_jobs=-1)
f1_rf = cross_val_score(model, X, y, scoring='f1', cv=cv, n_jobs=-1)
auc_rf = cross_val_score(model, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)
precision_rf = cross_val_score(model, X, y, scoring='precision', cv=cv, n_jobs=-1)
recall_rf = cross_val_score(model, X, y, scoring='recall', cv=cv, n_jobs=-1)


print('Accuracy: %.3f (%.3f)' % (np.mean(scores_rf), np.std(scores_rf)))
print('F1: %.3f (%.3f)' % (np.mean(f1_rf), np.std(f1_rf)))
print('AUC: %.3f (%.3f)' % (np.mean(auc_rf), np.std(auc_rf)))
print('Precision: %.3f (%.3f)' % (np.mean(precision_rf), np.std(precision_rf)))
print('Recall: %.3f (%.3f)' % (np.mean(recall_rf), np.std(recall_rf)))


# fit the model on the training set
model = DecisionTreeClassifier()
model.fit(X, y)

# get decision tree model (accuracy, precision, recall, f1, auc)):
dt = DecisionTreeClassifier()
scores_dt = cross_val_score(dt, X, y, scoring='accuracy', cv=cv, n_jobs=-1)
f1_dt = cross_val_score(dt, X, y, scoring='f1', cv=cv, n_jobs=-1)
auc_dt = cross_val_score(dt, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)
precision_dt = cross_val_score(dt, X, y, scoring='precision', cv=cv, n_jobs=-1)
recall_dt = cross_val_score(dt, X, y, scoring='recall', cv=cv, n_jobs=-1)


print('Accuracy: %.3f (%.3f)' % (np.mean(scores_dt), np.std(scores_dt)))
print('F1: %.3f (%.3f)' % (np.mean(f1_dt), np.std(f1_dt)))
print('AUC: %.3f (%.3f)' % (np.mean(auc_dt), np.std(auc_dt)))
print('Precision: %.3f (%.3f)' % (np.mean(precision_dt), np.std(precision_dt)))
print('Recall: %.3f (%.3f)' % (np.mean(recall_dt), np.std(recall_dt)))


# do same for SVM
svm = SVC()

# train/test a random forest model using stratified 10-fold cross-validation using RepeatedStratifiedKFold
scores_svm = cross_val_score(svm, X, y, scoring='accuracy', cv=cv, n_jobs=-1)
f1_svm = cross_val_score(svm, X, y, scoring='f1', cv=cv, n_jobs=-1)
auc_svm = cross_val_score(svm, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)
precision_svm = cross_val_score(svm, X, y, scoring='precision', cv=cv, n_jobs=-1)
recall_svm = cross_val_score(svm, X, y, scoring='recall', cv=cv, n_jobs=-1)

print('Accuracy: %.3f (%.3f)' % (np.mean(scores_svm), np.std(scores_svm)))
print('F1: %.3f (%.3f)' % (np.mean(f1_svm), np.std(f1_svm)))
print('AUC: %.3f (%.3f)' % (np.mean(auc_svm), np.std(auc_svm)))
print('Precision: %.3f (%.3f)' % (np.mean(precision_svm), np.std(precision_svm)))
print('Recall: %.3f (%.3f)' % (np.mean(recall_svm), np.std(recall_svm)))

# fit the model on the training set
svm.fit(X, y)


# create model with XGBoost
xgb = xgboost.XGBClassifier()

# train/test a random forest model using stratified 10-fold cross-validation using RepeatedStratifiedKFold
scores_xgb = cross_val_score(xgb, X, y, scoring='accuracy', cv=cv, n_jobs=-1)
f1_xgb = cross_val_score(xgb, X, y, scoring='f1', cv=cv, n_jobs=-1)
auc_xgb = cross_val_score(xgb, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)
precision_xgb = cross_val_score(xgb, X, y, scoring='precision', cv=cv, n_jobs=-1)
recall_xgb = cross_val_score(xgb, X, y, scoring='recall', cv=cv, n_jobs=-1)

print('Accuracy: %.3f (%.3f)' % (np.mean(scores_xgb), np.std(scores_xgb)))
print('F1: %.3f (%.3f)' % (np.mean(f1_xgb), np.std(f1_xgb)))
print('AUC: %.3f (%.3f)' % (np.mean(auc_xgb), np.std(auc_xgb)))
print('Precision: %.3f (%.3f)' % (np.mean(precision_xgb), np.std(precision_xgb)))
print('Recall: %.3f (%.3f)' % (np.mean(recall_xgb), np.std(recall_xgb)))


xgb.fit(X, y)

# use shap to explain the predictions and the feature importance of the model
xgb_shap = shap.TreeExplainer(xgb).shap_values(X)

# plot the shap values
f = pl.figure(figsize=(4,6))
shap.summary_plot(
    xgb_shap, X,  plot_type="bar", show=False, max_display=10
)
pl.xlabel("mean(|SHAP value|)")
pl.savefig("../analysis/XGB_summary_bar_upec_5.pdf", dpi=1000)

# plot shap local values
f = pl.figure(figsize=(4,6))
shap.summary_plot(
    xgb_shap, X,  plot_type="dot", auto_size_plot=False, show=False, max_display=10
)
pl.xlabel("SHAP value (impact on model output)")
pl.savefig("../analysis/XGB_summary_dot_upec_5.pdf", dpi=1000)


# create lightgbm model
lgbm = LGBMClassifier()

# train/test a random forest model using stratified 10-fold cross-validation using RepeatedStratifiedKFold
scores_lgbm = cross_val_score(lgbm, X, y, scoring='accuracy', cv=cv, n_jobs=-1)
f1_lgbm = cross_val_score(lgbm, X, y, scoring='f1', cv=cv, n_jobs=-1)
auc_lgbm = cross_val_score(lgbm, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)
precision_lgbm = cross_val_score(lgbm, X, y, scoring='precision', cv=cv, n_jobs=-1)
recall_lgbm = cross_val_score(lgbm, X, y, scoring='recall', cv=cv, n_jobs=-1)


print('Accuracy: %.3f (%.3f)' % (np.mean(scores_lgbm), np.std(scores_lgbm)))
print('F1: %.3f (%.3f)' % (np.mean(f1_lgbm), np.std(f1_lgbm)))
print('AUC: %.3f (%.3f)' % (np.mean(auc_lgbm), np.std(auc_lgbm)))
print('Precision: %.3f (%.3f)' % (np.mean(precision_lgbm), np.std(precision_lgbm)))
print('Recall: %.3f (%.3f)' % (np.mean(recall_lgbm), np.std(recall_lgbm)))


# fit the model on the training set
lgbm.fit(X, y)

# use shap to explain the predictions and the feature importance of the model
lgbm_shap = shap.TreeExplainer(lgbm).shap_values(X)

# plot the shap values
f = pl.figure(figsize=(4,6))
shap.summary_plot(
    lgbm_shap, X,  plot_type="bar", show=False, max_display=8
)
pl.xlabel("mean(|SHAP value|)")
pl.savefig("../analysis/LGBM_summary_bar_upec.pdf", dpi=1000)


# same for logistic regression
lr = LogisticRegression()

# train/test a random forest model using stratified 10-fold cross-validation using RepeatedStratifiedKFold
scores_lr = cross_val_score(lr, X, y, scoring='accuracy', cv=cv, n_jobs=-1)
f1_lr = cross_val_score(lr, X, y, scoring='f1', cv=cv, n_jobs=-1)
auc_lr = cross_val_score(lr, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)
precision_lr = cross_val_score(lr, X, y, scoring='precision', cv=cv, n_jobs=-1)
recall_lr = cross_val_score(lr, X, y, scoring='recall', cv=cv, n_jobs=-1)

# fit the model on the training set
lr.fit(X, y)


print('Accuracy: %.3f (%.3f)' % (np.mean(scores_lr), np.std(scores_lr)))
print('F1: %.3f (%.3f)' % (np.mean(f1_lr), np.std(f1_lr)))
print('AUC: %.3f (%.3f)' % (np.mean(auc_lr), np.std(auc_lr)))
print('Precision: %.3f (%.3f)' % (np.mean(precision_lr), np.std(precision_lr)))
print('Recall: %.3f (%.3f)' % (np.mean(recall_lr), np.std(recall_lr)))




# create boxplot of the model performance. make boxes steel-blue and median black. order by median accuracy
data = [scores_xgb, scores_rf,  scores_lr, scores_svm, scores_dt]
fig = pl.figure(figsize=(5,5))
ax = fig.add_subplot(111)
bp = ax.boxplot(data, patch_artist=True, showfliers=False, showmeans=False, meanline=False,
                medianprops=dict(linestyle='-', linewidth=1, color='black'))
ax.set_xticklabels([ 'XGB', 'RF', 'LR', 'SVM', 'DT'])
ax.set_ylabel('Accuracy')

pl.savefig("../analysis/accuracy_boxplot5.pdf", dpi=400)
pl.show()

# create boxplot of the model performance in f1. make boxes steel-blue and median black. order by median accuracy
data = [f1_xgb, f1_rf,  f1_lr, f1_svm, f1_dt]
fig = pl.figure(figsize=(5,5))
ax = fig.add_subplot(111)
bp = ax.boxplot(data, patch_artist=True, showfliers=False, showmeans=False, meanline=False,
                medianprops=dict(linestyle='-', linewidth=1, color='black'))
ax.set_xticklabels(['XGB', 'RF', 'LR', 'SVM', 'DT'])
ax.set_ylabel('F1')
pl.savefig("../analysis/boxplot_f1_5.pdf", dpi=400)
pl.show()

# do same for AUC
data = [auc_xgb, auc_rf,  auc_lr, auc_svm, auc_dt]
fig = pl.figure(figsize=(5,5))
ax = fig.add_subplot(111)
bp = ax.boxplot(data, patch_artist=True, showfliers=False, showmeans=False, meanline=False,
                medianprops=dict(linestyle='-', linewidth=1, color='black'))
ax.set_xticklabels(['XGB', 'RF', 'LR', 'SVM', 'DT'])
ax.set_ylabel('AUC')
#ax.set_ylim([0.6, 1])
pl.savefig("../analysis/boxplot_auc5.pdf", dpi=400)
pl.show()

# same for precision
data = [precision_xgb, precision_rf,  precision_lr, precision_svm, precision_dt]
fig = pl.figure(figsize=(6,5))
ax = fig.add_subplot(111)
bp = ax.boxplot(data, patch_artist=True, showfliers=False, showmeans=False, meanline=False,
                medianprops=dict(linestyle='-', linewidth=1, color='black'))
ax.set_xticklabels(['XGB', 'RF', 'LR', 'SVM', 'DT'])
ax.set_ylabel('Precision')
pl.savefig("../analysis/boxplot_precision5.pdf", dpi=400)
pl.show()

# same for recall
data = [recall_xgb, recall_rf,  recall_lr, recall_svm, recall_dt]
fig = pl.figure(figsize=(6,5))
ax = fig.add_subplot(111)
bp = ax.boxplot(data, patch_artist=True, showfliers=False, showmeans=False, meanline=False,
                medianprops=dict(linestyle='-', linewidth=1, color='black'))
ax.set_xticklabels(['XGB', 'RF', 'LR', 'SVM', 'DT'])
ax.set_ylabel('Recall')
pl.savefig("../analysis/boxplot_recall5.pdf", dpi=400)
pl.show()


# create a roc curve for each model

# get package for roc curve
from sklearn.metrics import roc_curve

fpr_xgb, tpr_xgb, thresholds_xgb = roc_curve(y, xgb.predict_proba(X)[:,1])
fpr_rf, tpr_rf, thresholds_rf = roc_curve(y, model.predict_proba(X)[:,1])
fpr_lr, tpr_lr, thresholds_lr = roc_curve(y, lr.predict_proba(X)[:,1])
fpr_svm, tpr_svm, thresholds_svm = roc_curve(y, svm.predict_proba(X)[:,1])
fpr_dt, tpr_dt, thresholds_dt = roc_curve(y, dt.predict_proba(X)[:,1])


# evaluate all the models on the dataset X2 and y2, and generate the same plots as above
# get the data
