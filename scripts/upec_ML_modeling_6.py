"""
Title:      Test machine learning models on the PNA dataset
Author:     Jakob Jung
Date:       22-03-2023
Details:    I read in the data with the PNAs and all curated features. I then split the data into a training and test
            set. I then train a machine learning model on the training set and test it on the test set. The data
            consists of MIC as output and >50 predictors. I then save the model and the test set predictions.
            I use scikit-learn to train the models (random forest ,SVMs, decision trees,XGBoost etc.). To improve it,
            I apply automated machine learning for feature selection and hyperparameter tuning. I will make use
            of the shap packages to explain the predictions and the feature importance of the tree-based models.
            The problem is a regression problem, so I will use regression models.
"""
import joblib
# Import packages for regression models
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# import ListedColormap
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import seaborn as sns

import shap
import xgboost as xgb
# import pickle:
import pickle

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.metrics import make_scorer
# import GradientBoostingRegressor
from sklearn.ensemble import GradientBoostingRegressor
# import lightgbm regressor
import lightgbm as lgb
# Import packages for automated machine learning
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor
# import lr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from BorutaShap import BorutaShap
from sklearn.preprocessing import MaxAbsScaler
from sklearn.linear_model import LassoLarsCV
# get svm regressor
from sklearn.svm import SVR
# get hyperopt
from hyperopt import hp, fmin, tpe, STATUS_OK, Trials, space_eval

# step 1: read in the data
# read in the data
pna_data = pd.read_csv("../data/pnas_predictors_mic_upec.tsv", delimiter="\t", header=0)
pna_data.head()


# step 2: split the data into a training and test set
# split the data into a training and test set
X = pna_data.drop(["pna_name", "gene_name", "pna_sequence", "target_seq", "upec_locus_tag", "MIC_UPEC", "MIC_K12",
                   "locus_tag"],
                  axis=1)
y = np.log2(pna_data["MIC_UPEC"])


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
print("The training set has {} samples and the test set has {} samples.".format(X_train.shape[0], X_test.shape[0]))


# make a function for above code, plotting the distribution of the MIC values in the training and test set next to one
# another. MIC values need to be made a factor variable first
def plot_mic_distribution(y_train, y_test, plot_name, file_name):
    y_train_fact = y_train.astype('category')
    y_test_fact = y_test.astype('category')
    # create a df to be able to use seaborn it is easier
    y_train_df = pd.DataFrame(y_train_fact)
    y_test_df = pd.DataFrame(y_test_fact)
    # rename the column
    y_train_df.columns = ["MIC"]
    y_test_df.columns = ["MIC"]
    # add a column with the dataset name
    y_train_df["dataset"] = "training set"
    y_test_df["dataset"] = "test set"
    # combine the two dfs
    y_df = pd.concat([y_train_df, y_test_df])

    # use the magma colormap from matplotlib, and use 2 colors from it
    cmap = matplotlib.cm.get_cmap('inferno')
    # make the bar plot
    fig, ax = plt.subplots(figsize=(5, 5))
    # make the bar plot
    sns.countplot(x="MIC", hue="dataset", data=y_df, palette=[cmap(0.15), cmap(0.80)])
    # rename the xticks
    ax.set_xticklabels(["1.25", "2.5", "5", "10", ">10"])
    # set the x and y labels
    ax.set_xlabel("MIC values")
    ax.set_ylabel("Number of samples")
    # remove legend title
    ax.legend(title="")
    # set the title
    ax.set_title(plot_name)
    # save the plot as svg
    plt.savefig("../analysis/" + file_name + ".svg")


# make a bar plot showing the distribution of the MIC values in the training and test set next to one another. MIC
# values need to be made a factor variable first
plot_mic_distribution(y_train, y_test, "Distribution of the MICs", "mic_distribution_training_test_set")


# step 3: train different  machine learning models on the training set and test it on the test set.
# I will use scikit-learn to train the models (random forest ,SVMs, decision trees,XGBoost, linear regression, etc.).
# I will make use of the shap packages to explain the predictions and the feature importance of the tree-based models.
# The problem is a regression problem (predicting MIC), so I will use regression models.

# step 3.1: train a random forest model
# Make a custom metric function which calculates the % of predictions lying within +-1 dilution of the true MIC!
# this metric will later be used to evaluate the models


def within_1_dil(y_true, y_actual):
    # calculate the % of predictions lying within +-1 dilution of the true MIC
    # calculate the difference between the true and predicted MICs
    diff = abs(y_true - y_actual)
    # calculate the % of predictions lying within +-1 dilution of the true MIC
    within_1_dil = len(diff[diff <= 1]) / len(diff)
    return within_1_dil


# make custom scorer for tpot:
within_1_dil_scorer = make_scorer(within_1_dil, greater_is_better=True)

# run the function on the X and y with the r2 score
# make r2 to score function
r2_scorer = make_scorer(r2_score)

# same for RMSE
# make RMSE scorer
rmse_scorer = make_scorer(mean_squared_error, squared=False)


# make a score function that gives true if predicted MIC is below 10 and actual MIC is below 10, or if predicted MIC is
# above 10 and actual MIC is above 10. This is to calculate the accuracy of the model in predicting whether the MIC is
# below or above 10
def within_10_precision(y_true, y_pred):
    """
    Function to calculate the accuracy of the model in predicting whether the MIC is below or above 10
    :param y_true: true y values
    :param y_pred: predicted y values
    :return: accuracy score
    """
    # make y_pred and y_true arrays snd exp2 them
    y_pred = np.exp2(y_pred)
    y_true = np.exp2(y_true)
    # get true positives
    tp = sum((y_pred <= 10) & (y_true <= 10))
    # get false positives
    fp = sum((y_pred <= 10) & (y_true > 10))
    # calculate the accuracy score
    accuracy = tp / (tp + fp)
    return accuracy


# do same for within 10
within_10_precision_score = make_scorer(within_10_precision)


def within_5_precision(y_true, y_pred):
    """
    Function to calculate the accuracy of the model in predicting whether the MIC is below or above 10
    :param y_true: true y values
    :param y_pred: predicted y values
    :return: accuracy score
    """
    # make y_pred and y_true arrays snd exp2 them
    y_pred = np.exp2(y_pred)
    y_true = np.exp2(y_true)
    # get true positives
    tp = sum((y_pred <= 5) & (y_true <= 5))
    # get false positives
    fp = sum((y_pred <= 5) & (y_true > 5))
    # calculate the accuracy score
    accuracy = tp / (tp + fp)
    return accuracy


# do same for within 5
within_5_precision_score = make_scorer(within_5_precision)


# do same with recall
def within_10_recall(y_true, y_pred):
    """
    Function to calculate the accuracy of the model in predicting whether the MIC is below or above 10
    :param y_true: true y values
    :param y_pred: predicted y values
    :return: accuracy score
    """
    # make y_pred and y_true arrays snd exp2 them
    y_pred = np.exp2(y_pred)
    y_true = np.exp2(y_true)
    # get true positives
    tp = sum((y_pred <= 10) & (y_true <= 10))
    # get false negatives
    fn = sum((y_pred > 10) & (y_true <= 10))
    # calculate the accuracy score
    rec = tp / (tp + fn)
    return rec


# do same for within 10
within_10_recall_score = make_scorer(within_10_recall)


# do same for f1 score
def within_10_f1(y_true, y_pred):
    """
    Function to calculate the accuracy of the model in predicting whether the MIC is below or above 10
    :param y_true: true y values
    :param y_pred: predicted y values
    :return: accuracy score
    """
    # make y_pred and y_true arrays snd exp2 them
    y_pred = np.exp2(y_pred)
    y_true = np.exp2(y_true)
    # get true positives
    tp = sum((y_pred <= 10) & (y_true <= 10))
    # get false positives
    fp = sum((y_pred <= 10) & (y_true > 10))
    # get false negatives
    fn = sum((y_pred > 10) & (y_true <= 10))
    # calculate the accuracy score
    prec = tp / (tp + fp)
    rec = tp / (tp + fn)
    f1 = 2 * ((prec * rec) / (prec + rec))
    return f1


# do same for within 10
within_10_f1_score = make_scorer(within_10_f1)


# do samefor recall within 5
def within_5_recall(y_true, y_pred):
    """
    Function to calculate the accuracy of the model in predicting whether the MIC is below or above 10
    :param y_true: true y values
    :param y_pred: predicted y values
    :return: accuracy score
    """
    # make y_pred and y_true arrays snd exp2 them
    y_pred = np.exp2(y_pred)
    y_true = np.exp2(y_true)
    # get true positives
    tp = sum((y_pred <= 5) & (y_true <= 5))
    # get false negatives
    fn = sum((y_pred > 5) & (y_true <= 5))
    # calculate the accuracy score
    rec = tp / (tp + fn)
    return rec


# do same for within 5
within_5_recall_score = make_scorer(within_5_recall)


# do within 5 f1 score
def within_5_f1(y_true, y_pred):
    """
    Function to calculate the accuracy of the model in predicting whether the MIC is below or above 10
    :param y_true: true y values
    :param y_pred: predicted y values
    :return: accuracy score
    """
    # make y_pred and y_true arrays snd exp2 them
    y_pred = np.exp2(y_pred)
    y_true = np.exp2(y_true)
    # get true positives
    tp = sum((y_pred <= 5) & (y_true <= 5))
    # get false positives
    fp = sum((y_pred <= 5) & (y_true > 5))
    # get false negatives
    fn = sum((y_pred > 5) & (y_true <= 5))
    # calculate the accuracy score
    prec = tp / (tp + fp)
    rec = tp / (tp + fn)
    f1 = 2 * ((prec * rec) / (prec + rec))
    return f1


# do same for within 5
within_5_f1_score = make_scorer(within_5_f1)


# Now create multiple models and evaluate them using cross validation. Use Linear regression, Lassolars, SVR, random forest,
# XGBoost, GradientBoostingRegressor, LightGBM
# define the models
models = []
models.append(('LR', LinearRegression()))
models.append(('LASSO', LassoLarsCV(normalize=True)))
models.append(('RF', RandomForestRegressor()))
models.append(('XGB', xgb.XGBRegressor()))
models.append(('GBoost', GradientBoostingRegressor()))
models.append(('LightGBM', lgb.LGBMRegressor()))




names_model = ["Linear regression", "LassoLarsCV", "Random forest", "XGBoost", "GBoost", "LightGBM",
               "SVR", "GPR"]


# create function running multiple models cv and returning the results
def run_multiple_cv_runs(models, model_names, X, y, cv_splits, num_runs, cv_metric=within_1_dil_scorer):
    results = []
    for m in range(len(models)):
        model_name = model_names[m]
        model = models[m]
        scores_list = []
        print(model_name)
        for i in range(num_runs):
            kf = KFold(n_splits=cv_splits, shuffle=True, random_state=i)
            scores = cross_val_score(model, X, y, cv=kf, scoring=cv_metric)
            # append the scores to the scores list if not nan
            scores_list.append(scores[~np.isnan(scores)])

        # Calculate the median of the scores for each run
        median_scores = [median_score.mean() for median_score in scores_list]

        # Store the median scores in a DataFrame
        results.append({'Model': model_name, 'Median_Scores': median_scores})

    df_results = pd.DataFrame(results)
    # expand df with the median scores so that median scores are expanded to 10 rows
    df_results = df_results.explode('Median_Scores')

    return df_results


# list all scoring metrics
scoring_metrics = [within_1_dil_scorer, r2_scorer, rmse_scorer, within_10_precision_score,
                   within_10_recall_score, within_10_f1_score, within_5_precision_score, within_5_recall_score,
                   within_5_f1_score]
scoring_metrics_names = ["within-1-tier accuracy", "R2 score", "Root mean squared error", "Precision (MIC <= 10)",
                         "Recall (MIC <= 10)", "F1 score (MIC <= 10)", "Precision (MIC <= 5)", "Recall (MIC <= 5)",
                         "F1 score (MIC <= 5)"]

# predict the MIC values using RF
# train rf model pn training set
rf = RandomForestRegressor()
rf.fit(X_train, y_train)
# make predictions on the test set
y_pred_rf = rf.predict(X_test)



# go through all metrics and run 10 fold cv on all models with that metric. then put dfs into one dataframe,
# # adding a column with the model name
# df_results_all_metrics = pd.DataFrame()
# for i in range(len(scoring_metrics)):
#     print(scoring_metrics_names[i])
#     df_results = run_multiple_cv_runs([model[1] for model in models],
#                                       names_model, X, y, 10, 10,
#                                       cv_metric=scoring_metrics[i])
#     df_results["scoring_metric"] = scoring_metrics_names[i]
#     # add to df_results_all_metrics
#     df_results_all_metrics = pd.concat([df_results_all_metrics, df_results])
#
# # save the df_results_all_metrics as csv, so I do not have to run the above code again
# df_results_all_metrics.to_csv("../data/df_results_all_metrics.csv", index=False)








# read in the df_results_all_metrics
df_results_all_metrics = pd.read_csv("../data/df_results_all_metrics.csv", header=0)

# now plot the results as boxplots and
# create a function of boxplot of the results, add points


def plot_multiple_cv_runs(df_results, plot_name, file_name, location="lower left", width=20, height=6):
    # change the style of the plot
    plt.figure(figsize=(width, height))
    # create a color palette. use the inferno palette, but reverse it
    cmap = plt.get_cmap('viridis', 512)
    c1 = matplotlib.colors.to_hex(cmap(1 / 6))
    c2 = matplotlib.colors.to_hex(cmap(2 / 6))
    c3 = matplotlib.colors.to_hex(cmap(3 / 6))
    c4 = matplotlib.colors.to_hex(cmap(4 / 6))
    c5 = matplotlib.colors.to_hex(cmap(5 / 6))
    c6 = matplotlib.colors.to_hex(cmap(6 / 6))
    # change style of the plot
    sns.set_style("whitegrid")
    # add points
    sns.swarmplot(x="scoring_metric", y="Median_Scores", data=df_results, hue="Model", palette=[c1, c2, c3, c4, c5, c6],
                  dodge=True, alpha=0.75, size=4, edgecolor="grey", linewidth=1, legend=False)
    # add the boxplot, change fill color and line color
    sns.boxplot(x="scoring_metric", y="Median_Scores", data=df_results, hue="Model",
                palette=[c1, c2, c3, c4, c5, c6], linewidth=1)

    plt.xticks(rotation=45, ha="right")
    # increase the font size of axis ticks
    plt.tick_params(axis='both', which='major', labelsize=14)
    # remove y axis label
    plt.xlabel("")
    plt.ylabel("")
    # increase size of y axis label
    plt.ylabel("")
    # add title
    plt.title(plot_name, fontsize=16, fontweight="bold", pad=20)
    # add legend bottom left
    plt.legend(loc=location, fontsize=14)
    # save the plot as svg
    plt.rcParams['savefig.facecolor'] = 'white'
    plt.savefig("../analysis/" + file_name + ".svg", bbox_inches='tight')

    return plt


# drop recall mic <= 5 from df_results_all_metrics
df_results_all_metrics_pl = df_results_all_metrics[df_results_all_metrics["scoring_metric"] != "Recall (MIC <= 5)"]
# also drop f1 score mic <= 5
df_results_all_metrics_pl = df_results_all_metrics_pl[df_results_all_metrics_pl["scoring_metric"] != "F1 score (MIC <= 5)"]
# also drop precision <= 5
df_results_all_metrics_pl = df_results_all_metrics_pl[df_results_all_metrics_pl["scoring_metric"] != "Precision (MIC <= 5)"]
# also drop r2 score
df_results_all_metrics_pl = df_results_all_metrics_pl[df_results_all_metrics_pl["scoring_metric"] != "R2 score"]


# plot the results
plt_1tier = plot_multiple_cv_runs(df_results_all_metrics_pl, "10-fold cross validation results raw models",
                                  "10foldcv_results_all_metrics_raw_models", width=15, height=6, location="upper left")
# plot all metrics
plt_all_metrics = plot_multiple_cv_runs(df_results_all_metrics, "10-fold cross validation results raw models",
                                        "10foldcv_results_all_metrics_raw_models_all_metrics")


# step 3.2: Check for feature importance in RF model with BorutaShap
# define the model
rf = RandomForestRegressor()

# # define BorutaShap feature selection method
# feat_selector = BorutaShap(importance_measure="shap", classification=False)
#
# # find all relevant features
# feat_selector.fit(X=X, y=y, n_trials=100, random_state=42)
#
#
# # get plot of the feature importance
# feat_selector.plot(which_features='all', figsize=(16, 5))
# # add ylim from 10^-0.8 to 10^1
# plt.ylim(10**-0.8, 10**1)
# # increase the font size of axis ticks
# plt.tick_params(axis='both', which='major', labelsize=12)
# # add lines for each feature (vertical lines)
# for i in range(len(X.columns)+4):
#    plt.axvline(i, color="grey", linestyle="--", linewidth=1, alpha=0.3)
# # save the plot as svg
# plt.savefig("../analysis/feature_importances_boruta_shap.svg", bbox_inches='tight')

# get the selected features
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

# now generate new X with only the selected features
X_selected = X[selected_features]





# # now run 10 fold cv on the RF model with only the selected features
# # go through all metrics and run 10 fold cv on all models with that metric. then put dfs into one dataframe,
# # adding a column with the model name
# df_results_all_metrics_selected_features = pd.DataFrame()
# for i in range(len(scoring_metrics)):
#     print(scoring_metrics_names[i])
#     df_results = run_multiple_cv_runs([model[1] for model in models],
#                                                        names_model, X_selected, y, 10, 10,
#                                                        cv_metric=scoring_metrics[i])
#     df_results["scoring_metric"] = scoring_metrics_names[i]
#     # add to df_results_all_metrics
#     df_results_all_metrics_selected_features = pd.concat([df_results_all_metrics_selected_features, df_results])
#
# # save the df_results_all_metrics as csv, so I do not have to run the above code again
# df_results_all_metrics_selected_features.to_csv("../data/df_results_all_metrics_selected_features.csv", index=False)







# read in the df_results_all_metrics
df_results_all_metrics_selected_features = pd.read_csv("../data/df_results_all_metrics_selected_features.csv", header=0)


# drop recall mic <= 5 from df_results_all_metrics
df_results_all_metrics_selected_features_pl = df_results_all_metrics_selected_features[
    df_results_all_metrics_selected_features["scoring_metric"] != "Recall (MIC <= 5)"]
# also drop r2 score
df_results_all_metrics_selected_features_pl = df_results_all_metrics_selected_features_pl[
    df_results_all_metrics_selected_features_pl["scoring_metric"] != "R2 score"]
# also drop f1 score mic <= 5
df_results_all_metrics_selected_features_pl = df_results_all_metrics_selected_features_pl[
    df_results_all_metrics_selected_features_pl["scoring_metric"] != "F1 score (MIC <= 5)"]
# also drop precision <= 5
df_results_all_metrics_selected_features_pl = df_results_all_metrics_selected_features_pl[
    df_results_all_metrics_selected_features_pl["scoring_metric"] != "Precision (MIC <= 5)"]


# plot the results
plot_multiple_cv_runs(df_results_all_metrics_selected_features_pl,
                      "10-fold cross validation results raw models with selected features",
                      "10foldcv_results_all_metrics_selected_features", width=15, height=6, location="upper left")

plot_multiple_cv_runs(df_results_all_metrics_selected_features,
                      "10-fold cross validation results raw models with selected features",
                      "10foldcv_results_all_metrics_selected_features_all_metrics")

# step 3.3: do hyperparameter tuning for all the models. Use grid search and 10 fold cv
# split data
X_train, X_test, y_train, y_test = train_test_split(X_selected, y, test_size=0.2, random_state=42)
# define. restrict the models to my models
models = []
models.append(('LR', LinearRegression()))
models.append(('LASSO', LassoLarsCV()))
models.append(('RF', RandomForestRegressor()))
models.append(('XGB', xgb.XGBRegressor()))
models.append(('GBoost', GradientBoostingRegressor()))
models.append(('LightGBM', lgb.LGBMRegressor()))

names_model = ["Linear regression", "LassoLarsCV", "Random forest", "XGBoost", "GradientBoostingRegressor", "LightGBM"]

# define the parameter grids for the models. rmrmber I have around 550 data points and 14 features
# linear regression
param_grid_lr = {
    'fit_intercept': [True, False],
    'copy_X': [True, False],
    'positive': [True, False],
    'n_jobs': [None, 1, 2, 4, 8]
}
# lasso lars cv
param_grid_lasso = {
    'eps': [1e-5, 1e-4, 1e-3, 1e-2, 2.220446049250313e-16],
    'max_iter': [500, 1000, 1500]
}

# random forest
param_grid_rf = {
    'n_estimators': [100, 200, 300, 400],
    'max_features': ['auto', 'sqrt', 1.0],
    'max_depth': [None, 10, 20, 30, 40],
    'min_samples_split': [2, 5],
    'min_samples_leaf': [1, 2],
    'bootstrap': [True, False]
}
# XGBoost
param_grid_xgb = {
    'n_estimators': [100, 200, 300],
    'max_depth': [10, 20, 30],
    'learning_rate': [0.01, 0.1, 0.2],
    'subsample': [0.7, 0.8, 0.9],

}
# GradientBoostingRegressor
param_grid_gboost = {
    'n_estimators': [100, 200, 300],
    'max_depth': [10, 20, 30],
    'learning_rate': [0.01, 0.1, 0.2],
    'subsample': [0.7, 0.8, 0.9]
}
# LightGBM
param_grid_lgb = {
    'n_estimators': [100, 200, 300],
    'max_depth': [10, 20, 30],
    'learning_rate': [0.01, 0.1, 0.2],
    'subsample': [0.7, 0.8, 0.9]
}


# define the grid search function
def grid_search_cv(model, param_grid, X_train, y_train, cv_splits, num_runs, cv_metric=within_1_dil_scorer):
    # define the grid search
    grid_search = GridSearchCV(estimator=model, param_grid=param_grid, cv=cv_splits, n_jobs=-1, verbose=1,
                               scoring=cv_metric)
    # run the grid search
    grid_search.fit(X_train, y_train)
    # get the best parameters
    best_params = grid_search.best_params_
    # get the best score
    best_score = grid_search.best_score_
    # get the best model
    best_model = grid_search.best_estimator_
    # get the cv results
    cv_results = grid_search.cv_results_
    return best_params, best_score, best_model, cv_results


# write above model comparison as function
def compare_optimized_non_optimized_model(model_1, model_2):
    """
    Function to compare the performance of an optimized and non-optimized model
    :param model_1: optimized model
    :param model_2: non-optimized
    :return: print statement with the results
    """
    # fit the models
    model_1.fit(X_train, y_train)
    model_2.fit(X_train, y_train)
    # make predictions on the test set
    y_pred_model_1 = model_1.predict(X_test)
    y_pred_model_2 = model_2.predict(X_test)
    # calculate the within 1 dilution accuracy
    within_1_dil_model_1 = within_1_dil(y_test, y_pred_model_1)
    within_1_dil_model_2 = within_1_dil(y_test, y_pred_model_2)
    # print the results
    print("The within 1 dilution accuracy of the optimized model is {} and the within 1 dilution accuracy of "
          "the non-optimized model is {}".format(within_1_dil_model_1, within_1_dil_model_2))
    # check for differences in hyperparameters of 2 models
    # get the hyperparameters of the optimized model
    hyperparameters_model_1 = model_1.get_params()
    # get the hyperparameters of the non-optimized model
    hyperparameters_model_2 = model_2.get_params()
    # get the hyperparameters that are different
    hyperparameters_different = {k: hyperparameters_model_1[k] for k in hyperparameters_model_1 if
                                 hyperparameters_model_1[k] != hyperparameters_model_2[k]}
    # print the hyperparameters that are different
    print("The hyperparameters that are different between the optimized and non-optimized model are: {}".format(
        hyperparameters_different))
    # print them for the optimized model
    print("The hyperparameters of the optimized model are: {}".format(hyperparameters_model_1))
    # print them for the non-optimized model
    print("The hyperparameters of the non-optimized model are: {}".format(hyperparameters_model_2))


# use the parameter grids and put them into a list
param_grids = [param_grid_lr, param_grid_lasso, param_grid_rf, param_grid_xgb, param_grid_gboost, param_grid_lgb]






# # run the grid search for all models
# best_params_list = []
# best_score_list = []
# best_model_list = []

# for i in range(len(models)):
#     print(names_model[i])
#     best_params, best_score, best_model, cv_results = grid_search_cv(models[i][1], param_grids[i],
#                                                                     X_train, y_train, 10, 10)
#     models = []
#     models.append(('LR', LinearRegression()))
#     models.append(('LASSO', LassoLarsCV()))
#     models.append(('RF', RandomForestRegressor()))
#     models.append(('XGB', xgb.XGBRegressor()))
#     models.append(('GBoost', GradientBoostingRegressor()))
#     models.append(('LightGBM', lgb.LGBMRegressor()))
#
#     # fit normal model
#     model_old = models[i][1]
#
#     # get the score of the normal model
#     compare_optimized_non_optimized_model(best_model, model_old)
#     # append to lists
#     best_params_list.append(best_params)
#     best_score_list.append(best_score)
#     best_model_list.append(best_model)
#
#






# write the best model list directly here:
# [LinearRegression(), LassoLarsCV(eps=1e-05, normalize=True), RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt',
#                       min_samples_leaf=2), RandomForestRegressor(random_state=42), XGBRegressor(base_score=None, booster=None, callbacks=None,
#              colsample_bylevel=None, colsample_bynode=None,
#              colsample_bytree=0.8, early_stopping_rounds=None,
#              enable_categorical=False, eval_metric=None, feature_types=None,
#              gamma=None, gpu_id=None, grow_policy=None, importance_type=None,
#              interaction_constraints=None, learning_rate=0.1, max_bin=None,
#              max_cat_threshold=None, max_cat_to_onehot=None,
#              max_delta_step=None, max_depth=30, max_leaves=None,
#              min_child_weight=None, missing=nan, monotone_constraints=None,
#              n_estimators=100, n_jobs=None, num_parallel_tree=None,
#              predictor=None, random_state=None, ...), GradientBoostingRegressor(max_depth=20, n_estimators=300, subsample=0.8), LGBMRegressor(learning_rate=0.01, max_depth=10, n_estimators=300, subsample=0.7)]

best_model_list = [LinearRegression(), LassoLarsCV(eps=1e-05, normalize=True),
                     RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt', min_samples_leaf=2),
                  # RandomForestRegressor(random_state=42),
                        xgb.XGBRegressor(max_depth=30, n_estimators=100, learning_rate=0.1, subsample=0.9),
                        GradientBoostingRegressor(max_depth=20, n_estimators=300, subsample=0.8),
                        lgb.LGBMRegressor(learning_rate=0.01, max_depth=10, n_estimators=300, subsample=0.7)]



best_score_list = [0.7633333333333334, 0.7611111111111111, 0.7988383838383839, 0.7764141414141414, 0.7921717171717171,
                   0.7853535353535354]

# get scores of raw models
raw_scores_list_test = []
for i in range(len(models)):
    model = models[i][1]
    # fit the model
    model.fit(X_train, y_train)
    # make predictions
    y_pred = model.predict(X_test)
    # calculate the within 1 dilution accuracy
    within_1_dil = within_1_dil_scorer(model, X_test, y_test)
    # append to list
    raw_scores_list_test.append(within_1_dil)

best_score_list_test = []
# get scores of optimized models
for i in range(len(best_model_list)):
    model = best_model_list[i]
    # fit the model
    model.fit(X_train, y_train)
    # make predictions
    y_pred = model.predict(X_test)
    # calculate the within 1 dilution accuracy
    within_1_dil = within_1_dil_scorer(model, X_test, y_test)
    # append to list
    best_score_list_test.append(within_1_dil)

# OK, now generate a plot with raw scores and optimized accuracy scores, as a barplot. X axis is the model name,
# y axis is the accuracy score. color by whether the model is optimized or not using inferno palette
# create a dataframe with the scores
df_scores = pd.DataFrame(
    {"model": names_model, "score": raw_scores_list_test, "optimized": ["raw"] * len(raw_scores_list_test)})
# add the optimized scores
df_scores_optimized = pd.DataFrame(
    {"model": names_model, "score": best_score_list_test, "optimized": ["h-tuning"] * len(best_score_list_test)})
# append to df_scores
df_scores = pd.concat([df_scores, df_scores_optimized])

# create a barplot
plt.figure(figsize=(6, 6))
# make a color palette
cmap = plt.get_cmap('inferno', 512)
c1 = matplotlib.colors.to_hex(cmap(1 / 6))
c2 = matplotlib.colors.to_hex(cmap(4 / 6))
# change style of the plot
sns.set_style("whitegrid")
# make the barplot
sns.barplot(x="model", y="score", data=df_scores, hue="optimized", palette=[c1, c2])
# add legend for the colors
plt.legend(loc='upper left', fontsize=14)
# add title
plt.title("Test accuracy of raw models vs. optimized models", fontsize=16, fontweight="bold", pad=20)
# make y axis title, within 1 dilution accuracy
plt.ylabel("Within 1 tier accuracy", fontsize=15)
# remove x axis label
plt.xlabel("")
# make ylimits to 0.6-0.8
plt.ylim(0.6, 0.8)
# rotate x axis labels
plt.xticks(rotation=45, ha="right")
# save the plot as svg
plt.rcParams['savefig.facecolor'] = 'white'
plt.savefig("../analysis/10foldcv_results_raw_models_vs_optimized_models.svg", bbox_inches='tight')








# # ok now run the 10 fold cv on the optimized models and get the results, as I have done before
# # go through all metrics and run 10 fold cv on all models with that metric. then put dfs into one dataframe,
# # adding a column with the model name
# df_results_all_metrics_optimized_models = pd.DataFrame()
# for i in range(len(scoring_metrics)):
#     print(scoring_metrics_names[i])
#     df_results = run_multiple_cv_runs(best_model_list, names_model, X_selected, y, 10, 10,
#                                                        cv_metric=scoring_metrics[i])
#     df_results["scoring_metric"] = scoring_metrics_names[i]
#     # add to df_results_all_metrics
#     df_results_all_metrics_optimized_models = pd.concat([df_results_all_metrics_optimized_models, df_results])
#
# # save the df_results_all_metrics as csv, so I do not have to run the above code again
# df_results_all_metrics_optimized_models.to_csv("../data/df_results_all_metrics_optimized_models.csv", index=False)











# read in the df_results_all_metrics
df_results_all_metrics_optimized_models = pd.read_csv("../data/df_results_all_metrics_optimized_models.csv", header=0)


# drop recall mic <= 5 from df_results_all_metrics
df_results_all_metrics_optimized_models_pl = df_results_all_metrics_optimized_models[
    df_results_all_metrics_optimized_models["scoring_metric"] != "Recall (MIC <= 5)"]
# also drop r2 score
df_results_all_metrics_optimized_models_pl = df_results_all_metrics_optimized_models_pl[
    df_results_all_metrics_optimized_models_pl["scoring_metric"] != "R2 score"]
# also drop f1 score mic <= 5
df_results_all_metrics_optimized_models_pl = df_results_all_metrics_optimized_models_pl[
    df_results_all_metrics_optimized_models_pl["scoring_metric"] != "F1 score (MIC <= 5)"]
# also drop precision <= 5
df_results_all_metrics_optimized_models_pl = df_results_all_metrics_optimized_models_pl[
    df_results_all_metrics_optimized_models_pl["scoring_metric"] != "Precision (MIC <= 5)"]


# plot the results
plot_multiple_cv_runs(df_results_all_metrics_optimized_models_pl, "10-fold cross validation results optimized models",
                      "10foldcv_results_all_metrics_optimized_models", width=15, height=6, location="upper left")
plot_multiple_cv_runs(df_results_all_metrics_optimized_models, "10-fold cross validation results optimized models",
                      "10foldcv_results_all_metrics_optimized_models_all_metrics")

# plot only within 1 dilution accuracy for all optimized models
df_results_all_metrics_optimized_models_pl_1tier = df_results_all_metrics_optimized_models[
    df_results_all_metrics_optimized_models["scoring_metric"] == "within-1-tier accuracy"]
# create boxplot
plt.figure(figsize=(6, 6))
# make boxplot, color steelblue
sns.boxplot(x="Model", y="Median_Scores", data=df_results_all_metrics_optimized_models_pl_1tier, color="#4682B4")
# add points
sns.swarmplot(x="Model", y="Median_Scores", data=df_results_all_metrics_optimized_models_pl_1tier, color="black",
                size=4, edgecolor="grey", linewidth=1)
# add title
plt.title("within-1-tier accuracy of models", fontsize=16, fontweight="bold", pad=20)
# make y axis title, within 1 dilution accuracy
plt.ylabel("Within 1 tier accuracy", fontsize=15)
# remove x axis label
plt.xlabel("")
# save the plot as svg
plt.rcParams['savefig.facecolor'] = 'white'
# make x axis labels rotate
plt.xticks(rotation=45, ha="right")
# make y limits 0.73-0.79
plt.ylim(0.73, 0.79)
# select ticks on y axis
plt.yticks(np.arange(0.73, 0.79, 0.01))
plt.savefig("../analysis/10foldcv_results_optimized_models_1tier.svg", bbox_inches='tight')





# now use hyperopt to do hyperparameter tuning for the random forest model, using cross validation for every run!
# define function to maximize for hyperopt (minimize MSE)
def rf_mse_cv(params, random_state=42, cv=10, X=X_selected, y=y):
    """
    Function to optimize the random forest model using hyperopt
    :param params: parameters to optimize
    :param random_state: random state
    :param cv: number of cross validation splits
    :return: mean squared error
    """
    # set the parameters
    params = {'n_estimators': params['n_estimators'],
              'max_features': params['max_features'],
              'max_depth': params['max_depth'],
              'min_samples_split': params['min_samples_split'],
              'min_samples_leaf': params['min_samples_leaf'],
              'bootstrap': params['bootstrap']}
    # define the model
    model = RandomForestRegressor(random_state=random_state, **params)
    # define the cross validation
    kf = KFold(n_splits=cv, shuffle=True, random_state=random_state)
    # calculate the scores
    scores = cross_val_score(model, X, y, cv=kf, scoring=rmse_scorer)
    # return the mean of the scores
    return np.mean(scores)


# define the parameter space
param_space = {
    "n_estimators": hp.choice('n_estimators', np.arange(50, 500)),
    "max_features": hp.choice("max_features", ["log2", "sqrt"] + list(np.arange(0.1, 1, 0.1))),
    # max depth, has to be an integer
    "max_depth": hp.choice('max_depth', np.arange(2, 21)),
    "min_samples_split": hp.uniform("min_samples_split", 0, 1),
    "min_samples_leaf": hp.choice('min_samples_leaf', np.arange(1, 10)),
    "bootstrap": hp.choice("bootstrap", [True, False])
}

# ok, now split data into 10 folds, run hyperopt on each fold, test model on test set for within 1 dilution accuracy
# define the number of folds
cv_splits = 10
# define the number of runs
num_runs = 50
# define the random state
random_state = 42
# define the number of iterations
n_iter = 50
# define the list to store the scores
scores = []
# define the list to store the models
models = []
# define the list to store the best parameters
best_params = []
# define the list to store the best scores
best_scores = []
# define the list to store the best models
best_models = []

# define the cross validation
kf = KFold(n_splits=cv_splits, shuffle=True, random_state=random_state)








#
# # # go through all folds
# for train_index, test_index in kf.split(X_selected):
#     # split the data
#     X_train, X_test = X_selected.iloc[train_index], X_selected.iloc[test_index]
#     y_train, y_test = y.iloc[train_index], y.iloc[test_index]
#     # trials will contain logging information
#     trials = Trials()
#     # runhyperopt
#     best = fmin(fn=rf_mse_cv,  # function to optimize
#                space=param_space,
#                algo=tpe.suggest,  # optimization algorithm, hyperotp will select its parameters automatically
#                max_evals=50,  # maximum number of iterations
#                trials=trials  # logging # fixing random state for the reproducibility
#                )
#
#     se_best = space_eval(param_space, best)
#     # test the model on the test set
#     # define the model
#     model = RandomForestRegressor(random_state=random_state, n_estimators=int(se_best["n_estimators"]),
#                                  max_features=se_best["max_features"],
#                                  max_depth=se_best["max_depth"],
#                                  min_samples_split=se_best["min_samples_split"],
#                                  min_samples_leaf=se_best["min_samples_leaf"],
#                                  bootstrap=se_best["bootstrap"])
#     # fit the model
#     model.fit(X_train, y_train)
#     # make predictions
#     y_pred = model.predict(X_test)
#     # calculate the within 1 dilution accuracy. prevent typeerror by converting y_test to array
#     y_test = np.array(y_test)
#     within_1_dil = within_1_dil_scorer(model, X_test, y_test)
#     # same with r2
#     r2 = r2_score(y_test, y_pred)
#     # append to list
#     scores.append(within_1_dil)
#     # append the model to the list
#     models.append(model)
#     # append the best score to the list
#     best_scores.append(trials.best_trial['result']['loss'])
#     print(best_scores)
#     print(scores)
#
# # go through models and get the tuned parameters for each. create a dictionary with the parameters as keys and the
# # values as lists
# tuned_params = {}
# for i in range(len(models)):
#     # get the parameters
#     params = models[i].get_params()
#     # go through the parameters and add them to the dictionary
#     for key, value in params.items():
#         if key in tuned_params:
#             tuned_params[key].append(value)
#         else:
#             tuned_params[key] = [value]
#
# # select only the features(keys) that are keys of param_space
# tuned_params = {key: value for key, value in tuned_params.items() if key in param_space.keys()}








# bbased on this, select the best parameters for each parameter based on most chosen or median
# parameters: {'bootstrap': [True, True, False, False, True, False, True, True, True, False],
# 'max_depth': [19, 10, 18, 17, 14, 17, 18, 14, 19, 16],
# 'max_features': [0.5, 0.30000000000000004, 0.5, 0.4, 0.7000000000000001, 'sqrt', 'sqrt', 0.7000000000000001, 'log2', 'sqrt'],
# 'min_samples_leaf': [4, 2, 1, 4, 6, 1, 2, 1, 1, 2],
# 'min_samples_split': [0.0020473307177855185, 0.0008738167622188075, 0.1173523294096549, 0.07464349184488561, 0.016262859765258104, 0.00792970153040256, 0.002372538526031177, 0.0037851064565116166, 0.017569201296589454, 0.004403612089165811],
# 'n_estimators': [428, 200, 308, 190, 407, 477, 178, 386, 451, 297]}

rf_hyper_tuned = RandomForestRegressor(bootstrap=True, max_depth=17, max_features="sqrt", min_samples_leaf=2,
                                       min_samples_split=0.006166656809784185, n_estimators=347, random_state=42)


rf_raw = RandomForestRegressor(random_state=42)








# # ok now run the 10 fold cv on the optimized models, including this one, as I have done before
# # go through all metrics and run 10 fold cv on all models with that metric. then put dfs into one dataframe,
# # adding a column with the model name
# df_results_all_metrics_optimized_model = pd.DataFrame()
# rf_other = RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt', min_samples_leaf=2)
#
# for i in range(len(scoring_metrics)):
#     print(scoring_metrics_names[i])
#     df_results = run_multiple_cv_runs([rf_raw, rf_hyper_tuned,  rf_other], ["RF_raw", "RF_hyperopt", "RF_grid"],
#                                       X_selected, y, 10, 10,
#                                       cv_metric=scoring_metrics[i])
#     df_results["scoring_metric"] = scoring_metrics_names[i]
#     # add to df_results_all_metrics
#     df_results_all_metrics_optimized_model = pd.concat([df_results_all_metrics_optimized_model, df_results])
#
# # save the df_results_all_metrics as csv, so I do not have to run the above code again
# df_results_all_metrics_optimized_model.to_csv("../data/df_results_all_metrics_RF_optimized.csv", index=False)






# read in the df_results_all_metrics
df_results_all_metrics_optimized_model = pd.read_csv("../data/df_results_all_metrics_RF_optimized.csv", header=0)

# drop recall mic <= 5 from df_results_all_metrics
df_results_all_metrics_optimized_model_pl = df_results_all_metrics_optimized_model[
    df_results_all_metrics_optimized_model["scoring_metric"] != "Recall (MIC <= 5)"]
# also drop r2 score
df_results_all_metrics_optimized_model_pl = df_results_all_metrics_optimized_model_pl[
    df_results_all_metrics_optimized_model_pl["scoring_metric"] != "R2 score"]
# also drop f1 score mic <= 5
df_results_all_metrics_optimized_model_pl = df_results_all_metrics_optimized_model_pl[
    df_results_all_metrics_optimized_model_pl["scoring_metric"] != "F1 score (MIC <= 5)"]

# plot the results
plot_multiple_cv_runs(df_results_all_metrics_optimized_model_pl, "10-fold cross validation results optimized RF",
                      "10foldcv_results_all_metrics_RF", width=15, height=6, location="upper left")




# step 3.4: check feature importance of the optimized models. use the RF model and run Shap analysis on it:
# define the model
rf = RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt', min_samples_leaf=2)

# fit the model
rf.fit(X_selected, y)

# save the trained random forest model as pickle file
# save the model to disk
filename = '../data/models_trained/rf_optimized_model.sav'
pickle.dump(rf, open(filename, 'wb'))


# create function plotting the shap values
# use treeexplainer to explain the predictions of the XGBoost model using shap values etc.
# use shap to explain the predictions and the feature importance of the model
# write as function to generate the 2 shap plots with any model


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
    plt.savefig("../analysis/" + file_name + "_summary_dot.svg", dpi=1200)

    return plt


# plot the shap values
plot_shap_values(rf, X_selected, y, "shap_values_rf_optimized_model")

print("start with interaction plots")

# now use SHAP to create interaction plots for the optimized RF model to see if there are any interactions between
# the features
# get the shap values
model_shap = shap.TreeExplainer(rf).shap_values(X_selected)
explainer = shap.TreeExplainer(rf)
explanation = explainer(X_selected)
shap_values = explanation.values

# get interaction values:
shap_interaction_values = explainer.shap_interaction_values(X_selected)
shap_interaction_values[0]


# lets have  look at shap depencence plots for the features.
# create plot
# create a colour map with orange as high, beige as middle and blue as low
# Define the colors for the colormap
colors = ['#4682B4', 'lightgray', '#B11226']  # Orange, Gray, Steel Blue

# Create a colormap using LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

plt.figure(figsize=(6, 6))
shap.dependence_plot("cas13_log2FC", shap_values, X_selected, interaction_index="expression_upec",
                     cmap=cmap)
# save the plot as svg
plt.savefig("../analysis/shap_dependence_plot_rf_optimized_model.svg", dpi=1200)

# create another dependence plot , the other way around
plt.figure(figsize=(6, 6))
shap.dependence_plot("expression_upec", shap_values, X_selected, interaction_index="cas13_log2FC",
                     cmap=cmap)
# save the plot as svg
plt.savefig("../analysis/shap_dependence_plot_rf_optimized_model_2.svg", dpi=1200)


# Ok, now try train the model on the entire X_imp and y dataset, and then predict the MIC values of the test set
# do this for all th models

# import test data:
tiling_pna_data = pd.read_csv("../data/tiling_test/tiling_data.tsv", delimiter="\t", header=0)
tiling_pna_data.head()
tiling_pna_data.columns
# get X_test and y_test
X_test = tiling_pna_data.drop(['MIC', 'pna_name', 'pna_name_full', 'pna_sequence', 'concentration',
                               'gene_name'#, 'PNA_molecular_weight'
                               ], axis=1)
# put X_test in same order as X_imp
X_test = X_test[selected_features]
y_test = np.log2(tiling_pna_data["MIC"])

# print the list of my improved model names
print(names_model)


def plot_predicted_vs_actual_mic_mod(model_name, y_pred, y_test, plot_name, file_name):
    # make y_test a factor variable
    y_test_new = y_test.astype('category')
    # make boxplot of the predicted MIC values for each actual MIC value
    palette = sns.color_palette("inferno", n_colors=7)
    palette.reverse()
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    sns.boxplot(x=y_test_new, y=y_pred, ax=ax, palette=palette)
    # add scatterplot of the predicted and actual MIC values
    sns.stripplot(x=y_test_new, y=y_pred, ax=ax, color="grey", dodge=True, alpha=0.7, size=7, edgecolor="white",
                  linewidth=1)
    # name labels [1.25, 2.5, 5, 10, 20]
    ax.set_xticklabels(["5", "10", "20", "<20"], fontsize=12)
    # set y axis labels as log2([1.25, 2.5, 5, 10, 20])
    ax.set_yticks(np.log2([5, 10, 20, 40]))
    # name labels [1.25, 2.5, 5, 10, 20]
    ax.set_yticklabels(["5", "10", "20", "<20"], fontsize=12)
    ax.set_xlabel('Tested MIC', fontsize=14)
    ax.set_ylabel('Predicted MIC', fontsize=14)
    ax.set_title("Predicted vs. tested MIC values", fontsize=16)


# create a function that takes in the model and the model name, fits the model to the entire X_imp and y dataset,
# and then predicts the MIC values of the test set
def fit_model_predict_test_set(model, model_names, X_imp, y, X_test, y_test):
    """
    Function that fits a model to the entire X_imp and y dataset, and then predicts the MIC values of the test set
    :param model: model to use
    :param model_name: name of the model
    :param X_imp: X data
    :param y: y data
    :param X_test: X test data
    :param y_test: y test data
    :return: predicted y values
    """
    # fit the model
    model.fit(X_imp, y)
    # make predictions on the test set
    y_pred = model.predict(X_test)
    # calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    print("The mean squared error of the {} model is: {}".format(model_name, mse))
    # calculate the R2 score
    r2 = r2_score(y_test, y_pred)
    # calculate
    print("The R2 score of the {} model is: {}".format(model_name, r2))
    # make boxplot of the predicted MIC values using the function
    plot_predicted_vs_actual_mic_mod(model_names, y_pred, y_test, "Predicted vs. tested MIC values",
                                     model_name + "_predicted_vs_actual_mic_boxplot")
    return y_pred


# do same but for a list of models, and save the plots as svg (of all plots)
def fit_models_predict_test_set(models, model_names, X_imp, X_test, y, y_test):
    """
    Function that fits a model to the entire X_imp and y dataset, and then predicts the MIC values of the test set
    :param model: model to use
    :param model_name: name of the model
    :param X_imp: X data
    :param y: y data
    :param X_test: X test data
    :param y_test: y test data
    :return: predicted y values
    """
    y_preds = []
    # create a plot (subplot), in twhich I can add the plots of the predicted vs actual MIC values
    fig, axs = plt.subplots(2, 4, figsize=(20, 10),
                            # make the space between the plots bigger
                            gridspec_kw={'hspace': 0.4, 'wspace': 0.4})

    # make boxplot of the predicted MIC values for each actual MIC value
    palette = [sns.color_palette("inferno", n_colors=10)[i] for i in [3, 5, 7, 9]]
    palette.reverse()

    # flatten the axs array
    axs = axs.flatten()
    for m in range(len(models)):
        # fit the model
        model = models[m]
        model_name = model_names[m]
        model.fit(X_imp, y)
        # make predictions on the test set
        y_pred = model.predict(X_test)
        # get the plot of the predicted vs actual MIC values and add it to the plot. as ax:
        ax = axs[m]
        # make y_test a factor variable
        y_test_new = y_test.astype('category')
        # make boxplot of the predicted MIC values for each actual MIC value
        sns.boxplot(x=y_test_new, y=y_pred, ax=ax, palette=palette)
        # add scatterplot of the predicted and actual MIC values
        sns.stripplot(x=y_test_new, y=y_pred, ax=ax, color="grey", dodge=True, alpha=0.7, size=7, edgecolor="white",
                      linewidth=1)
        sns.set_style("whitegrid")

        # name labels [5, 10, 20, <20]
        ax.set_xticklabels(["5", "10", "20", ">20"], fontsize=12)
        # set y axis labels as log2([5, 10, 20, <20])
        ax.set_yticks(np.log2([5, 10, 20, 40]))
        # name labels [5, 10, 20, <20]
        ax.set_yticklabels(["5", "10", "20", ">20"], fontsize=12)
        ax.set_xlabel('Tested MIC', fontsize=14)
        ax.set_ylabel('Predicted MIC', fontsize=14)
        ax.set_title(model_name, fontsize=16)
        # add box around plot
        plt.rcParams["axes.edgecolor"] = "grey"
        plt.rcParams["axes.linewidth"] = 1.5

        # make plot smaller
        plt.tight_layout()
        # add the predicted y values to the y_preds list
        y_preds.append(y_pred)
    # remove last plot (8th plot)
    #fig.delaxes(axs[7])
    # add a title to the entire plot
    fig.suptitle("Predicted vs. tested MIC values", fontsize=20, fontweight="bold", y=1.05)
    # add a legend to the plot
    fig.legend(loc="upper left", fontsize=12)
    plt.rcParams['savefig.facecolor'] = 'white'
    # save the plot as svg
    plt.savefig("../analysis/all_models_predicted_vs_actual_mic_boxplot.svg", bbox_inches='tight')
    return y_preds


best_model_list = [LinearRegression(), LassoLarsCV(eps=1e-05, normalize=True),
                    RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt', min_samples_leaf=2),
                    RandomForestRegressor(bootstrap=True, max_depth=17, max_features="sqrt", min_samples_leaf=2,
                                       min_samples_split=0.006166656809784185, n_estimators=347, random_state=42),
                    RandomForestRegressor(random_state=42),
                    xgb.XGBRegressor(colsample_bytree=0.8, learning_rate=0.1, max_depth=30, n_estimators=100,
                                    subsample=0.9),
                    GradientBoostingRegressor(max_depth=20, n_estimators=300, subsample=0.8),
                    lgb.LGBMRegressor(learning_rate=0.01, max_depth=10, n_estimators=300, subsample=0.7)]

# make columns of X_test same order as X_selected
X_test = X_test[selected_features]
X_selected = X_selected[selected_features]
names_model = ["Linear regression", "LassoLarsCV", "Random forest grid", "Random forest hyperopt", "Random forest raw",
               "XGBoost", "GradientBoostingRegressor", "LightGBM"]
# run function on all models
y_pred_tpot_rf = fit_models_predict_test_set(best_model_list, names_model, X_selected, X_test, y, y_test)

# predict only for the RF model
# train RF on X_selected and y, and then predict the MIC values of the test set
# define the model
rf = RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt', min_samples_leaf=2)
# fit the model
rf.fit(X_selected, y)
# make predictions on the test set
y_pred = rf.predict(X_test)

# make a table with pna names, gene names, pna sequences, predicted MIC values and actual MIC values
# create a dataframe with the predicted and actual MIC values
df_predicted_mic_RF = pd.DataFrame({"pna_name": tiling_pna_data["pna_name"],
                                    "gene_name": tiling_pna_data["gene_name"],
                                    "pna_sequence": tiling_pna_data["pna_sequence"],
                                    "predicted_mic": np.exp2(y_pred), "actual_mic": np.exp2(y_test)})
# print the dataframe
print(df_predicted_mic_RF)


# do same with tiling data with pure PNAS:
# import test data:
tiling_pna_data_pure = pd.read_csv("../data/tiling_test/tiling_data_pure_mics.tsv", delimiter="\t", header=0)
tiling_pna_data_pure.head()
# change RF_mean_cellsytems to gene_vulnerability_crispri and multiply by -1
tiling_pna_data_pure["gene_vulnerability_crispri"] = tiling_pna_data_pure["RF_mean_cellsytems"] * (-1) + 1
# drop RF_mean_cellsystems
tiling_pna_data_pure = tiling_pna_data_pure.drop("RF_mean_cellsytems", axis=1)
tiling_pna_data_pure.columns
# get X_test and y_test
X_test_pure = tiling_pna_data_pure.drop(['MIC', 'pna_name', 'pna_name_full', 'pna_sequence', 'concentration',
                                 'gene_name'#, 'PNA_molecular_weight'
                                         ], axis=1)
# put X_test in same order as X_imp
X_test_pure = X_test_pure[selected_features]
y_test_pure = np.log2(tiling_pna_data_pure["MIC"])

# train RF on X_selected and y, and then predict the MIC values of the test set
# define the model
rf = RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt', min_samples_leaf=2)
# fit the model
rf.fit(X_selected, y)
# make predictions on the test set
y_pred_pure = rf.predict(X_test_pure)

# make a table with pna names, gene names, pna sequences, predicted MIC values and actual MIC values
# create a dataframe with the predicted and actual MIC values
df_predicted_mic_pure = pd.DataFrame({"pna_name": tiling_pna_data_pure["pna_name"],
                                 "gene_name": tiling_pna_data_pure["gene_name"],
                                    "pna_sequence": tiling_pna_data_pure["pna_sequence"],
                                    "predicted_mic": np.exp2(y_pred_pure), "actual_mic": np.exp2(y_test_pure)})
# print the dataframe
print(df_predicted_mic_pure)

# calculate the within 1 dilution accuracy
within_1_dil = within_1_dil_scorer(rf, X_test_pure, y_test_pure)

# plot the predicted vs actual MIC values
# initiate plot
plt.figure(figsize=(5, 6))
# make y_test a factor variable
y_test_new = y_test_pure.astype('category')
# make boxplot of the predicted MIC values for each actual MIC value
palette = sns.color_palette("inferno", n_colors=7)
palette.reverse()
sns.boxplot(x=y_test_new, y=y_pred_pure, palette=palette)
# add scatterplot of the predicted and actual MIC values
sns.stripplot(x=y_test_new, y=y_pred_pure, color="grey", dodge=True, alpha=0.7, size=7, edgecolor="white",
                linewidth=1)
# set xticklabels to [1.25, 2.5, 5]
plt.xticks(np.log2([1, 2.1, 4.2]), ["1.25", "2.5", "5"], fontsize=12)
# set y axis labels as log2([1.25, 2.5, 5, 10, 20])
plt.yticks(np.log2([1.25, 2.5, 5, 10]), ["1.25", "2.5", "5", "10"], fontsize=12)
plt.xlabel('Tested MIC', fontsize=14)
plt.ylabel('Predicted MIC', fontsize=14)
plt.title("Predicted vs. tested MIC values", fontsize=16)
# save the plot as svg
plt.savefig("../analysis/pure_pnas_predicted_vs_actual_mic_boxplot.svg", dpi=1200)


# combine 3 differnt RF models (raw, less features, optimized) and make a boxplot of the predicted vs actual MIC values
# for each model
# combine the 3 cross validation results into one dataframe, while adding a column with the model name
# add a column with the model type
df_results_all_metrics_pl = df_results_all_metrics_pl.assign(type="Raw RF model")
df_results_all_metrics_selected_features = df_results_all_metrics_selected_features.assign(type="RF model less features")
df_results_all_metrics_optimized_model = df_results_all_metrics_optimized_model.assign(type="Optimized RF model")

df_results_all_metrics_rf = pd.concat([df_results_all_metrics_pl, df_results_all_metrics_selected_features,
                                       df_results_all_metrics_optimized_model])

# find which models are there;
df_results_all_metrics_rf["Model"].unique()
# select only the RF models
df_results_all_metrics_rf = df_results_all_metrics_rf[df_results_all_metrics_rf["Model"].isin(["Random forest", "RF_grid"])]
# make "Model" the type column
df_results_all_metrics_rf = df_results_all_metrics_rf.rename(columns={"Model": "mdl"})
df_results_all_metrics_rf = df_results_all_metrics_rf.rename(columns={"type": "Model"})

# ok now plot the results using the function
plot_multiple_cv_runs(df_results_all_metrics_rf, "10-fold cross validation results RF models",
                        "10foldcv_results_all_metrics_RF_combined")



# use model to predict all pnas of:
# rplS,glmS,glmU,lpxB,rplP,rplV

# first get the pnas of these genes
X_only_changing_genes = pna_data[pna_data["gene_name"].isin(["rplS", "glmS", "glmU", "lpxB", "rplP", "rplV"])]


# get the MIC values
y_actual = np.log2(X_only_changing_genes["MIC_UPEC"])
gene_names = X_only_changing_genes["gene_name"]
X_only_changing_genes = X_only_changing_genes[selected_features]
# rename columns of X_selected
X_selected = X_selected.rename(columns={'upec_tir_off_targets_2mm': "OT TIR 2mm", 'upec_tir_off_targets_1mm': "OT TIR 1mm",
                                        'upec_tir_off_targets_0mm': "OT TIR 0mm", 'upec_total_off_targets_2mm': "OT total 2mm",
                                        'upec_total_off_targets_1mm': "OT total 1mm", 'upec_total_off_targets_0mm': "OT total 0mm",
                                        'purine_percentage': "Purine %", 'sc_bases': "Self comp. bases", 'gene_vulnerability_crispri': "Gene vulnerability",
                                        'expression_upec': "Expression", "MFE_UPEC": "MFE"})
X_only_changing_genes = X_only_changing_genes.rename(columns={'upec_tir_off_targets_2mm': "OT TIR 2mm", 'upec_tir_off_targets_1mm': "OT TIR 1mm",
                                        'upec_tir_off_targets_0mm': "OT TIR 0mm", 'upec_total_off_targets_2mm': "OT total 2mm",
                                        'upec_total_off_targets_1mm': "OT total 1mm", 'upec_total_off_targets_0mm': "OT total 0mm",
                                        'purine_percentage': "Purine %", 'sc_bases': "Self comp. bases", 'gene_vulnerability_crispri': "Gene vulnerability",
                                        'expression_upec': "Expression", "MFE_UPEC": "MFE"})

# round vulnerability and expression to 2 decimals
X_selected["Gene vulnerability"] = X_selected["Gene vulnerability"].round(2)
X_only_changing_genes["Gene vulnerability"] = X_only_changing_genes["Gene vulnerability"].round(2)
X_selected["Expression"] = X_selected["Expression"].round(2)
X_only_changing_genes["Expression"] = X_only_changing_genes["Expression"].round(2)



# fit random forest model to the entire dataset
# define the model
rf = RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt', min_samples_leaf=2)
# fit the model
rf.fit(X_selected, y)
# make predictions on the test set
y_pred = rf.predict(X_only_changing_genes)

# create a table with the predicted and actual MIC values
# create a dataframe with the predicted and actual MIC values
df_predicted_mic = pd.DataFrame({"gene_name": gene_names, "predicted_mic": np.exp2(y_pred), "actual_mic": np.exp2(y_actual)})

# do exp2 of the predicted mic values

# print the dataframe
print(df_predicted_mic)

# use shap to explain the predictions for each gene and PNA. create a force plot for each PNA
# get the shap values
model_shap = shap.TreeExplainer(rf).shap_values(X_only_changing_genes)
explainer = shap.TreeExplainer(rf)
explanation = explainer(X_only_changing_genes)

plt.rcParams["axes.grid"] = False
# get force plot
shap.plots.force(explanation[0], matplotlib=True, show=False, text_rotation=45)
#plt.xlim(2, 3.5)
# change tick labels
plt.xticks(np.log2([5, 10]), [ "5", "10"])
# save the plot as svg
plt.savefig("../analysis/shap_force_plot_glmS_bad.svg", dpi=500)

shap.plots.force(explanation[1], matplotlib=True, show=False, text_rotation=45)
# change x limits
plt.xlim(1, 3.5)
# change tick labels
plt.xticks(np.log2([2.5, 5, 10]), [ "2.5", "5", "10"])
# save the plot as svg
plt.savefig("../analysis/shap_force_plot_glmS_good.svg", dpi=500)


shap.save_html("../analysis/shap_force_plot_glmS_bad.html",
               shap.force_plot(explainer.expected_value, model_shap[0,:], X_selected.iloc[0,:], show=False,
                               plot_cmap='DrDb', text_rotation=0))


# create a plot for the rplW and rpsh genes
# get the pnas of these genes
X_only_sec_genes = pna_data[pna_data["gene_name"].isin(["rplW", "rpsH"])]

# get the MIC values
y_actual = np.log2(X_only_sec_genes["MIC_UPEC"])

gene_names = X_only_sec_genes["gene_name"]
X_only_sec_genes = X_only_sec_genes[selected_features]
# rename columns of X_selected
X_selected = X_selected.rename(columns={'upec_tir_off_targets_2mm': "OT TIR 2mm", 'upec_tir_off_targets_1mm': "OT TIR 1mm",
                                        'upec_tir_off_targets_0mm': "OT TIR 0mm", 'upec_total_off_targets_2mm': "OT total 2mm",
                                        'upec_total_off_targets_1mm': "OT total 1mm", 'upec_total_off_targets_0mm': "OT total 0mm",
                                        'purine_percentage': "Purine %", 'sc_bases': "Self comp. bases", 'gene_vulnerability_crispri': "Gene vulnerability",
                                        'expression_upec': "Expression", "MFE_UPEC": "MFE"})
X_only_sec_genes = X_only_sec_genes.rename(columns={'upec_tir_off_targets_2mm': "OT TIR 2mm", 'upec_tir_off_targets_1mm': "OT TIR 1mm",
                                        'upec_tir_off_targets_0mm': "OT TIR 0mm", 'upec_total_off_targets_2mm': "OT total 2mm",
                                        'upec_total_off_targets_1mm': "OT total 1mm", 'upec_total_off_targets_0mm': "OT total 0mm",
                                        'purine_percentage': "Purine %", 'sc_bases': "Self comp. bases", 'gene_vulnerability_crispri': "Gene vulnerability",
                                        'expression_upec': "Expression", "MFE_UPEC": "MFE"})
# round vulnerability and expression to 2 decimals
X_selected["Gene vulnerability"] = X_selected["Gene vulnerability"].round(2)
X_only_sec_genes["Gene vulnerability"] = X_only_sec_genes["Gene vulnerability"].round(2)
X_selected["Expression"] = X_selected["Expression"].round(2)
X_only_sec_genes["Expression"] = X_only_sec_genes["Expression"].round(2)

# fit random forest model to the entire dataset
# define the model
rf = RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt', min_samples_leaf=2)
# fit the model
rf.fit(X_selected, y)
# make predictions on the test set
y_pred = rf.predict(X_only_sec_genes)

# create a table with the predicted and actual MIC values
# create a dataframe with the predicted and actual MIC values
df_predicted_mic = pd.DataFrame({"gene_name": gene_names, "predicted_mic": np.exp2(y_pred), "actual_mic": np.exp2(y_actual)})
# print the dataframe
print(df_predicted_mic)

# use shap to explain the predictions for each gene and PNA. create a force plot for each PNA
# get the shap values
model_shap = shap.TreeExplainer(rf).shap_values(X_only_sec_genes)
explainer = shap.TreeExplainer(rf)
explanation = explainer(X_only_sec_genes)

plt.rcParams["axes.grid"] = False
# get force plot
shap.plots.force(explanation[0], matplotlib=True, show=False, text_rotation=45)
#plt.xlim(2, 3.5)
# change tick labels
plt.xticks(np.log2([5, 10, 20]), [ "5", "10", ">10"])
# save the plot as svg
plt.savefig("../analysis/shap_force_plot_rplW_1_bad.svg", dpi=500)

shap.plots.force(explanation[1], matplotlib=True, show=False, text_rotation=45)
# change x limits
#plt.xlim(1, 3.5)
# change tick labels
plt.xticks(np.log2([5, 10, 20]), [ "5", "10", ">10"])
# save the plot as svg
plt.savefig("../analysis/shap_force_plot_rplW_2_bad.svg", dpi=500)

# now for rpsh
shap.plots.force(explanation[2], matplotlib=True, show=False, text_rotation=45)
# change x limits
#plt.xlim(1, 3.5)
# change tick labels
plt.xticks(np.log2([2.5, 5, 10]), [ "2.5", "5", "10"])
# save the plot as svg
plt.savefig("../analysis/shap_force_plot_rpsh_1_good.svg", dpi=500)

shap.plots.force(explanation[3], matplotlib=True, show=False, text_rotation=45)
# change x limits
#plt.xlim(1, 3.5)
# change tick labels
plt.xticks(np.log2([2.5, 5, 10]), [ "2.5", "5", "10"])
# save the plot as svg
plt.savefig("../analysis/shap_force_plot_rpsh_2_good.svg", dpi=500)

