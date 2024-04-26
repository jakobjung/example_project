"""
Title:      Test machine learning models on the PNA dataset but without expression and gene vulnerability data
Author:     Jakob Jung
Date:       22-03-2023
Details:    I do a similar analysis as in upec_ML_moddeling, but without the expression and gene vulnerability data.
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
                   "locus_tag", "mRNA_half_life", "gene_vulnerability_crispri", "expression_upec", "nr_pathways",
                   "nr_genes_operon", "downstream_genes_operon", "essential_genes_downstream", "upstream_genes_operon",
                   "essential_genes_upstream", "overlapping_genes", "string_interactions", "crispri_log2FC_rousset",
                   "crispri_log2FC_a_baum", "cas13_log2FC", "protein_half_life_elife"], axis=1)
# make sc bases relative to the length of the pna
X["sc_bases"] = X["sc_bases"] / 9

y = np.log2(pna_data["MIC_UPEC"])


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

names_model = ["Linear regression", "LassoLarsCV", "Random forest", "XGBoost", "GradientBoostingRegressor", "LightGBM"]


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



#
# # go through all metrics and run 10 fold cv on all models with that metric. then put dfs into one dataframe,
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
# df_results_all_metrics.to_csv("../data/mason_modeling/df_results_all_metrics.csv", index=False)





# read in the df_results_all_metrics
df_results_all_metrics = pd.read_csv("../data/mason_modeling/df_results_all_metrics.csv", header=0)

# now plot the results as boxplots and
# create a function of boxplot of the results, add points


def plot_multiple_cv_runs(df_results, plot_name, file_name):
    # change the style of the plot
    plt.figure(figsize=(20, 6))
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
    plt.legend(loc='lower left', fontsize=14)
    # save the plot as svg
    plt.rcParams['savefig.facecolor'] = 'white'
    plt.savefig("../analysis/mason_modeling/" + file_name + ".svg", bbox_inches='tight')

    return plt


# drop recall mic <= 5 from df_results_all_metrics
df_results_all_metrics_pl = df_results_all_metrics[df_results_all_metrics["scoring_metric"] != "Recall (MIC <= 5)"]
# also drop r2 score
df_results_all_metrics_pl = df_results_all_metrics_pl[df_results_all_metrics_pl["scoring_metric"] != "R2 score"]
# also drop f1 score mic <= 5
df_results_all_metrics_pl = df_results_all_metrics_pl[
    df_results_all_metrics_pl["scoring_metric"] != "F1 score (MIC <= 5)"]

# plot the results
plt_1tier = plot_multiple_cv_runs(df_results_all_metrics_pl, "10-fold cross validation results raw models",
                                  "10foldcv_results_all_metrics_raw_models")
# plot all metrics
plt_all_metrics = plot_multiple_cv_runs(df_results_all_metrics, "10-fold cross validation results raw models",
                                        "10foldcv_results_all_metrics_raw_models_all_metrics")

# step 3: feature selection
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
# plt.savefig("../analysis/mason_modeling/feature_importances_boruta_shap.svg", bbox_inches='tight')


# get the selected features
selected_features = ['Tm', 'sc_bases', 'PNA_molecular_weight', 'CAI', 'upec_tir_off_targets_1mm',
                     'MFE_UPEC', 'purine_percentage']

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
# df_results_all_metrics_selected_features.to_csv("../data/mason_modeling/df_results_all_metrics_selected_features.csv", index=False)
#





# read in the df_results_all_metrics
df_results_all_metrics_selected_features = pd.read_csv("../data/mason_modeling/df_results_all_metrics_selected_features.csv", header=0)


# drop recall mic <= 5 from df_results_all_metrics
df_results_all_metrics_selected_features_pl = df_results_all_metrics_selected_features[
    df_results_all_metrics_selected_features["scoring_metric"] != "Recall (MIC <= 5)"]
# also drop r2 score
df_results_all_metrics_selected_features_pl = df_results_all_metrics_selected_features_pl[
    df_results_all_metrics_selected_features_pl["scoring_metric"] != "R2 score"]
# also drop f1 score mic <= 5
df_results_all_metrics_selected_features_pl = df_results_all_metrics_selected_features_pl[
    df_results_all_metrics_selected_features_pl["scoring_metric"] != "F1 score (MIC <= 5)"]
# plot the results
plot_multiple_cv_runs(df_results_all_metrics_selected_features_pl,
                      "10-fold cross validation results raw models with selected features",
                      "10foldcv_results_all_metrics_selected_features")

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

# run grid search for all models
# import GridSearchCV
from sklearn.model_selection import GridSearchCV


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

# run the grid search for all models
best_params_list = []
best_score_list = []
best_model_list = []



#
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




# write the best model list directly here:
# [LinearRegression(), LassoLarsCV(eps=1e-05), RandomForestRegressor(max_depth=40, max_features='sqrt', min_samples_leaf=2,
#                       min_samples_split=5), XGBRegressor(max_depth=10, n_estimators=100),
#                       GradientBoostingRegressor(learning_rate=0.2, max_depth=10, n_estimators=200,
#                           subsample=0.7), LGBMRegressor(learning_rate=0.01, max_depth=10, subsample=0.7)]

best_model_list = [LinearRegression(), LassoLarsCV(eps=1e-05), RandomForestRegressor(max_depth=40, max_features='sqrt',
                                                                                     min_samples_leaf=2, min_samples_split=5),
                     xgb.XGBRegressor(max_depth=10, n_estimators=100), GradientBoostingRegressor(learning_rate=0.2,
                                                                                                  max_depth=10,
                                                                                                    n_estimators=200,
                                                                                                    subsample=0.7),
                        lgb.LGBMRegressor(learning_rate=0.01, max_depth=10, subsample=0.7)]


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
plt.savefig("../analysis/mason_modeling/10foldcv_results_raw_models_vs_optimized_models.svg", bbox_inches='tight')



# ok now run the 10 fold cv on the optimized models and get the results, as I have done before
# go through all metrics and run 10 fold cv on all models with that metric. then put dfs into one dataframe,
# adding a column with the model name



# df_results_all_metrics_optimized_models = pd.DataFrame()
# for i in range(len(scoring_metrics)):
#    print(scoring_metrics_names[i])
#    df_results = run_multiple_cv_runs(best_model_list, names_model, X_selected, y, 10, 10,
#                                                        cv_metric=scoring_metrics[i])
#    df_results["scoring_metric"] = scoring_metrics_names[i]
# # add to df_results_all_metrics
#    df_results_all_metrics_optimized_models = pd.concat([df_results_all_metrics_optimized_models, df_results])
#
# # save the df_results_all_metrics as csv, so I do not have to run the above code again
# df_results_all_metrics_optimized_models.to_csv("../data/mason_modeling/df_results_all_metrics_optimized_models.csv", index=False)






# read in the df_results_all_metrics
df_results_all_metrics_optimized_models = pd.read_csv("../data/mason_modeling/df_results_all_metrics_optimized_models.csv", header=0)


# drop recall mic <= 5 from df_results_all_metrics
df_results_all_metrics_optimized_models_pl = df_results_all_metrics_optimized_models[
    df_results_all_metrics_optimized_models["scoring_metric"] != "Recall (MIC <= 5)"]
# also drop r2 score
df_results_all_metrics_optimized_models_pl = df_results_all_metrics_optimized_models_pl[
    df_results_all_metrics_optimized_models_pl["scoring_metric"] != "R2 score"]
# plot the results
plot_multiple_cv_runs(df_results_all_metrics_optimized_models_pl, "10-fold cross validation results optimized models",
                      "10foldcv_results_all_metrics_optimized_models")
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
plt.savefig("../analysis/mason_modeling/10foldcv_results_optimized_models_1tier.svg", bbox_inches='tight')



# define the best performing model (RF!)
# rf = RandomForestRegressor(bootstrap=False, max_depth=20, max_features='sqrt', min_samples_leaf=2)
rf = RandomForestRegressor(max_depth=40, max_features='sqrt', min_samples_leaf=2, min_samples_split=5)
# fit the model
rf.fit(X_selected, y)

# save the trained random forest model as pickle file
# save the model to disk
filename = '../data/models_trained/rf_optimized_model_mason.sav'
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
    plt.savefig("../analysis/mason_modeling/" + file_name + "_summary_bar_mason.svg", dpi=1200)

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
    plt.savefig("../analysis/mason_modeling/" + file_name + "_summary_dot_mason.svg", dpi=1200)

    return plt


# plot the shap values
plot_shap_values(rf, X_selected, y, "shap_values_rf_optimized_model")

# get shap values for the optimized model
shap_values = shap.TreeExplainer(rf).shap_values(X_selected)

# make shap values into a dataframe
shap_values_df = pd.DataFrame(shap_values, columns=X_selected.columns)

# save them as csv
shap_values_df.to_csv("../data/mason_modeling/shap_values_rf_optimized_model.csv", index=False)

# save the raw X data as csv
X_selected.to_csv("../data/mason_modeling/X_selected.csv", index=False)





