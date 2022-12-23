#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to predict the MAGMA score from the Full and Age-parsed datasets, for a fixed choice of random seed rsval and fold kval.
This code can be run in parallel on several computers, choosing different pairs of values (rsval,kval).
In the present study, rsval = 1,2,...,200 and kval = 0,1,...,4.
"""

import numpy as np
import pandas as pd
import csv
from sklearn.ensemble import RandomForestRegressor # Regression algorithm exploited in Boruta feature selection
from boruta import BorutaPy # Library and function implementing Boruta feature selection
from sklearn.metrics import mean_absolute_error, mean_squared_error, mean_absolute_percentage_error, r2_score # Regression performance metrics
from sklearn.model_selection import StratifiedKFold # Library and function impementing Stratified k-fold Cross Validation
from xgboost import XGBRegressor # Library and function implementing XGBoost Regression algorithm

y = pd.read_csv('/input_files/y_SCZ.PGC3.kbp35.10.ZSTAT_PGC.core_more_genes.csv',header=None,index_col=0) 

'''
The input file y_SCZ.PGC3.kbp35.10.ZSTAT_PGC.core_more_genes.csv is a dataframe consisting of 21751 rows and 1 column; 
rows are indexed with the names of genes and the column corresponds to their MAGMA score.
'''

dataset = 'Full' # This variable should be set to 'Full' or 'Age-parsed', depending on the dataset used in the MAGMA prediction 

X_selected_features = pd.read_csv('/input_files/X_all_features_for_MAGMA_prediction_{0}_PGC.core_more_features_and_genes.csv'.format(dataset),index_col='Unnamed: 0')

'''
The instruction above imports one of the following input files:
1) X_all_features_for_MAGMA_prediction_Full_PGC.core_more_features_and_genes.csv if dataset=='Full';
2) X_all_features_for_MAGMA_prediction_Age-parsed_PGC.core_more_features_and_genes.csv if dataset=='Age-parsed'.

The input file X_all_features_for_MAGMA_prediction_Full_PGC.core_more_features_and_genes.csv is a dataframe consisting of 21751 rows and 212 columns;
rows are indexed with the names of genes;
columns correspond to the following features: 60 intrinsic gene attributes (including 21 chromosome variables in the one-hot encoding format), 3 kTotal connectivity variables, 146 KME connectivity variables, 3 median expression variables.

The input file X_all_features_for_MAGMA_prediction_Age-parsed_PGC.core_more_features_and_genes.csv is a dataframe consisting of 21751 rows and 491 columns;
rows are indexed with the names of genes;
columns correspond to the following features: 60 intrinsic gene attributes (including 21 chromosome variables in the one-hot encoding format), 11 kTotal connectivity variables, 409 KME connectivity variables, 11 median expression variables.

'''

features_seqnames = ['chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9'] # Names of the 21 chromosome columns

y_seqnames = []

for g_ID in X_selected_features.index:
    if list(X_selected_features[features_seqnames].loc[g_ID]) == [0 for feat in features_seqnames]:
        y_seqnames.append('chr1')
    else:
        for feat in features_seqnames:
            if X_selected_features.at[g_ID,feat] == 1:
                y_seqnames.append(feat)

# y_seqnames is a list of 21751 elements; each entry is a string indicating the chromosome to which the corresponding gene belongs

X_selected_features = X_selected_features[[col for col in X_selected_features.columns if col not in features_seqnames]]

'''
The instruction above removes the 21 chromnosomne columns from X_selected_features.
After this operation, the dataframe X_selected features has shape:
- (21751,191) if dataset=='Full'
- (21751,470) if dataset=='Age-parsed'.
'''

# ****** Parameters of the XGBoost Regression algorithm
n_trees_val = 100
max_depth_val = 2
n_jobs_val = 100

rs = rsval # rsval should be replaced by any value in range(1,201)

# ****** Stratified 5-fold Cross Validation
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=rs)



for k, (train, test) in enumerate(skf.split(X_selected_features.values,np.array(y_seqnames))): # creation of 5 folds stratified with respect to y_seqnames
    print(k)
    if k == kval: # kval should be replaced by any value in range(0,5)
        X_train, X_test, y_train, y_test = X_selected_features.values[train], X_selected_features.values[test], y.values[train], y.values[test]
        print('Shape of the training set:',X_train.shape)
        #print(y_train.shape)
        print('Shape of the test set:',X_test.shape)
        #print(y_test.shape)
        print('\n')
        #
        X_train_df = pd.DataFrame(X_train, index=X_selected_features.index[train], columns=X_selected_features.columns)
        X_test_df = pd.DataFrame(X_test, index=X_selected_features.index[test], columns=X_selected_features.columns)
        
        # Random Forest Regression algorithm exploited in Boruta feature selection
        estimator_forest = RandomForestRegressor(
           n_jobs = -1,
           max_depth = 5,
           n_estimators = 'auto',
           random_state = rs
        )

        # Boruta feature selection algorithm
        boruta = BorutaPy(
           estimator = estimator_forest,
           n_estimators = 'auto',
           max_iter = 200,
           random_state = rs
        )

        X_train_df_with_filled_nans = X_train_df.fillna(X_train_df.mean())

        # Fit Boruta on the training set (Boruta accepts np.array, not pd.DataFrame)
        print('shape np.array(X_train_rescaled_with_filled_nans): ',np.array(X_train_df_with_filled_nans).shape)
        print('shape y_train.reshape(-1): ',y_train.reshape(-1).shape)
        boruta.fit(np.array(X_train_df_with_filled_nans),y_train.reshape(-1))

        print('Boruta feature selection. Number of selected features: ',boruta.n_features_)
        X_filtered_Boruta_train_df = X_train_df[X_train_df.columns[boruta.support_]] # Training set containing only the features selected by Boruta
        X_filtered_Boruta_test_df = X_test_df[X_test_df.columns[boruta.support_]] # Test set containing only the features selected by Boruta
        list_features_Boruta = list(X_filtered_Boruta_train_df.columns)
        print('Selected features: ')
        print(list_features_Boruta)
        
        # Save the list of features selected by Boruta algorithm
        with open('/output_files/list_features_Boruta_fold_{0}_rs_{1}_{2}.csv'.format(k,rs,dataset),'w') as features_file:
            wr = csv.writer(features_file, dialect='excel')
            wr.writerows([list_features_Boruta])
        
        # XGBoost Regression algorithm        
        model_xgb_reg = XGBRegressor(num_parallel_tree=n_trees_val, max_depth=max_depth_val, importance_type='gain', n_jobs=n_jobs_val, random_state=rs, nthread=-1).fit(X_filtered_Boruta_train_df,y_train)
        y_train_pred = model_xgb_reg.predict(X_filtered_Boruta_train_df)
        y_test_pred = model_xgb_reg.predict(X_filtered_Boruta_test_df)
        print('Fold {0}'.format(k))
        print('R^2 test set: ',r2_score(list(y_test),list(y_test_pred)))
        print('R^2 train set: ',r2_score(list(y_train),list(y_train_pred)))
        print('MAE test set: ',mean_absolute_error(list(y_test),list(y_test_pred)))
        print('MAE training set: ',mean_absolute_error(list(y_train),list(y_train_pred)))
        print('RMSE test set: ',np.sqrt(mean_squared_error(list(y_test),list(y_test_pred))))
        print('RMSE training set: ',np.sqrt(mean_squared_error(list(y_train),list(y_train_pred))))
        # Save the list of MAGMA scores for genes in the test set
        np.savetxt('/output_files/MAGMA_y_test_fold_{0}_rs_{1}_Boruta_{2}_PGC.core_more_features_and_genes_skf_seqnames.csv'.format(k,rs,dataset),y_test)
        # Save the list of predicted MAGMA scores for genes in the test set 
        np.savetxt('/output_files/MAGMA_y_test_pred_fold_{0}_rs_{1}_Boruta_{2}_PGC.core_more_features_and_genes_skf_seqnames.csv'.format(k,rs,dataset),y_test_pred)
        print('\n')
