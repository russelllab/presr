#!/usr/bin/python3

##Header
__author__ = "Gurdeep Singh"
__credits__ = ["Gurdeep Singh"]
__license__ = "GPL"
__email__ = "gurdeep330@gmail.com"

import os, sys
from sklearn.datasets import load_iris
from sklearn.svm import LinearSVC
from sklearn.metrics import r2_score, roc_curve
from sklearn.linear_model import LogisticRegression
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score, roc_curve
import pandas as pd
import pylab
from sklearn.metrics import confusion_matrix, matthews_corrcoef, f1_score, precision_score, recall_score

def main(X, y, X_test, y_test, features, dic, new_test, exclude):
    ## stratified CV
    skf = StratifiedKFold(n_splits=5)
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10)
    X = np.array(X)
    y = np.array(y)

    ## To perform the randomizationt test (Salzberg test), enable the this line
    #np.random.shuffle(y)

    parameters = {'C': [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0],
                'solver': ['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'],
                'penalty': ['l1', 'l2', 'elasticnet', 'none'],
                'max_iter': [100, 250, 500, 1000, 1500, 2000]
                }

    model = LogisticRegression(class_weight='balanced')
    model = GridSearchCV(model, parameters, cv=rskf, scoring='roc_auc', n_jobs=-1)
    model.fit(X, y)

    breakLine = '-'.join(['-' for i in range(0, 50)])
    print (breakLine)
    ## Best model hyper-parameters
    print ('Best model found during the CV')
    print (model.best_params_)

    clf = LogisticRegression(class_weight='balanced', max_iter=model.best_params_['max_iter'], solver=model.best_params_['solver'], C=model.best_params_['C'], penalty=model.best_params_['penalty'])

    AUC= []; MCC= []; F1=[]; PRE=[]; REC=[]; SPE=[]
    for i in range(0,10):
        skf = StratifiedKFold(n_splits=5, shuffle=True)
        rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10)
        auc = []; mcc= []; f1=[]; pre=[]; rec=[]; spe=[]
        for train_index, test_index in skf.split(X, y):
            X_train, X_validation = X[train_index], X[test_index]
            y_train, y_validation = y[train_index], y[test_index]
            clf.fit(X_train, y_train)
            tn, fp, fn, tp = confusion_matrix(y_train, model.predict(X_train)).ravel()
            #print (tn, fp, fn, tp)
            auc.append(roc_auc_score(y_validation, model.predict_proba(X_validation)[:,1]))
            mcc.append(matthews_corrcoef(y_validation, model.predict(X_validation)))
            f1.append(f1_score(y_validation, model.predict(X_validation)))
            pre.append(precision_score(y_validation, model.predict(X_validation)))
            rec.append(recall_score(y_validation, model.predict(X_validation)))
            spe.append(recall_score(y_validation, model.predict(X_validation), pos_label=0))
        AUC.append(np.mean(auc))
        MCC.append(np.mean(mcc))
        F1.append(np.mean(f1))
        PRE.append(np.mean(pre))
        SPE.append(np.mean(spe))
        REC.append(np.mean(rec))

    ## Cross-validation results
    breakLine = '-'.join(['-' for i in range(0, 50)])
    print (breakLine)
    print ('Stratified CV results')
    print ('MET', 'AUC ', 'MCC ', 'F1  ', 'PRE ', 'REC ', 'SPE')
    print ('AVG', round(np.mean(AUC),2),round(np.mean(MCC),2),round(np.mean(F1),2),round(np.mean(PRE),2),round(np.mean(REC),2),round(np.mean(SPE),2))
    print ('STD', round(np.std(AUC),2),round(np.std(MCC),2),round(np.std(F1),2),round(np.std(PRE),2),round(np.std(REC),2),round(np.std(SPE),2))
    print ('Number of disease variants in the train set:', np.count_nonzero(y))
    print ('Number of neutral variants in the train set:', len(y) - np.count_nonzero(y))

    ## Fit the best model on the data
    clf = LogisticRegression(class_weight='balanced', max_iter=model.best_params_['max_iter'], solver=model.best_params_['solver'], C=model.best_params_['C'], penalty=model.best_params_['penalty'])
    clf.fit(X,y)

    ## Display feature_relevance
    data = []
    for value, name in zip(clf.coef_[0], features):
    #for value, name in zip(clf.feature_importances_, features):
        row = []
        row.append(name)
        row.append(value)
        row.append(abs(value))
        row.append(value/abs(value))
        data.append(row)
    df = pd.DataFrame(data, columns=['Features', 'value', 'abs', 'sign'])
    df = df.sort_values(by='abs')

    ## Run the best mdoel on the test set
    y_pred = clf.predict(X_test)

    ## Display test set metrics of the best model
    breakLine = '-'.join(['-' for i in range(0, 50)])
    print (breakLine)
    print ('Best model (PRESS) results on the test')
    print ('\t'.join(['MET', 'AUC ', 'MCC ', 'F1  ', 'PRE ', 'REC ', 'SPE']))
    print ('PRESS', end='\t')
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    print (round(matthews_corrcoef(y_test, y_pred),2), end='\t')
    print (round(roc_auc_score(y_test, clf.predict_proba(X_test)[:,1]),2), end='\t')
    predictor_auc = round(roc_auc_score(y_test, clf.predict_proba(X_test)[:,1]),2)
    print (round(f1_score(y_test, y_pred),2), end='\t')
    print (round(precision_score(y_test, y_pred),2), end='\t')
    print (round(recall_score(y_test, y_pred),2), end='\t')
    print (round(recall_score(y_test, y_pred, pos_label=0),2))
    print ('Number of disease variants in the test set:', np.count_nonzero(y_test))
    print ('Number of neutral variants in the test set:', len(y_test) - np.count_nonzero(y_test))

    breakLine = '-'.join(['-' for i in range(0, 50)])
    print (breakLine)
    print ('Results of other tools on the test')
    print ('\t'.join(['MET', 'AUC ', 'MCC ', 'F1  ', 'PRE ', 'REC ', 'SPE']))
    ## Display test set metrics of other tools
    for i in [1,3,5,7,9,11,13,15,17]:
        y_truth = []; y_pred=[]; y_bin=[]
        done = []
        for line in open('../../data/new_predictors.tsv', 'r'):
            if line[0] != '#':
                if 'None' not in line:
                    mutation = line.split('\t')[0]
                    #if True:
                    if mutation in new_test:
                        if mutation in exclude:
                        #if dic[mutation].label == None:
                            #truth = 1 if line.split('\t')[1]=='Disease' else 0
                            truth = dic[mutation].label
                            predict = float(line.split('\t')[i])
                            done.append(mutation)
                            y_truth.append(truth)
                            #print (mutation, truth)
                            y_pred.append(predict)
                            y_bin.append(float(line.split('\t')[i+1]))
            else:
                print (line.split('\t')[i].replace('\n',''), end = '\t')
                predictor = line.split('\t')[i].replace('\n','')

        print (round(matthews_corrcoef(y_truth, y_bin),2), end='\t')
        print (round(roc_auc_score(y_truth, y_pred),2), end='\t')
        predictor_auc = round(roc_auc_score(y_truth, y_pred),2)
        print (round(f1_score(y_truth, y_bin),2), end='\t')
        print (round(precision_score(y_truth, y_bin),2), end='\t')
        print (round(recall_score(y_truth, y_bin),2), end='\t')
        print (round(recall_score(y_truth, y_bin, pos_label=0),2))

    print ('\n\nContact: gurdeep[dot]singh[at]bioquant[dot]uni-heidelberg[dot]de')
    #Return the best model
    return clf
