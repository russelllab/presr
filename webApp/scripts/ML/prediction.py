import os, sys
from sklearn.datasets import load_iris
from sklearn.svm import LinearSVC
from sklearn.metrics import r2_score, roc_curve
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import confusion_matrix, matthews_corrcoef, f1_score, precision_score, recall_score
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pylab

def main(dic, clf, sfm, scaler):
    X = []
    name = []
    for mutation in dic:
        #if dic[mutation].label == None:
        if True:
            #print (mutation, dic[mutation].cons_orth)
            row = []
            ## Mutation-based features
            row.append(dic[mutation].phi_psi)
            row.append(dic[mutation].sec)
            row.append(dic[mutation].bur)
            row.append(dic[mutation].acc)
            #row.append(dic[mutation].naccess)
            row.append(dic[mutation].mechismo_intramolecular)
            row.append(dic[mutation].side_chain_contacts)
            row.append(dic[mutation].side_chain_residues)
            #row.append(dic[mutation].blosum)
            #row.append(dic[mutation].rob_iupred)
            row.append(dic[mutation].cons_orth)
            #row.append(dic[mutation].cons_para)
            row.append(dic[mutation].cons_spec_para)
            row.append(dic[mutation].cons_excl_para)
            ## position-based features
            #print (mutation, dic[mutation].side_chain_contacts)
            #row.append(dic[mutation].interface)
            #row.append(dic[mutation].interface_distance)
            row.append(dic[mutation].dynamine)
            row.append(dic[mutation].iupredl)
            row.append(dic[mutation].ptm)

            #row.append(dic[mutation].dssp)
            #print (row)
            if None not in row:
                #print (mutation, dic[mutation].label, row)
                X.append(row)
                name.append(mutation)

    print (len(X))
    X = scaler.transform(X)
    X = sfm.transform(X)
    #y_pred = clf.predict(X)
    y_pred = clf.predict_proba(X)[:,1]
    dic_predict = {}
    l = 'WT\tPOS\tMUT\tPRED\n'
    for mutation, pred in zip(name, y_pred):
        position = int(mutation[1:-1])
        mut_aa = mutation[-1]
        wt_aa = mutation[0]
        #print (mutation, mut_aa, wt_aa)
        if position not in dic_predict:
            dic_predict[position] = {}
        dic_predict[position][mut_aa] = float(pred)
        dic_predict[position]['wt'] = wt_aa
        l += wt_aa + '\t' + str(position)+ '\t' + mut_aa+ '\t' + str(pred) + '\n'

    open('web_data.tsv', 'w').write(l)


    #for mutation, pred in zip(dic, y_pred):
    for mutation in name:
        position = int(mutation[1:-1])
        mut_aa = mutation[-1]
        wt_aa = mutation[0]
        #print (mutation, mut_aa, wt_aa, pred)
        #print (mutation, mut_aa, wt_aa, dic_predict[position])
        #dic_predict[position][mut_aa] = float(dic_predict[position][mut_aa]) - float(dic_predict[position][wt_aa])
        dic_predict[position][mut_aa] = float(dic_predict[position][mut_aa])


    AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    l = 'WT\tPos\t' + '\t'.join(AA) + '\n'
    data = []
    positions = []
    for position in range(1, 595):
        if position in dic_predict:
            row = []
            positions.append(dic_predict[position]['wt']+str(position))
            l += dic_predict[position]['wt']+'\t'+str(position)
            #row.append(i)
            for aa in AA:
                row.append(dic_predict[position][aa])
                l+='\t'+str(round(dic_predict[position][aa],2))
            l += '\n'
            data.append(row)

    #print (data)
    #print (l)
    open('all_predictions.tsv', 'w').write(l)
    fig=pylab.figure()
    ax = fig.add_axes([0.05,0.05,0.9,0.9])
    sns.heatmap(data, xticklabels=AA, yticklabels=positions, vmin = 0.0, vmax = 1.0, cmap='Blues')
    ax.set_xlabel('Amino Acids')
    ax.set_ylabel('Positions')
    ax.tick_params(axis='y', labelsize=5, left=True, labelleft=True)
    ax.grid(True, linewidth=0.25, color='grey', linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    #plt.savefig('all_predictions.svg')
    plt.savefig('all_predictions_proba.svg')
    #plt.show()
