#!/usr/bin/python3

##Header
__author__ = "Gurdeep Singh"
__credits__ = ["Gurdeep Singh"]
__license__ = "GPL"
__email__ = "gurdeep330@gmail.com"

import os, sys
from sklearn.linear_model import LassoCV
from sklearn.feature_selection import SelectFromModel
import numpy as np
import pylab
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from Bio.SeqUtils import seq1

def main(dic, protein):
    ## Label the test variants and store them in exclude
    exclude = []
    for line in open('../../data/test.tsv', 'r'):
        mutation = line.split()[0]
        if mutation in dic:
            exclude.append(mutation)
            dic[mutation].label = 0 if 'Population' in line else 1
    exclude = list(set(exclude))

    ## Label the knwon neutral variants from literature review
    for line in open('../../data/210705_neutral_mutations.csv', 'r'):
        if line.split(',')[0] != 'Chromosome':
            consequence = line.split(',')[9].split('.')[1]
            mutation = seq1(consequence[:3])+consequence[3:-3]+seq1(consequence[-3:])
            if mutation not in exclude:
                if mutation not in dic:
                    dic[mutation] = protein(mutation)
                if int(line.split(',')[-3]) >= 1:
                    dic[mutation].label = 0
                #else:
                #    dic[mutation].label = 1

    ## Label variants from literature review
    for line in open('../../data/210702_total_predictiontools_FPFN.csv', 'r'):
        if line.split(',')[0] != 'Category':
            mutation = line.split(',')[5]
            if mutation not in exclude:
                l = 1 if line.split(',')[0]=='Disease' else 0
                if l == 1:
                    if mutation not in dic:
                        dic[mutation] = protein(mutation)
                    dic[mutation].label = l

    ## Label disease variants from UniProt
    for line in open('../../data/STXBP1_uniprot_disease.txt', 'r'):
        if line[0] != '#':
            mutation = line.split('/')[1].split()[0]
            if mutation not in exclude:
                if mutation not in dic:
                    dic[mutation] = protein(mutation)
                dic[mutation].label = 1

    count = []

    ## Label disease variants from Xian et al
    xian = []
    for line in open('../../data/210901_Xian2021_missense_list.tsv', 'r'):
        if line[0] != '#':
            mutation = line.split('\t')[12]
            if mutation in dic:
                if dic[mutation].label == None or dic[mutation].label == 0:
                    #print (mutation, dic[mutation].label)
                    count.append(mutation)
                    #exclude.append(mutation)
                    dic[mutation].label = 1
            else:
                dic[mutation] = protein(mutation)
                dic[mutation].label = 1

    ## Create a train matrix
    X=[]; y=[]; X_para=[]
    names = []
    for mutation in dic:
        if mutation not in exclude:
            if dic[mutation].label in [0, 1]:
                row = []
                ## Mutation-based features
                row.append(dic[mutation].phi_psi)
                row.append(dic[mutation].sec)
                row.append(dic[mutation].bur)
                row.append(dic[mutation].acc)
                row.append(dic[mutation].mechismo_intramolecular)
                row.append(dic[mutation].side_chain_contacts)
                row.append(dic[mutation].side_chain_residues)
                row.append(dic[mutation].cons_orth)
                row.append(dic[mutation].cons_spec_para)
                row.append(dic[mutation].cons_excl_para)

                ## position-based features
                row.append(dic[mutation].dynamine)
                row.append(dic[mutation].iupredl)
                row.append(dic[mutation].ptm)

                if None not in row:
                    names.append(mutation)
                    X.append(row)
                    y.append(dic[mutation].label)
            else:
                #print ('NO label found for', mutation)
                pass

    print ('Number of positives in train set:', np.count_nonzero(y))
    print ('Number of negatives in train set:', len(y) - np.count_nonzero(y))

    ## Scale the train matrix
    scaler = MinMaxScaler()
    scaler.fit(X)
    X = scaler.transform(X)

    extra = []
    for line in open('../../data/210901_Xian2021_missense_list.tsv', 'r'):
        if line.split('\t')[12] not in names and line.split('\t')[12] not in exclude:
            extra.append(line.split('\t')[12])

    ## Save train matrix
    l = ''
    for name, label in zip(names, y):
        l += name + '\t' + str(label) + '\n'
    open('train.txt', 'w').write(l)
    #sys.exit()

    features = ['Phi_Psi', 'Sec', 'Bur', 'Acc', 'Mechismo_intr', 'SC_Cont', 'SC_Res', 'Orthologs', 'Spec_Paralogs', 'Excl_Paralogs', 'Dynamine', 'IUPredl', 'PTM']

    ## Feature selection using lasso regression
    lasso = LassoCV().fit(X, y)

    ## Select impotant features (threshold >0)
    importance = np.abs(lasso.coef_)
    threshold = 0.001
    sfm = SelectFromModel(lasso, threshold=threshold).fit(X, y)
    X = sfm.transform(X)
    features = sfm.transform([features])[0]

    ## Create test matrix from variants stored in exclude
    X_test=[]; y_test=[]; new_test=[]; X_para=[]
    #for mutation, truth in zip(test, label):
    done = []
    for mutation in exclude:
        if mutation[0] != mutation[-1] and mutation not in done:
            done.append(mutation)
            truth = dic[mutation].label
            #print (mutation, dic[mutation].cons_orth)
            row = []
            ## Mutation-based features
            row.append(dic[mutation].phi_psi)
            row.append(dic[mutation].sec)
            row.append(dic[mutation].bur)
            row.append(dic[mutation].acc)
            row.append(dic[mutation].mechismo_intramolecular)
            row.append(dic[mutation].side_chain_contacts)
            row.append(dic[mutation].side_chain_residues)
            row.append(dic[mutation].cons_orth)
            row.append(dic[mutation].cons_spec_para)
            row.append(dic[mutation].cons_excl_para)

            ## position-based features
            row.append(dic[mutation].dynamine)
            row.append(dic[mutation].iupredl)
            row.append(dic[mutation].ptm)

            if None not in row:
                new_test.append(mutation)
                X_test.append(row)
                y_test.append(truth)
            else:
                print (mutation, row)
                print (dic[mutation].mechismo_intramolecular)

    X_test = np.array(X_test)
    y_test = np.array(y_test)
    ## Scale all the features from 0 to 1
    X_test = scaler.transform(X_test)

    ## Select the feature (based on lasso)
    X_test = sfm.transform(X_test)
    print ('Number of positives in test set:', np.count_nonzero(y_test))
    print ('Number of negatives in test set:', len(y_test) - np.count_nonzero(y_test))

    return X, y, X_test, y_test, features, new_test, exclude, sfm, scaler

if __name__ == "__main__":
    print ('calling make_fm')
    main(dic, protein)
