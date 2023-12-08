#!/usr/bin/python3

##Header
__author__ = "Gurdeep Singh"
__credits__ = ["Gurdeep Singh"]
__license__ = "GPL"
__email__ = "gurdeep330@gmail.com"

import os, sys, gzip
from sklearn.metrics import r2_score, roc_curve
from sklearn.linear_model import LogisticRegression
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix, matthews_corrcoef, f1_score, precision_score, recall_score
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pylab
import pickle

class protein:
    def __init__(self, mutation):
        self.mutation = mutation
        self.label = None
        # position-based features
        self.iupredl = None
        self.dynamine = None
        self.ptm = 0
        self.side_chain_contacts = None
        self.side_chain_residues = None
        #mutation-based features
        self.blosum = None
        self.cons_orth = None
        self.cons_excl_para = None
        self.cons_spec_para = None
        self.phi_psi = None
        self.sec = None
        self.bur = None
        self.acc = None
        self.mechismo_intramolecular = None

def main():
    ## Load BLOSUM scores
    dic = {}
    for line in open('../../data/STXBP1_blosum.txt', 'r'):
        mutation = line.split('/')[1].split()[0]
        value = float(line.split()[-1])
        if mutation not in dic:
            dic[mutation] = protein(mutation)
        dic[mutation].blosum = value

    ## Load mechismo intramolecular scores
    for line in open('/net/home.isilon/ag-russell/bq_rrussell/jobs/Munc18-1/STXBP1_mech_intra.txt', 'r'):
        if line[0] != '#':
            mutation = line.split('/')[1].split()[0]
            value = float(line.split()[-2])
            missing = line.split()[-1].replace('\n', '')
            if missing != 'missing':
                #print (mutation, value, missing)
                if mutation not in dic:
                    dic[mutation] = protein(mutation)
                dic[mutation].mechismo_intramolecular = value

    ## Load Orthologs scores
    for line in gzip.open('../../data/newEvolutionaryScores/out.3.scores.txt.gz', 'rt'):
        if line[0] != '#':
            mutation = line.split('/')[1].split()[0]
            value = float(line.split()[4])
            missing = line.split()[1]
            if missing != '?':
                #print (mutation, value)
                if mutation not in dic:
                    dic[mutation] = protein(mutation)
                dic[mutation].cons_orth = value

    ## Load human specific Paralogs scores
    for line in gzip.open('../../data/newEvolutionaryScores/P61764_spec_para.scores.txt.gz', 'rt'):
        mutation = line.split('/')[1].split()[0]
        value = float(line.split()[4])
        missing = line.split()[1]
        if missing != '?':
            #print (mutation, value)
            if mutation not in dic:
                dic[mutation] = protein(mutation)
            dic[mutation].cons_spec_para = value

    ## Load human exclusive paralogs (i.e. all homologs - all orthologs) scores
    for line in gzip.open('../../data/newEvolutionaryScores/outPara.1.scores.txt.gz', 'rt'):
        mutation = line.split('/')[1].split()[0]
        value = float(line.split()[4])
        missing = line.split()[1]
        if missing != '?':
            #print (mutation, value)
            if mutation not in dic:
                dic[mutation] = protein(mutation)
            dic[mutation].cons_excl_para = value

    ## Load DSSP scores
    for line in gzip.open('../../data/AF/AF-P61764-F1-model_v1.dssp-scores.gz', 'rt'):
        if line[0] != '#':
            mutation = line.split()[2] + line.split()[0] + line.split()[10]
            #print (mutation)
            phipsi = line.split()[18]
            sec = float(line.split()[22])
            bur = float(line.split()[26])
            acc = float(line.split()[30])
            #print (mutation, value)
            if mutation not in dic:
                dic[mutation] = protein(mutation)
            dic[mutation].phi_psi = phipsi
            dic[mutation].sec = sec
            dic[mutation].bur = bur
            dic[mutation].acc = acc

    ## position-based features
    def add_position_based_features(d, i):
        for mutation in dic:
            position = str(mutation[1:-1])
            if i == 0:
                dic[mutation].side_chain_contacts = d[position]
            elif i == 1:
                dic[mutation].side_chain_residues = d[position]
            elif i == 2:
                dic[mutation].dynamine = d[position]
            elif i == 3:
                dic[mutation].iupredl = d[position]

    dic_position = {}
    ncont={}; nres={}

    ## Load side-chain/residues scores
    for line in gzip.open('../../data/AF/AF-P61764-F1-model_v1.mech_intra.gz', 'rt'):
        if line.split()[0]=='CO':
            position = str(line.split()[1])
            ncont[position] = int(line.replace('\n', '').split()[3])
            nres[position] = int(line.replace('\n', '').split()[4])

    add_position_based_features(ncont, 0)
    add_position_based_features(nres, 1)

    ## Load DynaMine scores
    dynamine={}; count=1
    for line in open('../../data/STXBP1_dynamine.txt', 'r'):
        if line[0] not in ['*', '\n']:
            position = str(count)
            dynamine[position] = float(line.split()[1])
            count+=1

    add_position_based_features(dynamine, 2)

    ## Load IUPred scores
    iupred={}
    for line in open('../../data/STXBP1_iupred.txt', 'r'):
        if line[0] not in ['#', '\n']:
            position = str(line.split()[0])
            iupred[position] = float(line.split()[2])

    add_position_based_features(iupred, 3)

    ## Load PTM scores
    go = 0
    ptm=[]
    for line in gzip.open('../../data/Phosphorylation_site_dataset.gz', 'rt'):
        if len(line.split()) > 0:
            if line.split()[0] == 'GENE':
                go = 1
            elif go == 1 and 'STXBP1' in line and 'human' in line:
                #print (line)
                mod_res = str(line.split('\t')[4])
                mod_res = mod_res.split('-')[0]
                position = int(mod_res[1:])
                ptm.append(position)

    go = 0
    for line in gzip.open('../../data/Ubiquitination_site_dataset.gz', 'rt'):
        if len(line.split()) > 0:
            if line.split()[0] == 'GENE':
                go = 1
            elif go == 1 and 'STXBP1' in line and 'human' in line:
                #print (line)
                mod_res = str(line.split('\t')[4])
                mod_res = mod_res.split('-')[0]
                position = int(mod_res[1:])
                ptm.append(position)

    for mutation in dic:
        position = int(mutation[1:-1])
        if position in ptm:
            for i in range(position-10, position+11):
                if i in ptm:
                    dic[mutation].ptm += 1
                #print (i, position, dic[mutation].ptm)

    return dic, protein

if __name__ == "__main__":
    print ('Loading features into a dictionary')
    main()
