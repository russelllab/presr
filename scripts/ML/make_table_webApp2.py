#!/usr/bin/python3
import os, sys
from sklearn.linear_model import LassoCV
from sklearn.feature_selection import SelectFromModel
import numpy as np
import pandas as pd
import pylab
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from Bio.SeqUtils import seq1

class protein:
    def __init__(self, mutation):
        self.mutation = mutation
        self.label = None
        # position-based features
        self.iupredl = None
        self.dynamine = None
        self.dssp = None
        self.naccess = None
        self.ptm = 0
        self.interface = 0
        self.interface_distance = None
        self.side_chain_contacts = None
        self.side_chain_residues = None
        #mutation-based features
        self.blosum = None
        self.cons_orth = None
        self.cons_para_wt = None
        self.cons_para_mut = None
        self.phi_psi = None
        self.sec = None
        self.bur = None
        self.acc = None
        self.mechismo_intramolecular = None

dic = {}

for line in open('/net/home.isilon/ag-russell/bq_rrussell/jobs/Munc18-1/STXBP1_para.txt', 'r'):
    if line[0] != '#':
        mutation = line.split('/')[1].split()[0]
        position = mutation[:-1]
        wt_value = float(line.split()[2])
        mut_value = float(line.split()[3])
        missing = line.split()[1]
        if missing != '?':
            #print (mutation, value)
            if mutation not in dic:
                dic[mutation] = protein(mutation)
            dic[mutation].cons_para_pos = str(position)
            dic[mutation].cons_para_wt = str(wt_value)
            dic[mutation].cons_para_mut = str(mut_value)

data = []
for mutation in dic:
    if dic[mutation].cons_para_wt != None:
        row = []
        row.append(mutation)
        row.append(dic[mutation].cons_para_pos)
        row.append(dic[mutation].cons_para_wt)
        row.append(dic[mutation].cons_para_mut)
        data.append(row)

print (data)
df = pd.DataFrame(data, columns = ['Name', 'Pos', 'WT', 'MUT'])
df.to_csv('para2.tsv', sep='\t')
