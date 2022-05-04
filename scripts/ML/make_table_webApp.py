#!/usr/bin/python3

import os, sys, gzip
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
        self.cons_excl_para = None
        self.cons_spec_para = None
        #mutation-based features
        self.blosum = None
        self.cons_orth = None
        self.cons_para = None
        self.phi_psi = None
        self.sec = None
        self.bur = None
        self.acc = None
        self.mechismo_intramolecular = None

dic = {}

exclude = []
for line in open('../../data/test.tsv', 'r'):
    mutation = line.split()[0]
    if mutation in dic:
        exclude.append(mutation)
        dic[mutation].label = 0 if 'Population' in line else 1
exclude = list(set(exclude))

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

for line in open('../../data/210702_total_predictiontools_FPFN.csv', 'r'):
    if line.split(',')[0] != 'Category':
        mutation = line.split(',')[5]
        if mutation not in exclude:
            l = 1 if line.split(',')[0]=='Disease' else 0
            if l == 1:
                if mutation not in dic:
                    dic[mutation] = protein(mutation)
                dic[mutation].label = l

for line in open('/net/home.isilon/ag-russell/bq_rrussell/jobs/Munc18-1/STXBP1_uniprot_disease.txt', 'r'):
    if line[0] != '#':
        mutation = line.split('/')[1].split()[0]
        if mutation not in exclude:
            if mutation not in dic:
                dic[mutation] = protein(mutation)
            dic[mutation].label = 1

count = []

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

'''
for line in open('/net/home.isilon/ag-russell/bq_rrussell/jobs/Munc18-1/STXBP1_para.txt', 'r'):
    if line[0] != '#':
        mutation = line.split('/')[1].split()[0]
        value = float(line.split()[4])
        missing = line.split()[1]
        if missing != '?':
            #print (mutation, value)
            if mutation not in dic:
                dic[mutation] = protein(mutation)
            dic[mutation].cons_para = value
'''

for line in gzip.open('../../data/newEvolutionaryScores/P61764_spec_para.scores.txt.gz', 'rt'):
    mutation = line.split('/')[1].split()[0]
    value = float(line.split()[4])
    missing = line.split()[1]
    if missing != '?':
        #print (mutation, value)
        if mutation not in dic:
            dic[mutation] = protein(mutation)
        dic[mutation].cons_spec_para = value
'''
for line in open('/net/home.isilon/ag-russell/bq_rrussell/jobs/Munc18-1/STXBP1_ortho.txt', 'r'):
    if line[0] != '#':
        mutation = line.split('/')[1].split()[0]
        value = float(line.split()[4])
        missing = line.split()[1]
        if missing != '?':
            #print (mutation, value)
            if mutation not in dic:
                dic[mutation] = protein(mutation)
            dic[mutation].cons_orth = value
'''

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

data = []
for mutation in dic:
    if dic[mutation].cons_para != None:
        row = []
        row.append(mutation)
        row.append(dic[mutation].cons_para)
        if dic[mutation].label == 1:
            row.append('Diseased')
        elif dic[mutation].label == 0:
            row.append('Neutral')
        else:
            row.append('Unassigned')
        data.append(row)

print (data)
df = pd.DataFrame(data, columns = ['Mut', 'Score', 'Label'])
df.to_csv('para.tsv', sep='\t')

data = []
for mutation in dic:
    if dic[mutation].cons_spec_para != None:
        row = []
        row.append(mutation)
        row.append(dic[mutation].cons_spec_para)
        if dic[mutation].label == 1:
            row.append('Diseased')
        elif dic[mutation].label == 0:
            row.append('Neutral')
        else:
            row.append('Unassigned')
        data.append(row)

print (data)
df = pd.DataFrame(data, columns = ['Mut', 'Score', 'Label'])
df.to_csv('spec_para.tsv', sep='\t')

data = []
for mutation in dic:
    if dic[mutation].cons_orth != None:
        row = []
        row.append(mutation)
        row.append(dic[mutation].cons_orth)
        if dic[mutation].label == 1:
            row.append('Diseased')
        elif dic[mutation].label == 0:
            row.append('Neutral')
        else:
            row.append('Unassigned')
        data.append(row)

print (data)
df = pd.DataFrame(data, columns = ['Mut', 'Score', 'Label'])
df.to_csv('ortho.tsv', sep='\t')

def add_position_based_features(d, i):
    for mutation in dic:
        position = int(mutation[1:-1])
        if position in d:
            if i == 0:
                dic[mutation].side_chain_contacts = d[position]
            elif i == 1:
                print (mutation)
                dic[mutation].side_chain_residues = d[position]
            elif i == 2:
                if position in d:
                    dic[mutation].dynamine = d[position]
            elif i == 3:
                dic[mutation].iupredl = d[position]

dynamine={}; count=1
for line in open('../../data/STXBP1_dynamine.txt', 'r'):
    if line[0] not in ['*', '\n']:
        position = str(count)
        #print (count, line.split())
        dynamine[position] = float(line.split()[1])
        count+=1

add_position_based_features(dynamine, 2)
print (dynamine)

data = []
for mutation in dic:
    if dic[mutation].dynamine != None:
        row = []
        row.append(mutation)
        row.append(dic[mutation].dynamine)
        if dic[mutation].label == 1:
            row.append('Diseased')
        elif dic[mutation].label == 0:
            row.append('Neutral')
        else:
            row.append('Unassigned')
        data.append(row)

print (data)
df = pd.DataFrame(data, columns = ['Mut', 'Score', 'Label'])
df.to_csv('dyna.tsv', sep='\t')

ncont={}; nres={}

for line in gzip.open('../../data/AF/AF-P61764-F1-model_v1.mech_intra.gz', 'rt'):
    if line.split()[0]=='CO':
        position = int(line.split()[1])
        ncont[position] = int(line.replace('\n', '').split()[3])
        nres[position] = int(line.replace('\n', '').split()[4])

#add_position_based_features(ncont, 0)
add_position_based_features(nres, 1)

data = []
for mutation in dic:
    if dic[mutation].side_chain_residues != None:
        row = []
        row.append(mutation)
        row.append(dic[mutation].side_chain_residues)
        if dic[mutation].label == 1:
            row.append('Diseased')
        elif dic[mutation].label == 0:
            row.append('Neutral')
        else:
            row.append('Unassigned')
        data.append(row)

print (data)
df = pd.DataFrame(data, columns = ['Mut', 'Score', 'Label'])
df.to_csv('side_chain_residues.tsv', sep='\t')

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

data = []
for mutation in dic:
    if dic[mutation].phi_psi != None:
        row = []
        row.append(mutation)
        row.append(dic[mutation].phi_psi)
        if dic[mutation].label == 1:
            row.append('Diseased')
        elif dic[mutation].label == 0:
            row.append('Neutral')
        else:
            row.append('Unassigned')
        data.append(row)

print (data)
df = pd.DataFrame(data, columns = ['Mut', 'Score', 'Label'])
df.to_csv('phi_psi.tsv', sep='\t')
