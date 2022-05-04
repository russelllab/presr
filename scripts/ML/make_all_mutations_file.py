import os, sys
from Bio.SeqUtils import seq1

class muts:
    def __init__(self, mutation):
        self.mutation = mutation
        self.press = '-'
        self.xian = '-'
        self.train = '-'
        self.test = '-'
        self.literature1 = '-'
        self.literature2 = '-'
        self.uniprot = '-'

dic = {}
for line in open('all_predictions.tsv', 'r'):
    if line.split('\t')[0] != 'WT':
        for value, aa in zip(line.replace('\n', '').split('\t')[2:], AA):
            mutation = line.split('\t')[0] + line.split('\t')[1] + aa
            if mutation not in dic:
                dic[mutation] = muts(mutation)
                dic[mutation].press = value
    else:
        AA = line.replace('\n', '').split('\t')[2:]

for line in open('../../data/210702_total_predictiontools_FPFN.csv', 'r'):
    if line.split(',')[0] != 'Category':
        mutation = line.split(',')[5]
        l = 1 if line.split(',')[0]=='Disease' else 0
        if mutation not in dic:
            dic[mutation] = muts(mutation)
        dic[mutation].literature1 = l

for line in open('../../data/210705_neutral_mutations.csv', 'r'):
    if line.split(',')[0] != 'Chromosome':
        consequence = line.split(',')[9].split('.')[1]
        mutation = seq1(consequence[:3])+consequence[3:-3]+seq1(consequence[-3:])
        if int(line.split(',')[-3]) >= 1:
            if mutation not in dic:
                dic[mutation] = muts(mutation)
            dic[mutation].literature2 = 0

for line in open('/net/home.isilon/ag-russell/bq_rrussell/jobs/Munc18-1/STXBP1_uniprot_disease.txt', 'r'):
    if line[0] != '#':
        mutation = line.split('/')[1].split()[0]
        dic[mutation].uniprot = 1

for line in open('../../data/210901_Xian2021_missense_list.tsv', 'r'):
    if line[0] != '#':
        mutation = line.split('\t')[12]
        if mutation not in dic:
            dic[mutation] = muts(mutation)
        dic[mutation].xian = 1

for line in open('train.txt', 'r'):
    mutation = line.split('\t')[0]
    dic[mutation].train = line.replace('\n', '').split('\t')[1]

for line in open('../../data/test.tsv', 'r'):
    mutation = line.split()[0]
    value = line.replace('\n', '').split()[1]
    dic[mutation].test = 1 if value == 'Disease' else 0

l = 'Mutation\tWT\tPosition\tVAR\tPRESS\tLit1\tLit2\tUniProt\tXian etal\tTrain\tTest\n'
for mutation in dic:
    wt = mutation[0]
    var = mutation[-1]
    position = mutation[1:-1]
    l += mutation + '\t' + wt + '\t' + str(position) + '\t' + str(var) + '\t'
    l += str(dic[mutation].press) + '\t' + str(dic[mutation].literature1) + '\t' + str(dic[mutation].literature2) + '\t'
    l += str(dic[mutation].uniprot) + '\t' +str(dic[mutation].xian) + '\t' + str(str(dic[mutation].train)) + '\t'
    l += str(dic[mutation].test) + '\n'

#print (l)
open('all_mutations.tsv', 'w').write(l)
