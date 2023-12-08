#To extract positions where ortholog score is low and still predictes as a Disease Variant

dic = {}
for line in open('featureTable.tsv', 'r'):
    if line.split()[0] != 'Variant':
        variant = line.split('\t')[0]
        ortholog_score = float(line.split('\t')[4])
        #print (variant, ortholog_score)
        dic[variant] = ortholog_score

predictions = {}
for line in open('all_predictions.tsv', 'r'):
    if line.split('\t')[0]!='WT':
        wt = line.split('\t')[0] + line.split('\t')[1]
        for num, aa in enumerate(AA):
            variant = wt+aa
            prob = float(line.replace('\n','').split('\t')[num+2])
            predictions[variant] = prob
    else:
        AA = line.replace('\n','').split('\t')[2:]

for variant in predictions:
    if variant[0] != variant[-1]:
        wt = variant[:-1] + variant[0]
        if predictions[variant] >= 0.5 and predictions[wt]<0.45 and dic[variant] >= 0.65:
            #print (variant, wt, dic[variant])
            print (variant)
