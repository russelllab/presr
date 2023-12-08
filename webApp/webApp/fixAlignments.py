#inFile = 'spec_para.fasta'
#outFile = 'spec_para2.fasta'

inFile = 'ortho.fasta'
outFile = 'ortho2.fasta'

dic = {}
num = 0
for line in open(inFile, 'r'):
    if line[0]=='>':
        name = line
        if num == 0:
            stxbp1 = name
        num += 1
        dic[name] = ''
    else:
        dic[name] += line.replace('\n', '')

#print (dic)

l = stxbp1

positions = []
length = len(dic[stxbp1])
for i in range(length):
    if dic[stxbp1][i] == '-':
        positions.append(i)
    else:
        l += dic[stxbp1][i]

l += '\n'

for name in dic:
    if name != stxbp1:
        l += name
        for i in range(len(dic[name])):
            if i not in positions:
                l += dic[name][i]

        l += '\n'

open(outFile, 'w').write(l)
