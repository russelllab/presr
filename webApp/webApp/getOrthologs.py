#!/usr/bin/python3

lst = []
for line in open('P61764_bpsh.aln', 'r'):
    if len(line.split()) == 2:
        id = line.split()[0]
        if id != 'HUMAN_P61764':
            lst.append(id)

for line in open('P61764_excl_para.aln', 'r'):
    if len(line.split()) == 2:
        id = line.split()[0]
        if id != 'HUMAN_P61764':
            lst.append(id)

for line in open('P61764_spec_para.aln', 'r'):
    if len(line.split()) == 2:
        id = line.split()[0]
        if id != 'HUMAN_P61764':
            lst.append(id)

lst = list(set(lst))

l = ''
for line in open('P61764_all_homs.aln', 'r'):
    if len(line.split()) == 2:
        id = line.split()[0]
        if id not in lst:
            l += line
    else:
        l += line

#print (l)
open('P61764_ortho.aln', 'w').write(l)
