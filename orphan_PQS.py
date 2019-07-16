import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from math import floor

gquad_loci = np.loadtxt("/Users/iamqoqao/workspace/GQ_correlations/chr21_loci", delimiter = ",")
gencode = np.genfromtxt("/Users/iamqoqao/workspace/GQ_correlations/expr/gencode_chr21.txt", usecols=(0,1,2), skip_header=0, dtype=None)

L = 46709983

expr = np.zeros(L, dtype = int)

for k in gencode:
	expr[k[1]-200:k[2]+200] = 1

orph = np.zeros(L, dtype = int)
orph_ind = []

nonorph = np.zeros(L, dtype = int)
nonorph_ind = []

for i in gquad_loci:
	if expr[int(i)]==0:
		orph[int(i)] = 1
		orph_ind.append(int(i))
	else:
		nonorph[int(i)] = 1
		nonorph_ind.append(int(i))

g2 = []
i = 0

for j in range(len(orph_ind)):
	s = 0
	print("{} out of {}".format(i, len(orph_ind)))
	for r_n in orph_ind[j:]:
		s += orph[orph_ind[j]]*orph[r_n]
	g2.append(s)
	i += 1

print("Total PQS: {}, orphaned PQS: {}".format(len(gquad_loci),len(orph_ind)))

N = 10000
K = int(floor(L/N))
barso = []
barsn = []

for i in range(K-1):
	print("{} out of {}".format(i, K-1))
	countso = Counter(orph[i*N:(i+1)*N])
	barso.append(countso[1]/10000)
	countsn = Counter(nonorph[i*N:(i+1)*N])
	barsn.append(countsn[1]/10000)

countso = Counter(orph[K*N:])
barso.append(countso[1]/10000)
countsn = Counter(nonorph[K*N:])
barsn.append(countsn[1]/10000)

# Creates two subplots and unpacks the output array immediately
f, (ax1, ax2) = plt.subplots(2)
ax1.bar(range(K), barsn, color='green')
ax2.bar(range(K), barso, color='red')

plt.savefig('orphan.eps', format='eps')
plt.show()

#plt.scatter(orph_ind, g2, marker='o', s=1)
#plt.show()

#print(expr)
#np.savetxt("expr_test", expr)

#expr_binary = [0 if i == 0 else 1 for i in expr]

