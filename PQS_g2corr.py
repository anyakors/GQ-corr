import numpy as np
import os
import matplotlib.pyplot as plt
from Bio import SeqIO
import seaborn as sns
import time 

# enumerate all G in chromosome first
def find_occurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

def ifsequential(a, b):
	if b-a==1:
		return 1
	if b-a!=1:
		return 0

# define run of sequential G
def runs(arr):
	s = 0
	for j in range(len(arr)-1):
		if ifsequential(arr[j], arr[j+1]):
			s += 1
		else:
			break
	return s

def ifloop(a,b):
	return a-b < 8

# print iterations progress from StackOverflow
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

def find_gquad(occurences):
	# PQS definition: [G_{3+}L_{1-8}]_{4+}
	pqs = 0
	gquad_loci = []
	for i in range(len(occurences)):
		s = []
		loop = []
		j = 0
		try:
			s.append(runs(occurences[i+sum(s)+len(s):i+sum(s)+len(s)+4]))
			loop.append(occurences[i+sum(s)+len(s)]-occurences[i+sum(s)+len(s)-1])
			while s[-1]>=2 and loop[-1]<8:
				s.append(runs(occurences[i+sum(s)+len(s):i+sum(s)+len(s)+4]))
				loop.append(occurences[i+sum(s)+len(s)]-occurences[i+sum(s)+len(s)-1])
			if len(s)>=4:
				#print("PQS with G-runs length", [x+1 for x in s], "and loop length", [x-1 for x in loop], "at position", occurences[i])
				pqs += 1
				gquad_loci.append(occurences[i])
		except IndexError:
			print("indexError")
			break
	return gquad_loci

t0 = time.clock()

fasta_sequences = SeqIO.parse(open('/home/mookse/workspace/bioinfo/genome/hg38_by_chr/chr19.fa'),'fasta')

for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)

print(name)

L = len(sequence)

# just fill in an array of guanine occurences
occurences = find_occurrences(sequence, 'G')

# find and return PQS
gquad_loci = find_gquad(occurences)
# amount of PQS
N = len(gquad_loci)

print("{} PQS in {} bp".format(N, L))

# save PQS loci to a file
np.savetxt(name+'_PQS_loci', gquad_loci, delimiter=',')

# construct an occupancy array
occupancy = np.zeros(L, dtype = int)

for k in gquad_loci:
	occupancy[k] = 1

# g_2 correlation function
g2 = []

# further calculations are computationally heavy, ~O(n^2), so good to have a progress bar
printProgressBar(0, N, prefix = 'g_2 building:', suffix = 'complete', length = 99)
i = 0
for r in gquad_loci:
	s = 0
	for n in range (1, L-r):
		s += occupancy[r]*occupancy[n+r]
	g2.append(s/L)
	i += 1
	printProgressBar(i + 1, N, prefix = 'g_2 building:', suffix = 'complete', length = 99)

print('Finished in {} s'.format(time.clock() - t0))

norm = sum(g2)
g2 = [x/norm for x in g2]
np.savetxt(name+'_g2', g2, delimiter=',')

# plot the correlation function
plt.scatter(gquad_loci, g2, marker='o', s=1)
plt.xlabel('position in {}, bp'.format(name))
plt.ylim(0, max(g2))
plt.ylabel(r"$g_2(x)$")

plt.savefig(name+'_g2.eps', format='eps')
plt.show()