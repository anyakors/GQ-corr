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

# classify g-quad by loop length
def classify_loop(loop):
	if all([l<=3 for l in loop]):
		loop_type = 0
	elif all([l>3 for l in loop]):
		loop_type = 1
	else:
		loop_type = 2
	return loop_type

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
	loop_types = []
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
				#print("Loop is classified as {}".format(classify_loop([x-1 for x in loop])))
				pqs += 1
				loop_types.append(classify_loop([x-1 for x in loop]))
				gquad_loci.append(occurences[i])
		except IndexError:
			print("indexError")
			break
	return dict(gquad_loci=gquad_loci, loop_types=loop_types)

t0 = time.clock()

fasta_sequences = SeqIO.parse(open('/home/mookse/workspace/bioinfo/genome/hg38_by_chr/chr21'),'fasta')

for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)

L = len(sequence)

occurences = find_occurrences(sequence, 'G')

res = find_gquad(occurences)
gquad_loci = res["gquad_loci"]
loop_types = res["loop_types"]

loop0 = loop_types.count(0)
loop1 = loop_types.count(1)
loop2 = loop_types.count(2)

N = len(gquad_loci)

print("{} PQS in {} bp with loop0 {} loop1 {} loop2 {}".format(N, L, loop0, loop1, loop2))