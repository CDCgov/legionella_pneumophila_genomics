# Author : Shatavia Morrison @ cdc
# Python 2.7

#! /usr/bin/env python

import sys,os
from numpy import *
from collections import OrderedDict
from Bio import SeqIO
from Bio import Entrez
import subprocess,gzip,fnmatch,argparse



### Parameters in scripts to parse sites below 25x and consensus alignment

tempHold=[]

fastaFile = sys.argv[1]
maskInfo = sys.argv[2]
outputFile = sys.argv[3]


### the sites that need to be changed from nucleotide to N
def replaceNucs(seqString,coorList):
	seq = list(seqString)
	for v in range(0,len(coorList)):
		#print coorList[v]
		rangeLen = coorList[v][1] - coorList[v][0]
		for m in range(0, rangeLen):
			count = coorList[v][0] + m
			seq[count] ='N'
	return seq	
### Read fasta file with BioPython to read into temp directory	
sequence = list(SeqIO.parse(fastaFile, "fasta"))
for i in range(0,len(sequence)):
	tempHold.append(sequence[i].id)
	tempHold.append(sequence[i].seq)
### Read in the sites from bam file 
sitesChange=[]
with open(maskInfo,'r') as e:
	for line in e:
		coorInfo=[]
		tempInfo =  line.rstrip("\n")
		info = tempInfo.split("\t")
		start = int(info[1])
		end = int(info[2])
		depth = int(info[3])
		if depth < 26:
			coorInfo.append(start)
			coorInfo.append(end)
			sitesChange.append(coorInfo)


updateSeq = replaceNucs(tempHold[1], sitesChange)
#print updateSeq[3370453]
#print updateSeq[3393757]
#print updateSeq[3368758:3368929]


finalSeq = ''.join(updateSeq)

f = open(outputFile,"w")
f.write(">NC_002942.5"+"\n")
f.write(finalSeq)
print(outputFile)