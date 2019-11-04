# Author : Shatavia Morrison Heta Desai@ cdc 
# Python 2.7
# Purpose : Retrieve whole genome sequences with accession ids
# How to : MUST RUN FROM CMD LINE
#        Run: python Biopython_Retrieve_WGS.py accessionID.txt 

#! /usr/bin/env python

from __future__ import division
import sys,os
from numpy import *
from collections import OrderedDict
# Need to set up path for BioPython
sys.path.append("C:\Python27\\Lib\\site-packages")
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
import subprocess,gzip,fnmatch,argparse


bed = sys.argv[1]
isoName = sys.argv[2]
output = sys.argv[3]


# Variable to use throughout analysis
count = 0
presenceMtrx=[]
geneSeqMtrx=[]
allowed_bases = ["A","T","C","G"]

# Calculate the masking percentage of sequence
# https://pythonforbiologists.com/counting-bases-in-a-sequence/
def count_dna(seq, allowed_bases):
	seq = seq.upper()
	total_dna_bases =0
	for base in range(0,len(allowed_bases)):
		total_dna_bases = total_dna_bases + seq.count(allowed_bases[base])
	dna_fraction = total_dna_bases / len(seq)
	return(dna_fraction * 100,seq)



# Processing code
bed = args.bedFile
isoName = args.isolateName

orgName = bed.split("_")
orfCount = 0
orfPositions=[]
### LPS L. pneumophila Philadelphia LPS biosynthesis blast coordinates
bedInfo = list(SeqIO.parse(bed,"fasta"))
for i in range(0,len(bedInfo)):
	if "CDS"in bedInfo[i].id:
		tempSeq = bedInfo[i].seq
		tempID = bedInfo[i].id.split(":")
		rangeSplit = tempID[3].split("-")
		tempStart = int(rangeSplit[0])
		tempEnd = int(rangeSplit[1])
		if tempStart > 818059 and tempStart < 853640:
			if tempEnd > 818059 and tempEnd < 853640:
				orfCount += 1
				prcntInfo = count_dna(tempSeq, allowed_bases)
				geneSeqMtrx.append(str(prcntInfo[1]))
				if prcntInfo[0] > 80.0:
					presenceMtrx.append("1")
				else:
					presenceMtrx.append("0")
				for n in range(1,len(tempSeq)+1):
					orfPositions.append("orf"+str(orfCount)+":"+str(n))

seqAlign = ''.join(geneSeqMtrx)
charSeqMtrx = list(seqAlign)
nameHold = [isoName]
seqRegionFinal = charSeqMtrx + presenceMtrx + nameHold

for i in range(0,len(presenceMtrx)):
	orfPositions.append("G"+str(i+1))

orfPositions.append("isolateName")
orfPositions.append("SgType")

###For test Data
seqRegionFinal.append("18")

f = open(output,"w")
f.write(','.join(orfPositions)+"\n")
f.write(','.join(seqRegionFinal)+"\n")
f.close()



#LPSregion = open("LPSregionGenes_"+isoName+".txt","w")
#print ','.join(orfPositions)
#print ','.join(seqRegionFinal)

#LPSregion.write(','.join(orfPositions)+"\n")
#LPSregion.write(','.join(seqRegionFinal)+"\n")
#LPSregion.close()

