# Author : Shatavia Morrison @ cdc 
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

from collections import OrderedDict
# Input parameters for scripts

sites = "/scripts/RF_sites_20180622.txt"
sgPos = sys.argv[1]
output = sys.argv[2]



with open(sites) as file:
	sitesPos = [line.strip() for line in file]

#print sitesPos

with open(sgPos) as f:
	for line in f:
		lines = line.rstrip('\n')
		if lines.startswith("orf"):
			m = lines.replace(":",".")
			position = m.split(",")
		else:
			info = lines.split(",")
#sgSites=OrderedDict()
sgSites={}
for i in range(0,len(position)):
#	print position[i]
	if sgSites.has_key(position[i]):
		sgSites[position[i]].append(info[i])
	else:
		sgSites[position[i]] = [info[i]]

RFsites=[]
RFinfo=[]
for m in range(0,len(sitesPos)):
	if sgSites.has_key(sitesPos[m]):
#		print sitesPos[m], sgSites[sitesPos[m]][0]
		RFsites.append(sitesPos[m])
		RFinfo.append(sgSites[sitesPos[m]][0])

#print ','.join(RFsites)
#print ','.join(RFinfo)
f = open(output,"w")
f.write(','.join(RFsites)+"\n")
f.write(','.join(RFinfo)+"\n")
f.close()

