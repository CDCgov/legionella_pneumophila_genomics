# Author : Shatavia Morrison Heta Desai@ cdc
# Python 2.7

#! /usr/bin/env python

from __future__ import division
import sys,os
# Need to set up path for BioPython
import subprocess,gzip,fnmatch,argparse, shutil
# Input parameters for scripts

#parser = argparse.ArgumentParser(description='Get prediction with isolateName ', epilog="________________________________________________")
#parser.add_argument('-p','--predictFile',type=str,required=True,help="Prediction File")
#parser.add_argument('-i','--isolateName', type=str, required=True,help="Isolate Name")

#args = parser.parse_args()



# Processing code
results = sys.argv[1]
isoName = sys.argv[2]
output = sys.argv[3]
#infoPrediction = "/scripts/predictionInfo_2018_v1.txt"

with open(results,'r')as e:
        for i, line in enumerate(e):
                if i ==1:
                        #print line
                        info = line.split(",")
                        sgPredict= info[1].replace("\"","")
                        p = isoName+"\t"+sgPredict

with open(output,'a') as myfile:
        myfile.write("Individual Serogroup Prediction Accuracy and Error Rates Based on Pneumonia Response and Surveillance Laboratory Training Set – Sept 2018"+"\n")
        myfile.write("created by: S.Morrison - SMorrison@cdc.gov - 20190514"+"\n")
        myfile.write("\n")
        myfile.write("*Note: The accuracy and error rate listed above are provided as a guide and not necessarily apply to the predicted sequence."+"\n")
        myfile.write("\n")
        myfile.write("Please perform direct fluorescent antibody or slide agglutination for laboratory confirmatory results."+"\n")
        myfile.write("\n")
        myfile.write('This is not a CLIA approved test.'+"\n")
        myfile.write("\n")
        myfile.write("Serogroup (SG)  Error   Accuracy        SG Specific Sensitivity SG Specific Specificity"+"\n")
        myfile.write("1       0.00    1       100%    98.14%"+"\n")
        myfile.write("2       0.05    0.95    94.11%  99.38%"+"\n")
        myfile.write("3       0.15    0.85    84.61%  97.00%"+"\n")
        myfile.write("4       0.23    0.77    83.33%  95.83%"+"\n")
        myfile.write("5       0.38    0.62    61.11%  98.14%"+"\n")
        myfile.write("6       0.05    0.95    94.73%  94.40%"+"\n")
        myfile.write("7       0.00    1       100%    100%"+"\n")
        myfile.write("8       0.50    0.5     50.00%  91.35%"+"\n")
        myfile.write("9       0.00    1       100%    99.40%"+"\n")
        myfile.write("10      0.38    0.62    61.53%  98.20%"+"\n")
        myfile.write("12      1.00    0       0.00%   100.00%"+"\n")
        myfile.write("13      0.08    0.92    91.66%  100.00%"+"\n")
        myfile.write("14      1.00    0       0.00%   98.86%"+"\n")
        myfile.write("15      1.00    0       0.00%   100.00%"+"\n")
        myfile.write("\n")
        myfile.write("\n")
        myfile.write("IsolateName     Serogroup Pipeline Prediction"+"\n")
        myfile.write(p)
        myfile.close()
