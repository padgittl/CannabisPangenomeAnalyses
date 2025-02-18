#!/bin/python
import sys, os, re
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
import math


###############
# SUBROUTINES #
###############


def readTSVFile(tsvFile):
    dataDict = {}
    with open(tsvFile,'r') as F:
        for line in F:
            # Sample  Name    Length  EffectiveLength TPM     NumReads
            # EH23_Early_Flower       EH23a.chr1.v1.g000010.t1        3063    2849.169        2.464884        124.266
            if 'Sample' not in line:
                sampleID,transcriptID,length,effectiveLength,TPM,numReads = line.strip().split('\t')
                TPM = float(TPM)
                # print(sampleID,transcriptID,TPM)
                if sampleID not in dataDict:
                    dataDict[sampleID] = []
                dataDict[sampleID].append((transcriptID,TPM))
    for sampleID in dataDict:
        dataDict[sampleID].sort(key=lambda x: x[0], reverse=False)
    return(dataDict)
                
# https://academic.oup.com/plphys/article/180/4/1877/6117720
# "Transcripts with TPM values lower than 5 across all varieties were removed from subsequent analysis"
def createDataStructure(dataDict,tpmThreshold,outName):
    # supplementary material for this paper --> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7737656/
    # "genes were selected based on expression thresholds of 0.1 TPM in  20% of samples"
    sampleCountThreshold = 0.2
    OUT = open(outName + "_tpm" + str(tpmThreshold) + ".csv",'w')
    OUT.write("\"\",%s\n" % (','.join(list(dataDict.keys()))))
    geneDict = {}
    for sampleID in dataDict:
        for geneID,tpm in dataDict[sampleID]:
            if geneID not in geneDict:
                geneDict[geneID] = []
            geneDict[geneID].append(tpm)
    for geneID in geneDict:
        valueList = []
        sampleCount = len(geneDict[geneID])
        minSampleCount = sampleCount * sampleCountThreshold
        for tpm in geneDict[geneID]:
            if tpm >= tpmThreshold:
                valueList.append(tpm)
        if len(valueList) >= minSampleCount:
            joinedTPM = ",".join(str(i) for i in geneDict[geneID])
            OUT.write("%s,%s\n" % (geneID,joinedTPM))
            
                        
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <tsv file> <out file name> <tpm threshold>\n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

tsvFile = sys.argv[1]
outName = sys.argv[2]
tpmThreshold = sys.argv[3]

tpmThreshold = float(tpmThreshold)

dataDict = readTSVFile(tsvFile)
createDataStructure(dataDict,tpmThreshold,outName)
