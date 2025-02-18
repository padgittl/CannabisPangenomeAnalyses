import sys,re,os, math, statistics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

###############
# SUBROUTINES #
###############

def readGFF(gffFile):
    geneData = {}
    with open(gffFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                scaffoldID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                if feature == 'mRNA':
                    getTranscriptID = re.search('ID=(.+);Parent',attribute)
                    transcriptID = getTranscriptID.group(1)
                    getGeneID = re.search('(.+)\.t\d+',transcriptID)
                    geneID = getGeneID.group(1)
                    #print(transcriptID,geneID)
                    if transcriptID not in geneData:
                        geneData[transcriptID] = (int(start),int(end))
    return(geneData)


# geneID  mRNAStart       mRNAStop        uniprotGeneID   uniprotDescription      percentIdentity eValue  bitScore        queryCoverage
def readHitFile(hitFile):
    hitDict = {}
    with open(hitFile,'r') as F:
        for line in F:
            if 'mRNAStart' not in line:
                geneID,startPos,stopPos,uniprotID,functionDescription,percIdentity,eValue,bitScore,queryCov = line.strip().split('\t')
                if geneID not in hitDict:
                    hitDict[geneID] = functionDescription
    return(hitDict)


def readExpressionFile(expressionFile,geneData,assemblyID1,assemblyID2,tpmThreshold):
    sampleIDList = []
    tpmCountDict = {}
    geneIDList = []
    with open(expressionFile,'r') as F:
        for line in F:
            if '\"\"' in line:
                line = line.strip().split(',')
                for sampleID in line[1:]:
                    sampleID = sampleID.split('_')
                    if 'Csat' in sampleID[0]:
                        #print(sampleID)
                        sampleID = ' '.join(sampleID[2:])
                    else:
                        sampleID = ' '.join(sampleID[1:])
                    if 'Shoottips' in sampleID:
                        sampleID = 'Shoot tips'
                    if '12light' in sampleID:
                        sampleID = sampleID.replace('12light','12/12\nlights')
                    sampleIDList.append(sampleID)
                    # print(sampleID)
            else:
                line = line.strip().split(',')
                geneID = line[0]
                geneIDList.append(geneID)
                # print(geneID)
                if geneID not in tpmCountDict:
                    tpmCountDict[geneID] = []
                for tpmCount in line[1:]:
                    tpmCount = float(tpmCount)
                    tpmCountDict[geneID].append(tpmCount)

    filteredData = {}
    filteredSampleIDs = {}
    for geneID in geneData:
        start,stop = geneData[geneID]
        if geneID in tpmCountDict:
            for i in range(len(tpmCountDict[geneID])):
                tpmValue = tpmCountDict[geneID][i]
                sampleID = sampleIDList[i]
                # print(sampleID)
                items = sampleID.split()
                samplePlantName = items[0]
                plantID = items[1] + '_' + items[2]
                timeID = items[3]
                sexID = items[4]
                tissueID,suffix = items[5].split('.')
                newSampleID = sampleID.replace(' ','_')
                if 'flower' in sampleID:
                    if geneID not in filteredData:
                        filteredData[geneID] = {}
                    if sexID not in filteredData[geneID]:
                        filteredData[geneID][sexID] = []
                    filteredData[geneID][sexID].append(tpmValue)
                    if geneID not in filteredSampleIDs:
                        filteredSampleIDs[geneID] = {}
                    if sexID not in filteredSampleIDs[geneID]:
                        filteredSampleIDs[geneID][sexID] = []
                    filteredSampleIDs[geneID][sexID].append(newSampleID)
    avgTPMDict = {}
    for geneID in filteredData:
        for sexID in filteredData[geneID]:
            avgTPM = float(sum(filteredData[geneID][sexID])) / len(filteredData[geneID][sexID])
            if geneID not in avgTPMDict:
                avgTPMDict[geneID] = {}
            if sexID not in avgTPMDict[geneID]:
                avgTPMDict[geneID][sexID] = avgTPM
                
    return(tpmCountDict,sampleIDList,filteredData,filteredSampleIDs,avgTPMDict)


def assessDiffExp(geneData,filteredData,filteredSampleIDs,avgTPMDict):
    absTPMThreshold = 0.1
    tpmDifferenceThreshold = 5
    diffExpGenes = {}
    OUT_ALL_AVG = open('female_vs_male_all_average_exp.txt','w')
    OUT_ALL_AVG.write("geneID\tfemaleAvgTPM\tmaleAvgTPM\tabsDiff\n")
    print('female_vs_male_all_average_exp.txt')
    for geneID in geneData:
        if 'chr' in geneID:
            if geneID in filteredData:
                femaleAvgTPM = avgTPMDict[geneID]['female']
                maleAvgTPM = avgTPMDict[geneID]['male']
                absDiff = abs(femaleAvgTPM-maleAvgTPM)
                OUT_ALL_AVG.write("%s\t%s\t%s\t%s\n" % (geneID,femaleAvgTPM,maleAvgTPM,absDiff))
                if maleAvgTPM > femaleAvgTPM:
                    if absDiff >= tpmDifferenceThreshold:
                        # print(geneID,femaleAvgTPM,maleAvgTPM,absDiff,annotation)
                        if geneID not in diffExpGenes:
                            diffExpGenes[geneID] = {}
                        if 'male' not in diffExpGenes[geneID]:
                            diffExpGenes[geneID]['male'] = maleAvgTPM
                        if 'female' not in diffExpGenes[geneID]:
                            diffExpGenes[geneID]['female'] = femaleAvgTPM

    lowExpDiff = {}
    for geneID in geneData:
        if 'chr' in geneID:
            if geneID not in diffExpGenes:
                femaleAvgTPM = avgTPMDict[geneID]['female']
                maleAvgTPM = avgTPMDict[geneID]['male']
                femaleTPMSet = set(filteredData[geneID]['female'])
                maleTPMSet = set(filteredData[geneID]['male'])
                if len(femaleTPMSet) == 1:
                    if 0.0 not in femaleTPMSet:
                        print(geneID,femaleTPMSet)
                        sys.exit()
                    else:
                        if min(filteredData[geneID]['male']) > absTPMThreshold:
                            if geneID not in lowExpDiff:
                                lowExpDiff[geneID] = {}
                            if 'male' not in lowExpDiff[geneID]:
                                lowExpDiff[geneID]['male'] = maleAvgTPM
                            if 'female' not in lowExpDiff[geneID]:
                                lowExpDiff[geneID]['female'] = femaleAvgTPM

    OUT1 = open('male_vs_female_gene_exp_diffs.txt','w')
    OUT1.write("labelID,geneID,start,stop,femaleID,femaleAvgTPM,maleID,maleAvgTPM,tpmDiff,annotation\n")
    print("male_vs_female_gene_exp_diffs.txt")
    for geneID in lowExpDiff:
        if geneID in geneData:
            start,stop = geneData[geneID]
        else:
            print(geneID,'not in geneData')
            sys.exit()
        if geneID in summaryAnnotDict:
            annotation = summaryAnnotDict[geneID]
        else:
            if geneID in hitDict:
                annotation = hitDict[geneID]
            else:
                annotation = 'NoAnnotation'
        femaleAvgTPM = lowExpDiff[geneID]['female']
        maleAvgTPM = lowExpDiff[geneID]['male']
        diff = abs(femaleAvgTPM-maleAvgTPM)
        OUT1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("lowExpDiff",geneID,start,stop,'female',femaleAvgTPM,'male',maleAvgTPM,diff,annotation))

    OUT2 = open('male_vs_female_gene_min5TPM_exp_diffs.txt','w')
    OUT2.write("labelID,geneID,femaleID,femaleAvgTPM,maleID,maleAvgTPM,tpmDiff,annotation\n")
    print('male_vs_female_gene_min5TPM_exp_diffs.txt')
    for geneID in diffExpGenes:
        if geneID in geneData:
            start,stop = geneData[geneID]
        else:
            print(geneID,'not in geneData')
            sys.exit()
        if geneID in summaryAnnotDict:
            annotation = summaryAnnotDict[geneID]
        else:
            if geneID in hitDict:
                annotation = hitDict[geneID]
            else:
                annotation = 'NoAnnotation'
        femaleAvgTPM = diffExpGenes[geneID]['female']
        maleAvgTPM = diffExpGenes[geneID]['male']
        diff = abs(femaleAvgTPM-maleAvgTPM)
        OUT2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('min5TPMDiff',geneID,start,stop,'female',femaleAvgTPM,'male',maleAvgTPM,diff,annotation))

# https://stackoverflow.com/questions/39988048/what-is-the-origin-of-matplotlibs-symlog-a-k-a-symmetrical-log-scale
# https://support.bioconductor.org/p/102024/
def symmetric_logarithm(expValue):
    base = 2
    shift = 1
    # tested 03/30/2023 for consistency with other method
    # return(math.log(expValue+1,2))
    if expValue >= 0:
        return(math.log(expValue+shift,base)-math.log(shift,base))
    else:
        return(-math.log(-expValue+shift,base)+math.log(shift,base))


def readSummaryFile(summaryFile):
    summaryAnnotDict = {}
    with open(summaryFile,'r') as F:
        for line in F:
            if 'PERCENT_MASKED' not in line:
                line = line.strip().split('\t')
                geneID = line[0]
                transcriptID = line[1]
                annotation = line[-1]
                fractionMasked = line[8]
                fractionMasked = float(fractionMasked)
                # print(geneID,annotation)
                if transcriptID not in summaryAnnotDict:
                    summaryAnnotDict[transcriptID] = annotation
    return(summaryAnnotDict)


            
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <collinearity file> <expression file> <blast output file 1> <blast output file 2> <hit file> <summary tsv file> <gene gff file> <genome ID 1> <genome ID 2> <TPM threshold, e.g. 0.1> \n"
if len(sys.argv) != 11:
    print(usage)
    sys.exit()

collinearityFile = sys.argv[1]
expressionFile = sys.argv[2]
blast_output_file1 = sys.argv[3]
blast_output_file2 = sys.argv[4]
hitFile = sys.argv[5]
summaryFile = sys.argv[6]
gffFile = sys.argv[7]
assemblyID1 = sys.argv[8]
assemblyID2 = sys.argv[9]
tpmThreshold = sys.argv[10]

geneData = readGFF(gffFile)
summaryAnnotDict = readSummaryFile(summaryFile)
hitDict = readHitFile(hitFile)
tpmCountDict,sampleIDList,filteredData,filteredSampleIDs,avgTPMDict = readExpressionFile(expressionFile,geneData,assemblyID1,assemblyID2,tpmThreshold)
assessDiffExp(geneData,filteredData,filteredSampleIDs,avgTPMDict)


