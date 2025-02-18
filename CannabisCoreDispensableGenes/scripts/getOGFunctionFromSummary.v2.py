import os,re,sys,math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import pearsonr


def readGenomeIDFile(genomeIDFile):
    pangenomeIDs = {}
    with open(genomeIDFile,'r') as F:
        for line in F:
            genomeID = line.strip()
            if genomeID not in pangenomeIDs:
                pangenomeIDs[genomeID] = 1
    return(pangenomeIDs)


def readFileList(fileList):
    dataList = []
    with open(fileList,'r') as F:
        for line in F:
            fileName = line.strip()
            dataList.append(fileName)
    return(dataList)


def readSummaryFile(summaryDataList):
    summaryAnnotDict = {}
    for summaryFile in summaryDataList:
        fileID = summaryFile.strip()
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


def readOrthogroupTSV(orthogroupTSV,pangenomeIDs,summaryAnnotDict):
    functionDict = {}
    with open(orthogroupTSV,'r') as F:
        for line in F:
            orthoGroupInfo = line.strip().split('\t')
            if not line.startswith('Orthogroup'):
                orthogroupID = orthoGroupInfo[0]
                for item in orthoGroupInfo[1:]:
                    if ',' in item:
                        item = item.strip().split(',')
                        # filter empty lines
                        item = list(filter(str.strip, item))
                        for i in item:
                            geneID = i.strip()
                            # only cannabis gene annotations will have '.'
                            if '.' in geneID:
                                cols = geneID.split('.')
                                assemblyID = cols[0]
                                if assemblyID in pangenomeIDs:
                                    if geneID in summaryAnnotDict:
                                        annotation = summaryAnnotDict[geneID]
                                    else:
                                        print('not in summaryAnnotDict',geneID)
                                        sys.exit()
                                    if orthogroupID not in functionDict:
                                        functionDict[orthogroupID] = {}
                                    if annotation not in functionDict[orthogroupID]:
                                        functionDict[orthogroupID][annotation] = 1
                    else:
                        geneID = item
                        items = geneID.split('.')
                        assemblyID = items[0]
                        if assemblyID in pangenomeIDs:
                            if geneID in summaryAnnotDict:
                                annotation = summaryAnnotDict[geneID]
                            else:
                                print('not in summaryAnnotDict',geneID)
                                sys.exit()
                            if orthogroupID not in functionDict:
                                functionDict[orthogroupID] = {}
                            if annotation not in functionDict[orthogroupID]:
                                functionDict[orthogroupID][annotation] = 1
    ANNOTS = open('orthogroups_with_annotations.txt','w')
    print('orthogroups_with_annotations.txt')
    for orthogroupID in functionDict:
        joinedAnnot = ';'.join(list(functionDict[orthogroupID].keys()))
        #print(orthogroupID,joinedAnnot)
        ANNOTS.write("%s\t%s\n" % (orthogroupID,joinedAnnot))


############
# MAIN  ####
############


usage = "Usage: " + sys.argv[0] + " <orthogroup tsv file> <summary file list> <pangenome ID file> \n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

orthogroupTSV = sys.argv[1]
summaryFileList = sys.argv[2]
genomeIDFile = sys.argv[3]

pangenomeIDs = readGenomeIDFile(genomeIDFile)

summaryDataList = readFileList(summaryFileList)
summaryAnnotDict = readSummaryFile(summaryDataList)

readOrthogroupTSV(orthogroupTSV,pangenomeIDs,summaryAnnotDict)
