import os,re,sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np


def readFile(expInfoFile,labelID):
    syntenicCount = {}
    nonSyntenicCount = {}
    totalCount = {}
    with open(expInfoFile,'r') as F:
        for line in F:
            if 'femaleAvgTPM' not in line:
                expID,geneID,start,stop,syntenyID,femaleID,femaleAvgTPM,maleID,maleAvgTPM,absDiff,annotation = line.strip().split('\t')
                # print(geneID)
                items = geneID.split('.')
                genomeID = items[0]
                chrID = items[1]
                if genomeID not in totalCount:
                    totalCount[genomeID] = {}
                if geneID not in totalCount[genomeID]:
                    totalCount[genomeID][geneID] = (syntenyID,annotation)
                if syntenyID == 'NotSyntenic':
                    if genomeID not in nonSyntenicCount:
                        nonSyntenicCount[genomeID] = {}
                    if geneID not in nonSyntenicCount[genomeID]:
                        nonSyntenicCount[genomeID][geneID] = (syntenyID,annotation)
                else:
                    if genomeID not in syntenicCount:
                        syntenicCount[genomeID] = {}
                    if geneID not in syntenicCount[genomeID]:
                        syntenicCount[genomeID][geneID] = (syntenyID,annotation)

    for genomeID in totalCount:
        # print(labelID,'total',genomeID,len(totalCount[genomeID]))
        totalID = len(totalCount[genomeID])
        if genomeID in syntenicCount:
            syntenicCountID = len(syntenicCount[genomeID])
            percentSyntenic = float(syntenicCountID) / totalID * 100
            print("PERCENT SYNTENIC",labelID,genomeID,totalID,syntenicCountID,round(percentSyntenic,2))
    return(syntenicCount,nonSyntenicCount,totalCount)


def compileData(balancedData,biasedData,onlyData,tissueID,dataDict):
    for genomeID in balancedData:
        countID = len(balancedData[genomeID])
        newLabel = genomeID + ' ' + tissueID
        if newLabel not in dataDict:
            dataDict[newLabel] = []
        dataDict[newLabel].append(countID)
    for genomeID in biasedData:
        countID = len(biasedData[genomeID])
        newLabel = genomeID + ' ' + tissueID
        if newLabel not in dataDict:
            dataDict[newLabel] = []
        dataDict[newLabel].append(countID)
    for genomeID in onlyData:
        countID = len(onlyData[genomeID])
        newLabel = genomeID + ' ' + tissueID
        if newLabel not in dataDict:
            dataDict[newLabel] = []
        dataDict[newLabel].append(countID)
    return(dataDict)


def createPlot(dataDict):
    groups = list(dataDict.keys())
    subgroups = ['Balanced expression', 'Biased expression in male or female tissue (>= 5 TPM difference)', 'Expressed only in male or female tissue']
    fullDataList = []
    for genomeID in dataDict:
        dataList = []
        for countID in dataDict[genomeID]:
            dataList.append(countID)
        fullDataList.append(dataList)
    data = np.array(fullDataList)

    fig, ax = plt.subplots()

    bottom = np.zeros(len(groups))
    for i, subgroup in enumerate(subgroups):
        ax.bar(groups, data[:, i], bottom=bottom, label=subgroup)
        bottom += data[:, i]

    # Add labels and title
    plt.xticks(rotation=45, ha='right')
    ax.set_xlabel('Genome and tissue IDs')
    ax.set_ylabel('Number of genes')
    ax.legend(loc='lower center')
    plt.tight_layout()
    plt.savefig('male_female_barchart.png', dpi=600)
    plt.savefig('male_female_barchart.svg')
    print('male_female_barchart.png')
    plt.close() 
    '''
    plt.xlabel('Gene start positions (Mb)',size=16)
    plt.ylabel('Frequency',size=16)
    plt.legend()
    plt.tight_layout()
    plt.savefig('test.png', dpi=600)
    plt.close()
    '''

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <balanced flower> <balanced leaf> <male only flower> <male only leaf> <female only flower> <female only leaf> <male higher exp flower> <male higher exp leaf> <female higher exp flower> <female higher exp leaf> \n"
if len(sys.argv) != 11:
    print(usage)
    sys.exit()

balancedFlower = sys.argv[1]
balancedLeaf = sys.argv[2]
maleOnlyFlower = sys.argv[3]
maleOnlyLeaf = sys.argv[4]
femaleOnlyFlower = sys.argv[5]
femaleOnlyLeaf = sys.argv[6]
maleHigherExpFlower = sys.argv[7]
maleHigherExpLeaf = sys.argv[8]
femaleHigherExpFlower = sys.argv[9]
femaleHigherExpLeaf = sys.argv[10]

balancedFlower_syntenicCount,balancedFlower_nonSyntenicCount,balancedFlower_totalCount = readFile(balancedFlower,'balanced_flower')
balancedLeaf_syntenicCount,balancedLeaf_nonSyntenicCount,balancedLeaf_totalCount = readFile(balancedLeaf,'balanced_leaf')

maleOnlyFlower_syntenicCount,maleOnlyFlower_nonSyntenicCount,maleOnlyFlower_totalCount = readFile(maleOnlyFlower,'male_only_flower')
maleOnlyLeaf_syntenicCount,maleOnlyLeaf_nonSyntenicCount,maleOnlyLeaf_totalCount = readFile(maleOnlyLeaf,'male_only_leaf')

femaleOnlyFlower_syntenicCount,femaleOnlyFlower_nonSyntenicCount,femaleOnlyFlower_totalCount = readFile(femaleOnlyFlower,'female_only_flower')
femaleOnlyLeaf_syntenicCount,femaleOnlyLeaf_nonSyntenicCount,femaleOnlyLeaf_totalCount = readFile(femaleOnlyLeaf,'female_only_leaf')

maleHigherExpFlower_syntenicCount,maleHigherExpFlower_nonSyntenicCount,maleHigherExpFlower_totalCount = readFile(maleHigherExpFlower,'male_higher_exp_flower')
maleHigherExpLeaf_syntenicCount,maleHigherExpLeaf_nonSyntenicCount,maleHigherExpLeaf_totalCount = readFile(maleHigherExpLeaf,'male_higher_exp_leaf')

femaleHigherExpFlower_syntenicCount,femaleHigherExpFlower_nonSyntenicCount,femaleHigherExpFlower_totalCount = readFile(femaleHigherExpFlower,'female_higher_exp_flower')
femaleHigherExpLeaf_syntenicCount,femaleHigherExpLeaf_nonSyntenicCount,femaleHigherExpLeaf_totalCount = readFile(femaleHigherExpLeaf,'female_higher_exp_leaf')
'''
groups = ['AH3Ma Flower', 'AH3Mb Flower', 'AH3Ma Leaf', 'AH3Mb Leaf']
subgroups = ['Balanced', 'Only', 'Biased']
data = np.array([[10, 15], [20, 25], [30, 35]])
'''
dataDict = {}
dataDict = compileData(balancedFlower_totalCount, maleHigherExpFlower_totalCount, maleOnlyFlower_totalCount, 'Male Flower', dataDict)
dataDict = compileData(balancedLeaf_totalCount, maleHigherExpLeaf_totalCount, maleOnlyLeaf_totalCount, 'Male Leaf', dataDict)

dataDict = compileData(balancedFlower_totalCount, femaleHigherExpFlower_totalCount, femaleOnlyFlower_totalCount, 'Female Flower', dataDict)
dataDict = compileData(balancedLeaf_totalCount, femaleHigherExpLeaf_totalCount, femaleOnlyLeaf_totalCount, 'Female Leaf', dataDict)

createPlot(dataDict)
