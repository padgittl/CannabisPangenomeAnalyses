import os,re,sys,math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import pearsonr


# OG0028998       4       N3.HOG0044999   4       mj      AH3Ma   AH3Ma.chr6.v1.g381750.t1        0.0     True
# ### orthogroupID        ogProteinCount  populationID    genomeID        geneID  averageEntropy  annotation
def readEntropyFile(entropyFile,proteinCountMinimum):
    initDict = {}
    functionDict = {}
    initFunctionDict = {}
    orthogroupDict = {}
    with open(entropyFile,'r') as F:
        for line in F:
            if '#' not in line:
                orthogroupID,ogProteinCount,populationID,genomeID,geneID,averageEntropy,annotation = line.strip().split('\t')
                averageEntropy = float(averageEntropy)
                if orthogroupID not in orthogroupDict:
                    orthogroupDict[orthogroupID] = ogProteinCount
                
                if orthogroupID not in initDict:
                    initDict[orthogroupID] = {}
                if populationID not in initDict[orthogroupID]:
                    initDict[orthogroupID][populationID] = {}
                if averageEntropy not in initDict[orthogroupID][populationID]:
                    initDict[orthogroupID][populationID][averageEntropy] = int(ogProteinCount)

                if orthogroupID not in initFunctionDict:
                    initFunctionDict[orthogroupID] = {}
                if annotation not in initFunctionDict[orthogroupID]:
                    initFunctionDict[orthogroupID][annotation] = 1
    for orthogroupID in initFunctionDict:
        joinedAnnot = ','.join(list(initFunctionDict[orthogroupID].keys()))
        # print(joinedAnnot)
        if orthogroupID not in functionDict:
            functionDict[orthogroupID] = joinedAnnot

    entropyDict = {}
    for orthogroupID in initDict:
        for populationID in initDict[orthogroupID]:
            if len(initDict[orthogroupID][populationID]) > 1:
                print("unexpected number of entropies",orthogroupID,populationID,len(initDict[orthogroupID][populationID]))
                sys.exit()
            else:
                for averageEntropy in initDict[orthogroupID][populationID]:
                    ogProteinCount = initDict[orthogroupID][populationID][averageEntropy]
                    if ogProteinCount >= proteinCountMinimum:
                        if orthogroupID not in entropyDict:
                            entropyDict[orthogroupID] = {}
                        if populationID not in entropyDict[orthogroupID]:
                            entropyDict[orthogroupID][populationID] = []
                        entropyDict[orthogroupID][populationID].append(averageEntropy)
    return(entropyDict,orthogroupDict,functionDict,initDict)


def collectRepeatGroups(entropyDictWithTEs,entropyDictNoTEs):
    # entropyDict[orthogroupID][populationID].append(averageEntropy)
    onlyTEGroups = {}
    for orthogroupID in entropyDictWithTEs:
        if orthogroupID not in entropyDictNoTEs:
            #print(orthogroupID)
            if orthogroupID not in onlyTEGroups:
                onlyTEGroups[orthogroupID] = {}
            for populationID in entropyDictWithTEs[orthogroupID]:
                if populationID not in onlyTEGroups[orthogroupID]:
                    onlyTEGroups[orthogroupID][populationID] = []
                for averageEntropy in entropyDictWithTEs[orthogroupID][populationID]:
                    onlyTEGroups[orthogroupID][populationID].append(averageEntropy)
    return(onlyTEGroups)


def plot(entropyDict,pop1,pop2,entropyRatioDict,proteinCountMinimum,onlyTEGroups):
    dataFrame = {}
    newDictLabel = pop1 + '_vs_' + pop2
    xList = []
    yList = []
    # baseColor = '#a1c9f4'
    baseColor = '#1f77b4'
    teColor = 'red'
    colorList = []
    legendDict = {}
    for orthogroupID in entropyDict:
        if pop1 in entropyDict[orthogroupID] and pop2 in entropyDict[orthogroupID]:
            # print(orthogroupID,pop1,pop2,len(entropyDict[orthogroupID][pop1]),len(entropyDict[orthogroupID][pop2]))
            if orthogroupID not in entropyRatioDict:
                entropyRatioDict[orthogroupID] = {}
            if newDictLabel not in entropyRatioDict[orthogroupID]:
                entropyRatioDict[orthogroupID][newDictLabel] = []
            for averageEntropy1 in entropyDict[orthogroupID][pop1]:
                xList.append(averageEntropy1)
                entropyRatioDict[orthogroupID][newDictLabel].append(averageEntropy1)
            for averageEntropy2 in entropyDict[orthogroupID][pop2]:
                yList.append(averageEntropy2)
                entropyRatioDict[orthogroupID][newDictLabel].append(averageEntropy2)
            colorList.append(baseColor)
            '''
            if orthogroupID not in dataFrame:
                dataFrame[orthogroupID] = {}
            if pop1 not in dataFrame[orthogroupID]:
                dataFrame[orthogroupID][pop1] = averageEntropy1
            if pop2 not in dataFrame[orthogroupID]:
                dataFrame[orthogroupID][pop2] = averageEntropy2
            '''
    for orthogroupID in onlyTEGroups:
        if pop1 in onlyTEGroups[orthogroupID] and pop2 in onlyTEGroups[orthogroupID]:
            for averageEntropy1 in onlyTEGroups[orthogroupID][pop1]:
                xList.append(averageEntropy1)
            for averageEntropy2 in onlyTEGroups[orthogroupID][pop2]:
                yList.append(averageEntropy2)
            colorList.append(teColor)
    # https://stackoverflow.com/questions/18837262/convert-python-dict-into-a-dataframe
    # https://stackoverflow.com/questions/27005783/changing-color-and-marker-of-each-point-using-seaborn-jointplot
    statistic,pvalue = pearsonr(xList, yList)
    print(pop1,pop2,statistic,pvalue)
    sns_plt = sns.jointplot(x=xList, y=yList, alpha=0.1)
    sns_plt.ax_joint.cla()
    plt.sca(sns_plt.ax_joint)
    plt.scatter(xList, yList, c=colorList, alpha=0.1, edgecolors='none')
    coefficients = np.polyfit(xList, yList, 1)
    m, b = np.poly1d(coefficients)
    xArray = np.array(xList)
    # plt.plot(xList, p(xList), color="black", linestyle='--')
    # plt.plot(x, m * x + b, color='red')
    plt.plot(xList, m * xArray + b, color="black", linestyle='--')

    lineLabel = mlines.Line2D([], [], color='black', marker='s',linestyle="None",markersize=10, label='Line of best fit\n' + f"y = {m:.2f}x + {b:.2f}", fillstyle='full', markeredgecolor='black', markeredgewidth=0.0)
    teLegendInfo = mlines.Line2D([], [], color='red', marker='o',linestyle="None",markersize=10, label='TE-related OGs', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    noTELegendInfo = mlines.Line2D([], [], color='#1f77b4', marker='o',linestyle="None",markersize=10, label='All other OGs', fillstyle='full', markeredgecolor='#1f77b4', markeredgewidth=0.0)
    if 'withTEs' not in legendDict:
        legendDict['withTEs'] = teLegendInfo
    if 'noTEs' not in legendDict:
        legendDict['noTEs'] = noTELegendInfo
    if 'BestFitLine' not in legendDict:
        legendDict['BestFitLine'] = lineLabel
        
    legendList = list(legendDict.values())
    plt.legend(handles=legendList, frameon=False, ncol=2)
    
    ### explicitly setting both axes to entropy 1.8
    #plt.xlim(0,1.8)
    #plt.ylim(0,1.8)
    maxBuffer = 0.1
    maxXValue = max(xList)
    maxYValue = max(yList)
    if maxXValue > maxYValue:
        maxValue = maxXValue
    elif maxXValue < maxYValue:
        maxValue = maxYValue
    else:
        maxValue = maxXValue
    newMax = maxValue + (maxValue * maxBuffer)
    # sns.jointplot(x=xList, y=yList, alpha=0.1)
    plt.xlim(0,newMax)
    plt.ylim(0,newMax)
    plt.xlabel(pop1 + ' entropy values')
    plt.ylabel(pop2 + ' entropy values')
    
    plt.title('Gene orthogroups (OGs)')
    plt.tight_layout()
    plt.savefig(pop1 + '_' + pop2 + '_jointHist_minCount' + str(proteinCountMinimum) + '.png',dpi=600)
    plt.savefig(pop1 + '_' + pop2 + '_jointHist_minCount' + str(proteinCountMinimum) + '.svg')
    plt.close()
    return(entropyRatioDict)


############
# MAIN  ####
############


usage = "Usage: " + sys.argv[0] + " <entropy output file with TEs> <entropoy file no TEs> <minimum protein count for each population> \n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

entropyFileWithTEs = sys.argv[1]
entropyFileNoTEs = sys.argv[2]
proteinCountMinimum = sys.argv[3]

proteinCountMinimum = int(proteinCountMinimum)


populationDict = {'asian_hemp':'1', 'F1':'1',
                  'feral':'1', 'hc_hemp':'1',
                  'hemp':'1', 'mj':'1'}
'''
populationDict = {'F1':'1',
                  'hc_hemp':'1',
                  'hemp':'1',
                  'mj':'1'}
'''
entropyDictWithTEs,orthogroupDictWithTEs,functionDictWithTEs,initDictWithTEs = readEntropyFile(entropyFileWithTEs,proteinCountMinimum)
entropyDictNoTEs,orthogroupDictNoTEs,functionDictNoTEs,initDictNoTEs = readEntropyFile(entropyFileNoTEs,proteinCountMinimum)

onlyTEGroups = collectRepeatGroups(entropyDictWithTEs,entropyDictNoTEs)

entropyRatioDict = {}
print("pop1\tpop2\tstatistic\tpvalue")
for i in populationDict:
    for j in populationDict:
        if i < j:
            entropyRatioDict = plot(entropyDictNoTEs,i,j,entropyRatioDict,proteinCountMinimum,onlyTEGroups)


entropyRatioList = []
OUT = open('entropy_ratios_minCount' + str(proteinCountMinimum) + '.txt','w')
OUT.write("orthogroupID\tlabelID\tpop1\togProteinCount1\taverageEntropy1\tpop2\togProteinCount2\taverageEntropy2\tentropyDifference\tentropyRatio\tannotation\n")
for orthogroupID in entropyRatioDict:
    if orthogroupID in orthogroupDictNoTEs:
        ogProteinCount = orthogroupDictNoTEs[orthogroupID]
    else:
        print(orthogroupID,'not in orthogroupDictNoTEs')
        sys.exit()
    if orthogroupID in functionDictNoTEs:
        annotation = functionDictNoTEs[orthogroupID]
    else:
        print(orthogroupID,'not in functionDictNoTEs')
        sys.exit()
    for labelID in entropyRatioDict[orthogroupID]:
        pop1,pop2 = labelID.split('_vs_')
        if len(entropyRatioDict[orthogroupID][labelID]) > 2:
            print("problem with number of values",orthogroupID,labelID,len(entropyRatioDict[orthogroupID][labelID]),entropyRatioDict[orthogroupID][labelID])
            sys.exit()
        else:
            averageEntropy1 = entropyRatioDict[orthogroupID][labelID][0]
            averageEntropy2 = entropyRatioDict[orthogroupID][labelID][1]
            averageEntropyDifference = averageEntropy1 - averageEntropy2
            
            if orthogroupID in initDictNoTEs and pop1 in initDictNoTEs[orthogroupID] and averageEntropy1 in initDictNoTEs[orthogroupID][pop1]:
                ogProteinCount1 = initDictNoTEs[orthogroupID][pop1][averageEntropy1]
            else:
                print(orthogroupID,pop1,averageEntropy1,'not in initDictNoTEs')
                sys.exit()

            if orthogroupID in	initDictNoTEs and pop2 in initDictNoTEs[orthogroupID] and averageEntropy2 in initDictNoTEs[orthogroupID][pop2]:
                ogProteinCount2 = initDictNoTEs[orthogroupID][pop2][averageEntropy2]
            else:
                print(orthogroupID,pop2,averageEntropy2,'not in initDictNoTEs')
                sys.exit()
                
            if averageEntropy1 > averageEntropy2:
                if averageEntropy2 > 0:
                    entropyRatio = float(averageEntropy1)/averageEntropy2
                else:
                    entropyRatio = 'nan'
            elif averageEntropy2 > averageEntropy1:
                if averageEntropy1 > 0:
                    entropyRatio = float(averageEntropy2)/averageEntropy1
                else:
                    entropyRatio = 'nan'
            else:
                # print(averageEntropy1,averageEntropy2)
                if averageEntropy1 == averageEntropy2 == 0:
                    entropyRatio = 'nan'
                else:
                    entropyRatio = float(averageEntropy2)/averageEntropy1
            #print(orthogroupID,labelID,averageEntropy1,averageEntropy2,entropyRatio)
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (orthogroupID,labelID,pop1,ogProteinCount1,averageEntropy1,pop2,ogProteinCount2,averageEntropy2,averageEntropyDifference,entropyRatio,annotation))

