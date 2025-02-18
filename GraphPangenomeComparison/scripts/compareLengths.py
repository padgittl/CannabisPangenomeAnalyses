#!/bin/python
import sys, os, re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
# from statistics import mean, median, stdev
import statistics
import math
import scipy.stats 

###############
# SUBROUTINES #
###############

def readSyriSVFile(syriFile,nonLogDataPoints,logDataPoints):
    syriData = {}
    logData = {}
    medianDict = {}
    nonLogDict = {}
    with open(syriFile,'r') as F:
        # print(F)
        for line in F:
            # print(line)
            if 'svLength' not in line:
                genomeID,chrID,start,stop,varLen,logLen = line.strip().split('\t')
                if 'chr' in chrID:
                    start = int(start)
                    stop = int(stop)
                    varLen = int(varLen)
                    logLen = float(logLen)
                    nonLogDataPoints.append(varLen)
                    logDataPoints.append(logLen)
                    # print(varLen,logLen)
                    if chrID not in syriData:
                        syriData[chrID] = []
                    syriData[chrID].append(int(varLen))
                    # syriData[chrID].append(logLen)
                    if chrID not in nonLogDict:
                        nonLogDict[chrID] = []
                    nonLogDict[chrID].append(varLen)
                    if chrID not in logData:
                        logData[chrID] = []
                    logData[chrID].append(logLen)
    nonLogDict = dict(sorted(nonLogDict.items()))
    for chrID in nonLogDict:
        medianValue = statistics.median(nonLogDict[chrID])
        print("Syri",chrID,round(medianValue,3))
        if chrID not in medianDict:
            medianDict[chrID] = medianValue
    return(syriData,nonLogDataPoints,logDataPoints,medianDict,logData)


# https://seqan.readthedocs.io/en/seqan-v2.0.2/Tutorial/VcfIO.html#:~:text=In%20text%20representations%2C%20such%20as%20VCF%2C%201%2Dbased,data%20structures%2C%20SeqAn%20uses%200%2Dbased%20half%2Dopen%20intervals.
# https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
# In text representations, such as VCF, 1-BASED CLOSED intervals are used
# included means CLOSED; excluded means OPEN, this is from interval terminology
# Size = stop - start + 1
# https://github.com/samtools/bcftools/issues/661
# https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants/
# abs(strlen(ALT)-strlen(REF))<=50
# ##INFO=<ID=SVLEN,Number=-1,Type=Integer,Description="Difference in length between REF and ALT alleles">, take longer of two variants if comma separated, but first investigate

# ['AH3Ma.chrX', '10761', '>733598>733601', 'G', 'A', '60.0', '.', 'AC=1;AF=0.5;AN=2;AT=>733598>733599>733601,>733598>733600>733601;NS=2;LV=0', 'GT', '.', '.', '.', '.', '.', '.', '1', '.', '.', '.', '0', '.']
def readVCF(vcfFile,nonLogDataPoints,logDataPoints,graphDict,logData):
    with open(vcfFile,'r') as F:
        # assemblyID      chrID   pos     refAllele       altAllele       svLen   logSVLen
        for line in F:
            if not line.startswith('#') and 'refAllele' not in line:
                # line = line.strip().split('\t')
                # print(line)
                assemblyID,chrID,pos,refAllele,altAllele,svLen,logSVLen = line.strip().split('\t')
                fullChromID = assemblyID + '.' + chrID
                pos = int(pos)
                if 'chr' in fullChromID:
                    if chrID not in graphDict:
                        graphDict[chrID] = []
                    if chrID not in logData:
                        logData[chrID] = []
                    if ',' in refAllele:
                        print('unexpected comman in refAllele',refAllele)
                        sys.exit()
                    if ',' in altAllele:
                        altItems = altAllele.split(',')
                        # print(len(altItems),altItems)
                        itemList = []
                        for i in altItems:
                            i = i.strip()
                            # print(i,len(i))
                            itemList.append(len(i))
                            # print(assemblyID,chrID,pos,refAllele,altAllele)
                        # print(max(itemList),itemList)
                        maxValue = max(itemList)
                        diff = abs(len(refAllele)-maxValue)
                        # if diff > 0:
                        if diff >= 50:
                            logLen = math.log(diff)
                            logData[chrID].append(logLen)
                            nonLogDataPoints.append(diff)
                            logDataPoints.append(logLen)
                            graphDict[chrID].append(diff)
                    else:
                        # print(len(altAllele),altAllele)
                        # graphDict[assemblyID][chrID].append(len(altAllele))
                        diff = abs(len(refAllele)-len(altAllele))
                        #if diff > 0:
                        if diff >= 50:
                            logLen = math.log(diff)
                            logData[chrID].append(logLen)
                            nonLogDataPoints.append(diff)
                            logDataPoints.append(logLen)
                            graphDict[chrID].append(diff)
    return(graphDict,logData,nonLogDataPoints,logDataPoints)

#Chrom   Start   End     Name
#chr1    13839866        13852176        AH3Ma
#chr6    38673417        38685697        DUP     AH3Ma
#refGenomeID     chrID   queryStart      queryStop       varLen  logLen
#bed files are 0-based
#AH3Ma   chr6    38673417        38685697        12280   9.415727201701133
def readMGCSVFile(mgcFile,nonLogDataPoints,logDataPoints):
    mgcData = {}
    logData = {}
    with open(mgcFile,'r') as F:
        # print(F)
        for line in F:
            # print(line)
            if 'svLength' not in line and 'refGenomeID' not in line:
                genomeID,chrID,start,stop,varLen,logLen = line.strip().split('\t')
                start = int(start)
                stop = int(stop)
                varLen = int(varLen)
                logLen = float(logLen)
                nonLogDataPoints.append(varLen)
                logDataPoints.append(logLen)
                if chrID not in logData:
                    logData[chrID] = []
                logData[chrID].append(logLen)
                if chrID not in mgcData:
                    mgcData[chrID] = []
                mgcData[chrID].append(int(varLen))
    return(mgcData,logData,nonLogDataPoints,logDataPoints)


def compare(data1,data2,labelID):
    chrIDs = {'chr1':'1', 'chr2':'1', 'chr3':'1', 'chr4':'1', 'chr5':'1', 'chr6':'1', 'chr7':'1', 'chr8':'1', 'chr9':'1', 'chrX':'1'}
    for chrID in chrIDs:
        if chrID in data1 and chrID in data2:
            ks_statistic,pValue = scipy.stats.ks_2samp(data1[chrID], data2[chrID], alternative='two-sided', method='auto', axis=0, nan_policy='propagate', keepdims=False)
            #wilcoxon_result = scipy.stats.wilcoxon(syriData[chrID], graphData[chrID], alternative='two-sided')
            print(labelID,chrID,round(ks_statistic,2),pValue)
            #print(chrID,wilcoxon_result)
            

def getStdDev(dataList):
    dataStdDev = round(np.std(dataList), 2)
    return(dataStdDev)


def createDataArray(specificDataDict):
    fullDataArray = []
    fullStdDevList = []
    idList = []
    for infoID in specificDataDict:
        idList.append(infoID)
        #print(infoID)
        dataList = specificDataDict[infoID]
        #print(exonID,dataList,type(dataList))
        dataStdDev = getStdDev(dataList)
        # createIndividualBoxplot(exonID,dataList,dataStdDev)
        # print(exonID,dataStdDev)
        fullDataArray.append(np.asarray(dataList))
        fullStdDevList.append(dataStdDev)
    # https://stackoverflow.com/questions/19931975/sort-multiple-lists-simultaneously
    idList,fullDataArray,fullStdDevList = map(list, zip(*sorted(zip(idList, fullDataArray, fullStdDevList), reverse=False)))
    return(fullDataArray,fullStdDevList,idList)


# https://python-charts.com/distribution/violin-plot-matplotlib/
# https://matplotlib.org/stable/gallery/statistics/violinplot.html
# https://matplotlib.org/stable/gallery/statistics/customized_violin.html#sphx-glr-gallery-statistics-customized-violin-py
def createPlot(pggb_fullDataArray,pggb_fullStdDevList,pggb_fullIDList,
           mgc_fullDataArray,mgc_fullStdDevList,mgc_fullIDList,
           syri_fullDataArray,syri_fullStdDevList,syri_fullIDList,
	   nonLogMax,logMax):
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(10,8), sharex=True)
    ax1.set_ylabel('PGGB\nvariants', fontsize=10)
    ax2.set_ylabel('MGC\nvariants', fontsize=10)
    ax3.set_ylabel('All Syri\nstructural variants', fontsize=10)
    ax1.set_title('Gaussian kernel density of\nlog-transformed variant lengths')

    labelList = []
    for infoID in pggb_fullIDList:
        label = infoID
        labelList.append(label)

    parts1 = ax1.violinplot(pggb_fullDataArray,
                            showmeans=False, showmedians=False,
                            showextrema=False, widths=0.95)
    parts2 = ax2.violinplot(mgc_fullDataArray,
                            showmeans=False, showmedians=False,
                            showextrema=False, widths=0.95)
    parts3 = ax3.violinplot(syri_fullDataArray,
                            showmeans=False, showmedians=False,
                            showextrema=False, widths=0.95)
    ax3.set_xlabel('Chromosome IDs')
    # nonLogMax,logMax
    for ax in [ax1, ax2, ax3]:
        ax.set_ylim(0,logMax)
        # ax.set_ylim(0,nonLogMax)
    
    for ax in [ax1, ax2, ax3]:
        set_axis_style(ax, labelList)

    for pc in parts1["bodies"]:
        set_part_properties(pc)

    for pc in parts2["bodies"]:
        set_part_properties(pc)

    for pc in parts3["bodies"]:
        set_part_properties(pc)

    plt.savefig('graph_vs_syri_length_comparison.png', bbox_inches='tight', dpi=600)
    #plt.savefig(outName + '.svg', bbox_inches='tight')
    plt.close()


def set_axis_style(ax, labels):
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    plt.xticks(rotation=45, ha="right")
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def set_part_properties(pc):
    pc.set_facecolor("none")
    pc.set_edgecolor('black')
    pc.set_linewidth(1)
    pc.set_alpha(1)

    
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <pggb vcf file 1> <pggb vcf file 2> <syri file> <mgc file> \n"
if len(sys.argv) != 5:
    print(usage)
    sys.exit()

vcfFile1 = sys.argv[1]
vcfFile2 = sys.argv[2]
syriFile = sys.argv[3]
mgcFile = sys.argv[4]

nonLogDataPoints = []
logDataPoints = []

syriData,nonLogDataPoints,logDataPoints,medianDict,syriLogData = readSyriSVFile(syriFile,nonLogDataPoints,logDataPoints)
# print(nonLogDataPoints,logDataPoints)

pggbGraphData = {}
pggbLogData = {}

pggbGraphData,pggbLogData,nonLogDataPoints,logDataPoints = readVCF(vcfFile1,nonLogDataPoints,logDataPoints,pggbGraphData,pggbLogData)
pggbGraphData,pggbLogData,nonLogDataPoints,logDataPoints = readVCF(vcfFile2,nonLogDataPoints,logDataPoints,pggbGraphData,pggbLogData)

# print(nonLogDataPoints,logDataPoints)

mgcData,mgcLogData,nonLogDataPoints,logDataPoints = readMGCSVFile(mgcFile,nonLogDataPoints,logDataPoints)

# print(nonLogDataPoints,logDataPoints)
# syriData[chrID].append(logLen), logData[chrID].append(logLen)
# graphDict[chrID].append(diff), logData[chrID].append(logLen)
# mgcData[chrID].append(int(varLen)), logData[chrID].append(logLen)

compare(mgcData,pggbGraphData,'MGC vs PGGB')
compare(syriData,pggbGraphData,'Syri vs PGGB')
compare(syriData,mgcData,'Syri vs MGC')

'''
nonLogMaxValues = []
# print(nonLogDataPoints)
for i in nonLogDataPoints:
    maxValue = max(i)
    #print(maxValue)
    nonLogMaxValues.append(maxValue)
'''
nonLogMax = max(nonLogDataPoints)
nonLogMax = nonLogMax+(nonLogMax*0.1)

'''
logMaxValues = []
for i in logDataPoints:
    maxValue = max(i)
    # print(maxValue)
    logMaxValues.append(maxValue)
'''
logMax = max(logDataPoints)
logMax = logMax+(logMax*0.1)

pggb_fullDataArray,pggb_fullStdDevList,pggb_fullIDList = createDataArray(pggbLogData)
mgc_fullDataArray,mgc_fullStdDevList,mgc_fullIDList = createDataArray(mgcLogData)
syri_fullDataArray,syri_fullStdDevList,syri_fullIDList = createDataArray(syriLogData)

createPlot(pggb_fullDataArray,pggb_fullStdDevList,pggb_fullIDList,
           mgc_fullDataArray,mgc_fullStdDevList,mgc_fullIDList,
           syri_fullDataArray,syri_fullStdDevList,syri_fullIDList,
           nonLogMax,logMax)

