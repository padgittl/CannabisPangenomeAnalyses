import sys,re,os, math, statistics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

###############
# SUBROUTINES #
###############

def collectHaplotypeGenes(geneData,mbh_genes_only,mcscanx_genes_only):
    initialGenePairs = {}
    preference = {}
    singleIDs = {}
    OUT = open('combinedGenePairs.txt','w')
    print("combinedGenePairs.txt")
    for geneID1 in geneData:
        if geneID1 in mcscanx_genes_only:
            geneID2 = mcscanx_genes_only[geneID1]
            # print(geneID1,geneID2)
            IDs = sorted([geneID1,geneID2])
            if (IDs[0],IDs[1]) not in initialGenePairs:
                initialGenePairs[(IDs[0],IDs[1])] = 'collinear'
            if geneID1 not in preference:
                preference[geneID1] = []
            preference[geneID1].append(geneID2)
            if geneID2 not in preference:
                preference[geneID2] = []
            preference[geneID2].append(geneID1)
            # print(geneID1,geneID2)
        else:
            if geneID1 in mbh_genes_only:
                geneID2 = mbh_genes_only[geneID1]
                IDs = sorted([geneID1,geneID2])
                if geneID1 not in preference and geneID2 not in preference:
                    # if geneID1 not in preference:
                    # print(geneID1,geneID2)
                    if (IDs[0],IDs[1]) not in initialGenePairs:
                        initialGenePairs[(IDs[0],IDs[1])] = 'mbh'

    otherHits = {}
    otherInfo = {}
    for geneID1,geneID2 in initialGenePairs:
        methodID = initialGenePairs[(geneID1,geneID2)]
        if geneID2 in otherHits:
            if otherHits[geneID2] == 'collinear':
                continue
        else:
            if methodID == 'collinear':
                otherHits[geneID2] = 'collinear'
                otherInfo[geneID2] = (geneID1,'collinear')

    genePairs = {}
    for geneID1,geneID2 in initialGenePairs:
        methodID = initialGenePairs[(geneID1,geneID2)]
        if geneID2 in otherInfo:
            otherGeneID,otherMethodID = otherInfo[geneID2]
            if otherGeneID == geneID1:
                # print(geneID1,geneID2,methodID)
                OUT.write("%s\t%s\t%s\n" % (geneID1,geneID2,methodID))
                if (geneID1,geneID2) not in genePairs:
                    genePairs[(geneID1,geneID2)] = methodID
        else:
            # print(geneID1,geneID2,methodID)
            OUT.write("%s\t%s\t%s\n" % (geneID1,geneID2,methodID))
            if (geneID1,geneID2) not in genePairs:
                genePairs[(geneID1,geneID2)] = methodID
        
    for geneID1 in preference:
        if len(preference[geneID1]) > 1:
            for geneID2 in preference[geneID1]:
                if (geneID1,geneID2) in initialGenePairs:
                    methodID = initialGenePairs[(geneID1,geneID2)]
                    #if geneID2 not in otherInfo:
                    # this is going to capture examples that are in both mbh and collinear
                    #print("multiples",geneID1,len(preference[geneID1]),preference[geneID1],methodID)
                    #sys.exit()
    return(genePairs)


def assessExpression(genePairs,avgTPMDict,hitDict,summaryAnnotDict,tissueSpecificDict):
    differenceThreshold = 5
    pairs = {}
    OUT = open('average_haplotype_biased_expression.txt','w')
    OUT_IDs_hapA = open('average_haplotype_biased_expression_geneIDs_' + assemblyID1 + '.txt','w')
    OUT_IDs_hapB = open('average_haplotype_biased_expression_geneIDs_' + assemblyID1 + '.txt','w')
    print('average_haplotype_biased_expression.txt')
    print('average_haplotype_biased_expression_geneIDs_' + assemblyID1 + '.txt')
    print('average_haplotype_biased_expression_geneIDs_' + assemblyID2 + '.txt')
    for geneID1,geneID2 in genePairs:
        pairID = geneID1 + "_" + geneID2
        if geneID1 in avgTPMDict and geneID2 in avgTPMDict:
            avgTPM1 = avgTPMDict[geneID1]
            avgTPM2 = avgTPMDict[geneID2]
            diff = abs(avgTPM1-avgTPM2)
            if geneID1 in hitDict and geneID2 in hitDict:
                annotation1 = hitDict[geneID1]
                annotation2 = hitDict[geneID2]
            else:
                if geneID1 in summaryAnnotDict and geneID2 in summaryAnnotDict:
                    annotation1 = summaryAnnotDict[geneID1]
                    annotation2 = summaryAnnotDict[geneID2]
                else:
                    print(geneID1,geneID2,'not in summaryAnnotDict')
                    sys.exit()
            if pairID not in pairs:
                pairs[pairID] = (avgTPM1,avgTPM2,diff,annotation1,annotation2)
            if diff >= differenceThreshold:
                expressionID = 'BIASED'
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('average',expressionID,diff,geneID1,avgTPM1,geneID2,avgTPM2,annotation1,annotation2))
                new_geneID1 = '"' + geneID1 + '"'
                new_geneID2 = '"' + geneID2 + '"'
                if avgTPM1 > avgTPM2:
                    OUT_IDs_hapA.write("%s\n" % (new_geneID1))
                elif avgTPM2 > avgTPM1:
                    OUT_IDs_hapB.write("%s\n" % (new_geneID2))
                else:
                    print("something else going on with avg TPM",geneID1,geneID2)
                    sys.exit()
            if diff < differenceThreshold:
                expressionID = 'BALANCED'
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('average',expressionID,diff,geneID1,avgTPM1,geneID2,avgTPM2,annotation1,annotation2))
                    
    for tissueID in tissueSpecificDict:
        TISSUE_OUT = open(tissueID + '_haplotype_biased_expression.txt','w')
        TISSUE_OUT_hapA = open(tissueID + '_haplotype_biased_expression_' + assemblyID1 + '.txt','w')
        TISSUE_OUT_hapB = open(tissueID + '_haplotype_biased_expression_' + assemblyID2 + '.txt','w')
        print(tissueID + '_haplotype_biased_expression.txt')
        print(tissueID + '_haplotype_biased_expression_' + assemblyID1 + '.txt')
        print(tissueID + '_haplotype_biased_expression_' + assemblyID2 + '.txt')
        for geneID1,geneID2 in genePairs:
            #print(geneID1,geneID2)
            pairID = geneID1 + "_" + geneID2
            avgTPM1,avgTPM2,avgDiff,annotation1,annotation2 = pairs[pairID]
            if geneID1 in tissueSpecificDict[tissueID] and geneID2 in tissueSpecificDict[tissueID]:
                if geneID1 in goDict and geneID2 in goDict:
                    goTerms1 = goDict[geneID1]
                    goTerms2 = goDict[geneID2]
                else:
                    print(geneID1,geneID2,'no go terms')
                    sys.exit()
                tpmValue1 = tissueSpecificDict[tissueID][geneID1]
                tpmValue2 = tissueSpecificDict[tissueID][geneID2]
                diff = abs(tpmValue1-tpmValue2)
                if diff >= differenceThreshold:
                    expressionID = 'BIASED'
                    TISSUE_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tissueID,expressionID,diff,geneID1,tpmValue1,avgTPM1,geneID2,tpmValue2,avgTPM2,avgDiff,annotation1,annotation2,goTerms1,goTerms2))
                    new_geneID1 = '"' + geneID1 + '"'
                    new_geneID2 = '"' + geneID2 + '"'
                    if tpmValue1 > tpmValue2:
                        TISSUE_OUT_hapA.write("%s\n" % (new_geneID1))
                    elif tpmValue2 > tpmValue1:
                        TISSUE_OUT_hapB.write("%s\n" % (new_geneID2))
                    else:
                        print("something else going on with avg TPM",geneID1,geneID2)
                        sys.exit()
                if diff < differenceThreshold:
                    expressionID = 'BALANCED'
                    TISSUE_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tissueID,expressionID,diff,geneID1,tpmValue1,avgTPM1,geneID2,tpmValue2,avgTPM2,avgDiff,annotation1,annotation2,goTerms1,goTerms2))
                    

# is primary
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
                    # print(transcriptID,geneID)
                    if transcriptID not in geneData:
                        geneData[transcriptID] = (int(start),int(end))
    return(geneData)


def read_blast_output_file(blast_output_file):
    bestHits = {}
    bestScore = {}
    eValueThreshold = 0.001 # 1e-3
    with open(blast_output_file,'r') as F:
        for line in F:
            queryID,subjectID,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = line.strip().split('\t')
            qItems = queryID.split('.')
            qGenomeID = qItems[0]
            qChrID = qItems[1]
            sItems = subjectID.split('.')
            sGenomeID = sItems[0]
            sChrID = sItems[1]
            # print(queryID,subjectID)
            if queryID != subjectID:
                if qGenomeID != sGenomeID:
                    if qChrID == sChrID:
                        # print(queryID,subjectID)
                        pident = float(pident)
                        length = int(length)
                        evalue = float(evalue)
                        bitscore = float(bitscore)
                        if evalue < eValueThreshold:
                            if queryID in bestHits:
                                if bitscore > bestScore[queryID]:
                                    bestScore[queryID] = bitscore
                                    bestHits[queryID] = (subjectID,pident,evalue,bitscore)
                            else:
                                bestScore[queryID] = bitscore
                                bestHits[queryID] = (subjectID,pident,evalue,bitscore)
    return(bestHits)


def printMutualBestHits(bestHits1,bestHits2,assemblyID1,assemblyID2):
    mbh_blast_results = {}
    mbh_genes_only = {}
    OUT = open(assemblyID1 + '_vs_' + assemblyID2 + '_mutual_best_hits.txt','w')
    for queryID1 in bestHits1:
        subjectID1,pIdent1,eValue1,bitScore1 = bestHits1[queryID1]
        if subjectID1 in bestHits2:
            subjectID2,pIdent2,eValue2,bitScore2 = bestHits2[subjectID1]
            IDs = sorted([queryID1,subjectID1])
            labelID = IDs[0] + "_vs_" + IDs[1]
            maxPIden = max(pIdent1,pIdent2)
            minEValue = min(eValue1,eValue2)
            maxBitScore = max(bitScore1,bitScore2)
            if subjectID2 == queryID1:
                #if queryID1 != subjectID1:
                OUT.write("%s\t%s\n" % (IDs[0],IDs[1]))
                if labelID not in mbh_blast_results:
                    mbh_blast_results[labelID] = (IDs[0],IDs[1],maxPIden,minEValue,maxBitScore)
                if IDs[0] not in mbh_genes_only:
                    mbh_genes_only[IDs[0]] = IDs[1]
                if IDs[1] not in mbh_genes_only:
                    mbh_genes_only[IDs[1]] = IDs[0]
    return(mbh_blast_results,mbh_genes_only)


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


def parseCollinearityFile(collinearityFile):
    genome1BlockDict = {}
    genome2BlockDict = {}
    genome1_vs_genome2BlockDict = {}
    mcscanx_genes_only = {}
    checkDuplicates = {}
    blockCount = {}
    geneToBlock = {}
    geneCount = {}

    with open(collinearityFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                fullNumberBlock,geneID1,geneID2,eValue = line.strip().split('\t')
                cols1 = geneID1.split('.')
                assemblyID1 = cols1[0]
                chrID1 = cols1[1]
                cols2 = geneID2.split('.')
                assemblyID2 = cols2[0]
                chrID2 = cols2[1]
                numberBlock,extra = fullNumberBlock.split('-')
                numberBlock = int(numberBlock)
                IDs = sorted([geneID1,geneID2])
                labelID = IDs[0] + "_vs_" + IDs[1]
                if assemblyID1 != assemblyID2:
                    if chrID1 == chrID2:
                        if numberBlock not in blockCount:
                            blockCount[numberBlock] = 0
                        blockCount[numberBlock] += 1
                        if geneID1 not in geneToBlock:
                            geneToBlock[geneID1] = {}
                        if numberBlock not in geneToBlock[geneID1]:
                            geneToBlock[geneID1][numberBlock] = 1
                        if geneID2 not in geneToBlock:
                            geneToBlock[geneID2] = {}
                        if numberBlock not in geneToBlock[geneID2]:
                            geneToBlock[geneID2][numberBlock] = 1

    bestHits1 = {}
    topBlock1 = {}
    bestHits2 = {}
    topBlock2 = {}
    with open(collinearityFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                fullNumberBlock,geneID1,geneID2,eValue = line.strip().split('\t')
                cols1 = geneID1.split('.')
                assemblyID1 = cols1[0]
                chrID1 = cols1[1]
                cols2 = geneID2.split('.')
                assemblyID2 = cols2[0]
                chrID2 = cols2[1]
                numberBlock,extra = fullNumberBlock.split('-')
                numberBlock = int(numberBlock)
                IDs = sorted([geneID1,geneID2])
                labelID = IDs[0] + "_vs_" + IDs[1]
                if assemblyID1 != assemblyID2:
                    if chrID1 == chrID2:
                        # blockCount[numberBlock]
                        blockCountID = blockCount[numberBlock]
                        if geneID1 in bestHits1:
                            if blockCountID > bestHits1[geneID1]:
                                bestHits1[geneID1] = blockCountID
                                topBlock1[geneID1] = (geneID2,numberBlock,blockCountID)
                        else:
                            bestHits1[geneID1] = blockCountID
                            topBlock1[geneID1] = (geneID2,numberBlock,blockCountID)

                        if geneID2 in bestHits2:
                            if blockCountID > bestHits2[geneID2]:
                                bestHits2[geneID2] = blockCountID
                                topBlock2[geneID2] = (geneID1,numberBlock,blockCountID)
                        else:
                            bestHits2[geneID2] = blockCountID
                            topBlock2[geneID2] = (geneID1,numberBlock,blockCountID)

    final_mcscanx_genes_only = {}
    for queryID1 in topBlock1:
        subjectID1,numberBlock1,blockCountID1 = topBlock1[queryID1]
        if subjectID1 in topBlock2:
            subjectID2,numberBlock2,blockCountID2 = topBlock2[subjectID1]
            IDs = sorted([queryID1,subjectID1])
            if subjectID2 == queryID1:
                if IDs[0] not in final_mcscanx_genes_only:
                    final_mcscanx_genes_only[IDs[0]] = IDs[1]
                if IDs[1] not in final_mcscanx_genes_only:
                    final_mcscanx_genes_only[IDs[1]] = IDs[0]
    return(final_mcscanx_genes_only)


def readExpressionFile(expressionFile,geneData,assemblyID1,assemblyID2,tpmThreshold):
    sampleIDList = []
    tpmCountDict = {}
    geneIDList = []
    with open(expressionFile,'r') as F:
        for line in F:
            if '\"\"' in line:
                line = line.strip().split(',')
                for sampleID in line[1:]:
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
                    # https://stackoverflow.com/questions/39988048/what-is-the-origin-of-matplotlibs-symlog-a-k-a-symmetrical-log-scale
                    # symLogValue = symmetric_logarithm(tpmCount)
                    # tpmCountDict[geneID].append(symLogValue)
                    tpmCountDict[geneID].append(tpmCount)
    # geneData[transcriptID] = geneID
    filteredData = {}
    filteredSampleIDs = {}
    tissueSpecificDict = {}
    for geneID in geneData:
        start,stop = geneData[geneID]
        if geneID in tpmCountDict:
            for i in range(len(tpmCountDict[geneID])):
                tpmValue = tpmCountDict[geneID][i]
                sampleID = sampleIDList[i]
                # print(sampleID)
                if geneID not in filteredData:
                    filteredData[geneID] = []
                filteredData[geneID].append(tpmValue)
                if geneID not in filteredSampleIDs:
                    filteredSampleIDs[geneID] = []
                filteredSampleIDs[geneID].append(sampleID)
                if sampleID not in tissueSpecificDict:
                    tissueSpecificDict[sampleID] = {}
                if geneID not in tissueSpecificDict[sampleID]:
                    tissueSpecificDict[sampleID][geneID] = tpmValue
    avgTPMDict = {}
    for geneID in filteredData:
        avgTPM = float(sum(filteredData[geneID])) / len(filteredData[geneID])
        # print(geneID,avgTPM)
        if geneID not in avgTPMDict:
            avgTPMDict[geneID] = avgTPM
    return(tpmCountDict,sampleIDList,filteredData,filteredSampleIDs,avgTPMDict,tissueSpecificDict)


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
    fullSummaryDict = {}
    goDict = {}
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
                if transcriptID not in summaryAnnotDict:
                    summaryAnnotDict[transcriptID] = annotation
    df = pd.read_csv(fileID, sep="\t", low_memory=False, on_bad_lines='warn')
    # for (colname,colval) in df.iteritems():
    for (colname,colval) in df.items():
        if colname not in fullSummaryDict:
            fullSummaryDict[colname] = []
        fullSummaryDict[colname].append(colval.values)
    combinedList = []
    for i in range(len(fullSummaryDict['TID'])):
        for j in range(len(fullSummaryDict['TID'][i])):
            combinedList.append((fullSummaryDict['TID'][i][j],str(fullSummaryDict['GOS'][i][j])))
            # print(fullSummaryDict['TID'][i][j],fullSummaryDict['GOS'][i][j])
            # print(fullSummaryDict)
            if fullSummaryDict['TID'][i][j] not in goDict:
                goDict[fullSummaryDict['TID'][i][j]] = fullSummaryDict['GOS'][i][j]
    return(goDict,summaryAnnotDict)

            
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <expression file> <hit file> <summary tsv file> <gene gff file> <haplotype gene file> <collinearity file> <blast file 1> <blast file 2> <genome ID 1> <genome ID 2> <TPM threshold, e.g. 0.1> \n"
if len(sys.argv) != 12:
    print(usage)
    sys.exit()

expressionFile = sys.argv[1]
hitFile = sys.argv[2]
summaryFile = sys.argv[3]
gffFile = sys.argv[4]
haplotypeGeneFile = sys.argv[5]
collinearityFile = sys.argv[6]
blast_output_file1 = sys.argv[7]
blast_output_file2 = sys.argv[8]
assemblyID1 = sys.argv[9]
assemblyID2 = sys.argv[10]
tpmThreshold = sys.argv[11]

final_mcscanx_genes_only = parseCollinearityFile(collinearityFile)
bestHits1 = read_blast_output_file(blast_output_file1)
bestHits2 = read_blast_output_file(blast_output_file2)
mbh_blast_results,mbh_genes_only = printMutualBestHits(bestHits1,bestHits2,assemblyID1,assemblyID2)

geneData = readGFF(gffFile)
goDict,summaryAnnotDict = readSummaryFile(summaryFile)
hitDict = readHitFile(hitFile)

tpmCountDict,sampleIDList,filteredData,filteredSampleIDs,avgTPMDict,tissueSpecificDict = readExpressionFile(expressionFile,geneData,assemblyID1,assemblyID2,tpmThreshold)

genePairs = collectHaplotypeGenes(geneData,mbh_genes_only,final_mcscanx_genes_only)
assessExpression(genePairs,avgTPMDict,hitDict,summaryAnnotDict,tissueSpecificDict)

