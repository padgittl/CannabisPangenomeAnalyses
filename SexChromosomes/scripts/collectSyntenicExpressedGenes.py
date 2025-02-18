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

def readSyntenyFile(syntenyFile):
    prepSyntenicGenes = {}
    syntenicGenes = {}
    with open(syntenyFile,'r') as F:
        for line in F:
            if 'AH3Ma' in line or 'AH3Mb' in line:
                line = line.strip().split('\t')
                if 'NotSyntenic' not in line:
                    genomeID,syntenyID,start,end,geneID,ogID,femaleAvgTPM,maleAvgTPM,absDiff,uniprotAnnotation,annotation,goTerms = line
                    if geneID not in prepSyntenicGenes:
                        prepSyntenicGenes[geneID] = []
                    prepSyntenicGenes[geneID].append(syntenyID)
                else:
                    geneID = line[4]
                    # print(geneID)
                    syntenyID = 'NotSyntenic'
                    if geneID not in prepSyntenicGenes:
                        prepSyntenicGenes[geneID] = []
                    prepSyntenicGenes[geneID].append(syntenyID)
    for geneID in prepSyntenicGenes:
        if len(prepSyntenicGenes[geneID]) > 1:
            print(geneID,prepSyntenicGenes[geneID],'present in > 1 OG')
            sys.exit()
        else:
            for syntenyID in prepSyntenicGenes[geneID]:
                if geneID not in syntenicGenes:
                    syntenicGenes[geneID] = syntenyID
    return(syntenicGenes)


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


def read_blast_output_file(blast_output_file):
    # AH3Ma.chr1.v1.g082810.t2        AH3Ma.chr1.v1.g082810.t2        100.0   326     0       0       1       326     1       326     9.96e-238       646.0
    bestHits = {}
    bestScore = {}
    eValueThreshold = 0.001 # 1e-3
    with open(blast_output_file,'r') as F:
        for line in F:
            queryID,subjectID,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = line.strip().split('\t')
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
            # print(queryID1,subjectID1)
            IDs = sorted([queryID1,subjectID1])
            labelID = IDs[0] + "_vs_" + IDs[1]
            maxPIden = max(pIdent1,pIdent2)
            minEValue = min(eValue1,eValue2)
            maxBitScore = max(bitScore1,bitScore2)
            if subjectID2 == queryID1:
                OUT.write("%s\t%s\n" % (IDs[0],IDs[1]))
                if labelID not in mbh_blast_results:
                    mbh_blast_results[labelID] = (IDs[0],IDs[1],maxPIden,minEValue,maxBitScore)
                if IDs[0] not in mbh_genes_only:
                    mbh_genes_only[IDs[0]] = IDs[1]
                if IDs[1] not in mbh_genes_only:
                    mbh_genes_only[IDs[1]] = IDs[0]
    return(mbh_blast_results,mbh_genes_only)


def readExpressionFile(expressionFile,geneData,assemblyID1,assemblyID2,tpmThreshold,specificTissueID):
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
                if specificTissueID in sampleID:
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


def assessDiffExp(geneData,filteredData,filteredSampleIDs,avgTPMDict,syntenicGenes,specificTissueID,absTPMThreshold):
    tpmDifferenceThreshold = 5
    diffExpGenes = {}
    for geneID in geneData:
        if 'chr' in geneID:
            if geneID in filteredData:
                # print(geneID)
                femaleAvgTPM = avgTPMDict[geneID]['female']
                maleAvgTPM = avgTPMDict[geneID]['male']
                absDiff = abs(femaleAvgTPM-maleAvgTPM)
                if absDiff >= tpmDifferenceThreshold:
                    # print(geneID,femaleAvgTPM,maleAvgTPM,absDiff)
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
                        # this is checking for cases where female expression is zero and male expression > 0
                        if maleAvgTPM > absTPMThreshold:
                            if geneID not in lowExpDiff:
                                lowExpDiff[geneID] = {}
                            if 'male' not in lowExpDiff[geneID]:
                                lowExpDiff[geneID]['male'] = maleAvgTPM
                            if 'female' not in lowExpDiff[geneID]:
                                lowExpDiff[geneID]['female'] = femaleAvgTPM
                if len(maleTPMSet) == 1:
                    if 0.0 not in maleTPMSet:
                        print(geneID,maleTPMSet)
                        sys.exit()
                    else:
                        if femaleAvgTPM > absTPMThreshold:
                            if geneID not in lowExpDiff:
                                lowExpDiff[geneID] = {}
                            if 'male' not in lowExpDiff[geneID]:
                                lowExpDiff[geneID]['male'] = maleAvgTPM
                            if 'female'	not in lowExpDiff[geneID]:
                                lowExpDiff[geneID]['female'] = femaleAvgTPM

    OUT1_MALE = open(specificTissueID + '_male_only_expression.txt','w')
    OUT1_MALE.write("labelID\tgeneID\tstart\tstop\tsyntenyID\tfemaleID\tfemaleAvgTPM\tmaleID\tmaleAvgTPM\ttpmDiff\tannotation\n")
    print(specificTissueID + "_male_only_expression.txt")
    OUT1_FEMALE = open(specificTissueID + '_female_only_expression.txt','w')
    OUT1_FEMALE.write("labelID\tgeneID\tstart\tstop\tsyntenyID\tfemaleID\tfemaleAvgTPM\tmaleID\tmaleAvgTPM\ttpmDiff\tannotation\n")
    print(specificTissueID + "_female_only_expression.txt")
    for geneID in lowExpDiff:
        if geneID in syntenicGenes:
            syntenyID = syntenicGenes[geneID]
        else:
            print(geneID,'not in syntenicGenes')
            sys.exit()
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
        # print(geneID)
        if maleAvgTPM > femaleAvgTPM:
            OUT1_MALE.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("lowExpDiff",geneID,start,stop,syntenyID,'female',femaleAvgTPM,'male',maleAvgTPM,diff,annotation))
        if femaleAvgTPM > maleAvgTPM:
            OUT1_FEMALE.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("lowExpDiff",geneID,start,stop,syntenyID,'female',femaleAvgTPM,'male',maleAvgTPM,diff,annotation))
            
    OUT2_MALE = open(specificTissueID + '_male_higher_exp_min5TPM_diffs.txt','w')
    OUT2_MALE.write("labelID\tgeneID\tstart\tstop\tsyntenyID\tfemaleID\tfemaleAvgTPM\tmaleID\tmaleAvgTPM\ttpmDiff\tannotation\n")
    print(specificTissueID + '_male_higher_exp_min5TPM_diffs.txt')

    OUT2_FEMALE = open(specificTissueID + '_female_higher_exp_min5TPM_diffs.txt','w')
    OUT2_FEMALE.write("labelID\tgeneID\tstart\tstop\tsyntenyID\tfemaleID\tfemaleAvgTPM\tmaleID\tmaleAvgTPM\ttpmDiff\tannotation\n")
    print(specificTissueID + '_female_higher_exp_min5TPM_diffs.txt')
    for geneID in diffExpGenes:
        if geneID not in lowExpDiff:
            if geneID in syntenicGenes:
                syntenyID = syntenicGenes[geneID]
            else:
                print(geneID,'not in syntenicGenes')
                sys.exit()
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
            if maleAvgTPM > femaleAvgTPM:
                # this if-statement is to collect examples where > 5 TPM but either male/female also has avg 0.0 expression
                if femaleAvgTPM == 0.0:
                    OUT1_MALE.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("lowExpDiff",geneID,start,stop,syntenyID,'female',femaleAvgTPM,'male',maleAvgTPM,diff,annotation))
                else:
                    OUT2_MALE.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('Min5tpmdiff',geneID,start,stop,syntenyID,'female',femaleAvgTPM,'male',maleAvgTPM,diff,annotation))
            if femaleAvgTPM > maleAvgTPM:
                if maleAvgTPM == 0.0:
                    OUT1_FEMALE.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("lowExpDiff",geneID,start,stop,syntenyID,'female',femaleAvgTPM,'male',maleAvgTPM,diff,annotation))
                else:
                    OUT2_FEMALE.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('min5TPMDiff',geneID,start,stop,syntenyID,'female',femaleAvgTPM,'male',maleAvgTPM,diff,annotation))

    OUT_BALANCED = open(specificTissueID + '_balanced_male_female_expression.txt','w')
    print(specificTissueID + '_balanced_male_female_expression.txt')
    for geneID in geneData:
        if 'chr' in geneID:
            if geneID in filteredData:
                if geneID not in lowExpDiff:
                    if geneID in syntenicGenes:
                        syntenyID = syntenicGenes[geneID]
                    else:
                        print(geneID,'not in syntenicGenes')
                        sys.exit()
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
                    if geneID in avgTPMDict:
                        # print(geneID)
                        femaleAvgTPM = avgTPMDict[geneID]['female']
                        maleAvgTPM = avgTPMDict[geneID]['male']
                        absDiff = abs(femaleAvgTPM-maleAvgTPM)
                        print(geneID,femaleAvgTPM,maleAvgTPM,absDiff)
                        # absTPMThreshold
                        # this is requiring some minimal level of expression, otherwise expression for many genes will be 0
                        if femaleAvgTPM >= absTPMThreshold and maleAvgTPM >= absTPMThreshold:
                            if absDiff < tpmDifferenceThreshold:
                                OUT_BALANCED.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("balanced",geneID,start,stop,syntenyID,'female',femaleAvgTPM,'male',maleAvgTPM,absDiff,annotation))
                    else:
                        print(geneID,'not in avgTPMDict')
                        sys.exit()

            
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


def parseCollinearityFile(collinearityFile):
    genome1BlockDict = {}
    genome2BlockDict = {}
    genome1_vs_genome2BlockDict = {}
    mcscanx_genes_only = {}
    with open(collinearityFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                fullNumberBlock,geneID1,geneID2,eValue = line.strip().split('\t')
                cols1 = geneID1.split('.')
                assemblyID1 = cols1[0]
                cols2 = geneID2.split('.')
                assemblyID2 = cols2[0]
                # print(assemblyID1,assemblyID2)
                numberBlock,extra = fullNumberBlock.split('-')
                numberBlock = int(numberBlock)
                IDs = sorted([geneID1,geneID2])
                labelID = IDs[0] + "_vs_" + IDs[1]
                if assemblyID1 != assemblyID2:
                    if labelID not in genome1_vs_genome2BlockDict:
                        genome1_vs_genome2BlockDict[labelID] = 1
                    if IDs[0] not in mcscanx_genes_only:
                        mcscanx_genes_only[IDs[0]] = IDs[1]
                    if IDs[1] not in mcscanx_genes_only:
                        mcscanx_genes_only[IDs[1]] = IDs[0]
    return(genome1BlockDict,genome2BlockDict,genome1_vs_genome2BlockDict,mcscanx_genes_only)


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

usage = "Usage: " + sys.argv[0] + " <collinearity file> <expression file> <blast output file 1> <blast output file 2> <hit file> <summary tsv file> <gene gff file> <syntenic genes> <genome ID 1> <genome ID 2> <TPM threshold, e.g. 0.0> <minimum TPM threshold to be considered expressed>\n"
if len(sys.argv) != 13:
    print(usage)
    sys.exit()

collinearityFile = sys.argv[1]
expressionFile = sys.argv[2]
blast_output_file1 = sys.argv[3]
blast_output_file2 = sys.argv[4]
hitFile = sys.argv[5]
summaryFile = sys.argv[6]
gffFile = sys.argv[7]
syntenyFile = sys.argv[8]
assemblyID1 = sys.argv[9]
assemblyID2 = sys.argv[10]
tpmThreshold = sys.argv[11]
minAbsTPM = sys.argv[12]

minAbsTPM = float(minAbsTPM)

geneData = readGFF(gffFile)
summaryAnnotDict = readSummaryFile(summaryFile)
hitDict = readHitFile(hitFile)

syntenicGenes =	readSyntenyFile(syntenyFile)
flower_tpmCountDict,flower_sampleIDList,flower_filteredData,flower_filteredSampleIDs,flower_avgTPMDict = readExpressionFile(expressionFile,geneData,assemblyID1,assemblyID2,tpmThreshold,'flower')
assessDiffExp(geneData,flower_filteredData,flower_filteredSampleIDs,flower_avgTPMDict,syntenicGenes,'flower',minAbsTPM)
leaf_tpmCountDict,leaf_sampleIDList,leaf_filteredData,leaf_filteredSampleIDs,leaf_avgTPMDict = readExpressionFile(expressionFile,geneData,assemblyID1,assemblyID2,tpmThreshold,'leaf')
assessDiffExp(geneData,leaf_filteredData,leaf_filteredSampleIDs,leaf_avgTPMDict,syntenicGenes,'leaf',minAbsTPM)


