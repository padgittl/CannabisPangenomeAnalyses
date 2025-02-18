import sys, re, os
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np

###############
# SUBROUTINES #
###############

def readMaleFemaleAvgExpressionFile(maleFemaleAvgExpFile):
    expressionData = {}
    with open(maleFemaleAvgExpFile,'r') as F:
        for line in F:
            if 'femaleAvgTPM' not in line:
                geneID,femaleAvgTPM,maleAvgTPM,absDiff = line.strip().split('\t')
                if geneID not in expressionData:
                    expressionData[geneID] = (float(femaleAvgTPM),float(maleAvgTPM),float(absDiff))
    return(expressionData)


def readUniprotHitFile(uniprotHitFile):
    uniprotHits = {}
    with open(uniprotHitFile,'r') as F:
        for line in F:
            if 'uniprotGeneID' not in line:
                geneID,mRNAStart,mRNAStop,uniprotGeneID,uniprotDescription,percentIdentity,eValue,bitScore,queryCoverage = line.strip().split('\t')
                if geneID not in uniprotHits:
                    uniprotHits[geneID] = (uniprotGeneID,uniprotDescription)
    return(uniprotHits)


def readTerpeneInfoFile(terpeneInfoFile):
    terpeneInfoDict = {}
    with open(terpeneInfoFile,'r') as F:
        for line in F:
            if 'chemotypeID' not in line:
                assemblyID,chemotypeID,chrID,transcriptID,mRNAStart,mRNAStop,cdsLength,exonCount,strand,fractionMasked,percentIdentity,eValue,queryCoverage,uniprotAnnotation,pfamAnnotation = line.strip().split('\t')
                if transcriptID not in terpeneInfoDict:
                    terpeneInfoDict[transcriptID] = (uniprotAnnotation,pfamAnnotation)
    return(terpeneInfoDict)


def readGenomeList(genomeList):
    genomeDict = {}
    with open(genomeList,'r') as F:
        for line in F:
            genomeID = line.strip()
            if genomeID not in genomeDict:
                genomeDict[genomeID] = 1
    genomePrefixes = {}
    for genomeID in genomeDict:
        genomePrefix,extra = genomeID.split(genomeID[-1])
        if genomePrefix not in genomePrefixes:
            genomePrefixes[genomePrefix] = 1
    genomeCount = len(genomePrefixes)
    # print(genomeCount)
    return(genomeDict,genomeCount)


def readGFFWithOGsFile(gffWithOGsFile,goDict,summaryAnnotDict,genomeDict,genomeCount,terpeneInfoDict,uniprotHits):
    ogDict = {}
    geneCoords = {}
    verifyAssignment = {}
    with open(gffWithOGsFile,'r') as F:
        for line in F:
            if 'isArrayRep' not in line:
                genomeID,geneID,chromID,start,end,strand,order,ofID,pepLen,globOG,arrayID,isArrayRep,arrayOrd,synOG,inblkOG,og = line.strip().split('\t')
                og = 'ogID_' + og
                if geneID not in geneCoords:
                    geneCoords[geneID] = (int(start),int(end))
                if og not in ogDict:
                    ogDict[og] = []
                ogDict[og].append(geneID)
                if geneID not in verifyAssignment:
                    verifyAssignment[geneID] = []
                verifyAssignment[geneID].append(og)
    for geneID in verifyAssignment:
        if len(verifyAssignment[geneID]) > 1:
            print("PROBLEM: geneID is in multiple ogs (should only be one) -->",geneID,verifyAssignment[geneID])
            sys.exit()
    ogInfo = {}
    chromRepresentation = {}
    genomeRepresentation = {}
    for ogID in ogDict:
        # AH3Ma.chrX.v1.g406630.t2
        if len(ogDict[ogID]) > 1:
            for geneID in ogDict[ogID]:
                items = geneID.split('.')
                genomeID = items[0]
                chrID = items[1]
                genomePrefix,extra = genomeID.split(genomeID[-1])
                if ogID not in chromRepresentation:
                    chromRepresentation[ogID] = []
                chromRepresentation[ogID].append(chrID)
    for ogID in ogDict:
        if len(ogDict[ogID]) > 1:
            if ogID in chromRepresentation:
                # this is requiring that only one chromosome is represented
                if len(set(chromRepresentation[ogID])) == 1:
                    for geneID in ogDict[ogID]:
                        items = geneID.split('.')
                        genomeID = items[0]
                        chrID = items[1]
                        genomePrefix,extra = genomeID.split(genomeID[-1])
                        # print(ogID,genomePrefix,chrID,geneID)
                        if ogID not in genomeRepresentation:
                            genomeRepresentation[ogID] = {}
                        if genomePrefix not in genomeRepresentation[ogID]:
                            genomeRepresentation[ogID][genomePrefix] = 1
                        if ogID not in ogInfo:
                            ogInfo[ogID] = {}
                        if chrID not in ogInfo[ogID]:
                            ogInfo[ogID][chrID] = {}
                        if genomeID not in ogInfo[ogID][chrID]:
                            ogInfo[ogID][chrID][genomeID] = []
                        ogInfo[ogID][chrID][genomeID].append(geneID)
                else:
                    # print(ogID,ogDict[ogID])
                    if len(set(chromRepresentation[ogID])) == 2:
                        if 'chrX' in set(chromRepresentation[ogID]) and 'chrY' in set(chromRepresentation[ogID]):
                            # print(ogID,chromRepresentation[ogID],ogDict[ogID])
                            for geneID in ogDict[ogID]:
                                items = geneID.split('.')
                                genomeID = items[0]
                                chrID = items[1]
                                genomePrefix,extra = genomeID.split(genomeID[-1])
                                if ogID not in genomeRepresentation:
                                    genomeRepresentation[ogID] = {}
                                if genomePrefix not in genomeRepresentation[ogID]:
                                    genomeRepresentation[ogID][genomePrefix] = 1
                                if ogID not in ogInfo:
                                    ogInfo[ogID] = {}
                                if chrID not in ogInfo[ogID]:
                                    ogInfo[ogID][chrID] = {}
                                if genomeID not in ogInfo[ogID][chrID]:
                                    ogInfo[ogID][chrID][genomeID] = []
                                ogInfo[ogID][chrID][genomeID].append(geneID)
                                
    # goDict[fullSummaryDict['TID'][i][j]] = fullSummaryDict['GOS'][i][j]
    # summaryAnnotDict[transcriptID] = annotation
    OUT_BOTH = open('present_in_chrX_and_chrY.txt','w')
    OUT_BOTH.write("### ogID\tchrID\tgenomeID\tgeneID\tstart\tend\tannotation\tgoTerms\n")
    OUT_X = open('present_in_only_chrX.txt','w')
    OUT_X.write("### ogID\tchrID\tgenomeID\tgeneID\tstart\tend\tannotation\tgoTerms\n")
    OUT_Y = open('present_in_only_chrY.txt','w')
    OUT_Y.write("### ogID\tchrID\tgenomeID\tgeneID\tstart\tend\tannotation\tgoTerms\n")
    OUT_AUTO = open('present_in_only_autosomes.txt','w')
    OUT_AUTO.write("### ogID\tchrID\tgenomeID\tgeneID\tstart\tend\tannotation\tgoTerms\n")
    OUT_TERPS = open('terpene_synthase_related_syntenic_genes.txt','w')
    OUT_TERPS.write("### ogID\tchrID\tgenomeID\tgeneID\tstart\tend\tuniprot_annotation\tpfam_annotation\tsummary_annotation\tgoTerms\n")
    print("present_in_chrX_and_chrY.txt")
    print("present_in_only_chrX.txt")
    print("present_in_only_chrY.txt")
    print("present_in_only_autosomes.txt")
    print("terpene_synthase_related_syntenic_genes.txt")
    allGenomes = sorted(set({'AH3Ma':'1', 'AH3Mb':'1', 'BCMa':'1', 'BCMb':'1', 'GRMa':'1', 'GRMb':'1', 'KOMPa':'1', 'KOMPb':'1'}))
    femaleGenomes = sorted(set({'AH3Ma':'1', 'BCMa':'1', 'GRMb':'1', 'KOMPa':'1'}))
    maleGenomes = sorted(set({'AH3Mb':'1', 'BCMb':'1', 'GRMa':'1', 'KOMPb':'1'}))
    syntenic_groupings = {}
    prepGenomeSetDict = {}
    genomeSetDict = {}
    # ogInfo[ogID][chrID][genomeID].append(geneID)
    for ogID in ogInfo:
        for chrID in ogInfo[ogID]:
            for genomeID in ogInfo[ogID][chrID]:
                if ogID not in prepGenomeSetDict:
                    prepGenomeSetDict[ogID] = []
                prepGenomeSetDict[ogID].append(genomeID)
    for ogID in prepGenomeSetDict:
        if ogID not in genomeSetDict:
            genomeSetDict[ogID] = sorted(set(prepGenomeSetDict[ogID]))
        
    for ogID in ogInfo:
        chromSet = set(ogInfo[ogID])
        if ogID in genomeRepresentation:
            genomeSet = set(genomeRepresentation[ogID])
            numberGenomes = len(genomeSet)
        else:
            print(ogID,'not in genomeRepresentation')
            sys.exit()
        for chrID in ogInfo[ogID]:
            for genomeID in ogInfo[ogID][chrID]:
                # GRMb
                # KOMPa
                for geneID in ogInfo[ogID][chrID][genomeID]:
                    if geneID in uniprotHits:
                        uniprotGeneID,uniprotDescription = uniprotHits[geneID]
                        uniprotAnnotation = uniprotGeneID + ' ' + uniprotDescription
                    else:
                        uniprotAnnotation = 'NoAnnotation'
                    if geneID in summaryAnnotDict:
                        transcript_start,transcript_stop,annotation = summaryAnnotDict[geneID]
                    else:
                        print(geneID,'not in summaryAnnotDict')
                        sys.exit()
                    if geneID in goDict:
                        goTerms = goDict[geneID]
                    else:
                        goTerms = 'NoGOTerms'
                    if geneID in geneCoords:
                        # geneCoords[geneID] = (int(start),int(end))
                        start,end = geneCoords[geneID]
                    else:
                        print(geneID,'not in geneCoords')
                        sys.exit()

                    if ogID in genomeSetDict:
                        ogGenomes = genomeSetDict[ogID]
                        # print(ogID,ogGenomes,chromSet)
                        if 'chrY' in chromSet and 'chrX' in chromSet:
                            if ogID in genomeSetDict:
                                ogGenomes = genomeSetDict[ogID]
                                # print(ogID,ogGenomes,chromSet)
                                if ogGenomes == allGenomes:
                                    # print(ogGenomes,allGenomes)
                                    OUT_BOTH.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ogID,chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                                    if geneID not in syntenic_groupings:
                                        syntenic_groupings[geneID] = []
                                    syntenic_groupings[geneID].append(('bothXY',ogID,chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                            else:
                                print(ogID,'not in genomeSetDict')
                                sys.exit()
                                # these are going to have occupancy in all eight genomes
                                # print(ogID,chrID,genomeID,geneID,start,end,annotation,goTerms)
                                # ogID_549 chrX AH3Ma AH3Ma.chrX.v1.g431490.t2 51572183 51575745 protein modification by small protein removal NoGOTerms
                        if 'chrY' not in chromSet and 'chrX' in chromSet:
                            # femaleGenomes
                            if ogID in genomeSetDict:
                                ogGenomes = genomeSetDict[ogID]
                                # print(ogID,ogGenomes,femaleGenomes)
                                if ogGenomes == femaleGenomes:
                                    OUT_X.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ogID,chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                                    if geneID not in syntenic_groupings:
                                        syntenic_groupings[geneID] = []
                                    syntenic_groupings[geneID].append(('onlyChrX',ogID,chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                            else:
                                print(ogID,'not in genomeSetDict')
                                sys.exit()
                        if 'chrY' in chromSet and 'chrX' not in chromSet:
                            if ogID in genomeSetDict:
                                ogGenomes = genomeSetDict[ogID]
                                if ogGenomes == maleGenomes:
                                    # print(ogID,ogGenomes,maleGenomes)
                                    OUT_Y.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ogID,chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                                    if geneID not in syntenic_groupings:
                                        syntenic_groupings[geneID] = []
                                    syntenic_groupings[geneID].append(('onlyChrY',ogID,chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                            else:
                                print(ogID,'not in genomeSetDict')
                                sys.exit()
                        # autosomes
                        if 'chrX' not in chromSet and 'chrY' not in chromSet:
                            if ogID in genomeSetDict:
                                ogGenomes = genomeSetDict[ogID]
                                if ogGenomes == femaleGenomes:
                                    OUT_AUTO.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ogID,"FemaleHaplotypes",chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                                    if geneID not in syntenic_groupings:
                                        syntenic_groupings[geneID] = []
                                    syntenic_groupings[geneID].append(('syntenicAutosomeFemaleHaplotypes',ogID,chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                                elif ogGenomes == maleGenomes: 
                                    OUT_AUTO.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ogID,"MaleHaplotypes",chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                                    if geneID not in syntenic_groupings:
                                        syntenic_groupings[geneID] = []
                                    syntenic_groupings[geneID].append(('syntenicAutosomeMaleHaplotypes',ogID,chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                                elif ogGenomes == allGenomes:
                                    OUT_AUTO.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ogID,"AllHaplotypes",chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                                    if geneID not in syntenic_groupings:
                                        syntenic_groupings[geneID] = []
                                    syntenic_groupings[geneID].append(('syntenicAutosomeAllHaplotypes',ogID,chrID,genomeID,geneID,start,end,uniprotAnnotation,annotation,goTerms))
                            else:
                                print(ogID,'not in genomeSetDict')
                                sys.exit()
                        if geneID in terpeneInfoDict:
                            terp_uniprotAnnotation,terp_pfamAnnotation = terpeneInfoDict[geneID]
                            OUT_TERPS.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ogID,chrID,genomeID,geneID,start,end,terp_uniprotAnnotation,terp_pfamAnnotation,annotation,goTerms))
    return(syntenic_groupings)


def collectDataForBoxplot(summaryAnnotDict,syntenic_groupings,expressionData):
    prepData = {}
    orderedGenesForBoxplot = {}
    for geneID in summaryAnnotDict:
        transcriptStart,transcriptStop,annotation = summaryAnnotDict[geneID]
        items = geneID.split('.')
        assemblyID = items[0]
        chromID = items[1]
        if 'chr' in geneID:
            if geneID in syntenic_groupings:
                if len(syntenic_groupings[geneID]) > 1:
                    print(geneID,syntenic_groupings[geneID],'more than one entry; investigate')
                    sys.exit()
                else:
                    for groupingID,ogID,chrID,genomeID,gene_ID,start,end,uniprotAnnotation,annotation,goTerms in syntenic_groupings[geneID]:
                        if assemblyID not in prepData:
                            prepData[assemblyID] = {}
                        if chromID not in prepData[assemblyID]:
                            prepData[assemblyID][chromID] = {}
                        if ogID not in prepData[assemblyID][chromID]:
                            prepData[assemblyID][chromID][ogID] = (groupingID,start,end,ogID,chrID,genomeID,geneID,uniprotAnnotation,annotation,goTerms)
                            
    for assemblyID in prepData:
        for chromID in prepData[assemblyID]:
            for ogID in prepData[assemblyID][chromID]:
                groupingID,start,end,ogID,chrID,genomeID,geneID,uniprotAnnotation,annotation,goTerms = prepData[assemblyID][chromID][ogID]
                if assemblyID not in orderedGenesForBoxplot:
                    orderedGenesForBoxplot[assemblyID] = {}
                if chromID not in orderedGenesForBoxplot[assemblyID]:
                    orderedGenesForBoxplot[assemblyID][chromID] = []
                orderedGenesForBoxplot[assemblyID][chromID].append((groupingID,start,end,ogID,chrID,genomeID,geneID,uniprotAnnotation,annotation,goTerms))
                    
    for geneID in summaryAnnotDict:
        transcriptStart,transcriptStop,annotation = summaryAnnotDict[geneID]
        items = geneID.split('.')
        assemblyID = items[0]
        chromID = items[1]
        if 'chr' in geneID:
            if geneID not in syntenic_groupings:
                groupingID = 'NotSyntenic'
                if assemblyID not in orderedGenesForBoxplot:
                    orderedGenesForBoxplot[assemblyID] = {}
                if chromID not in orderedGenesForBoxplot[assemblyID]:
                    orderedGenesForBoxplot[assemblyID][chromID] = []
                orderedGenesForBoxplot[assemblyID][chromID].append((groupingID,transcriptStart,transcriptStop,geneID,annotation))
    # just get one representative per syntenic orthogroup -- some of them are not actually syntenic
    nonSyntenicGenesForPlot = {}
    syntenicGenesForPlot = {}
    for assemblyID in orderedGenesForBoxplot:
        for chromID in orderedGenesForBoxplot[assemblyID]:
            for item in orderedGenesForBoxplot[assemblyID][chromID]:
                if 'NotSyntenic' in item:
                    groupingID,start,end,geneID,annotation = item
                    if chromID not in nonSyntenicGenesForPlot:
                        nonSyntenicGenesForPlot[chromID] = {}
                    if assemblyID not in nonSyntenicGenesForPlot[chromID]:
                        nonSyntenicGenesForPlot[chromID][assemblyID] = []
                    nonSyntenicGenesForPlot[chromID][assemblyID].append(geneID)
                else:
                    groupingID,start,end,ogID,chrID,genomeID,geneID,uniprotAnnotation,annotation,goTerms = item
                    if chromID not in syntenicGenesForPlot:
                        syntenicGenesForPlot[chromID] = {}
                    if assemblyID not in syntenicGenesForPlot[chromID]:
                        syntenicGenesForPlot[chromID][assemblyID] = []
                    syntenicGenesForPlot[chromID][assemblyID].append(geneID)
    return(orderedGenesForBoxplot,syntenicGenesForPlot,nonSyntenicGenesForPlot)
                    

# summaryAnnotDict[transcriptID] = annotation
def assessSynteny(summaryAnnotDict,syntenic_groupings,expressionData):
    orderedGenes = {}
    for geneID in summaryAnnotDict:
        transcriptStart,transcriptStop,annotation = summaryAnnotDict[geneID]
        items = geneID.split('.')
        assemblyID = items[0]
        chromID = items[1]
        if 'chr' in geneID:
            if geneID in syntenic_groupings:
                if len(syntenic_groupings[geneID]) > 1:
                    print(geneID,syntenic_groupings[geneID],'more than one entry; investigate')
                    sys.exit()
                else:
                    for groupingID,ogID,chrID,genomeID,gene_ID,start,end,uniprotAnnotation,annotation,goTerms in syntenic_groupings[geneID]:
                        if assemblyID not in orderedGenes:
                            orderedGenes[assemblyID] = {}
                        if chromID not in orderedGenes[assemblyID]:
                            orderedGenes[assemblyID][chromID] = []
                        orderedGenes[assemblyID][chromID].append((groupingID,start,end,ogID,chrID,genomeID,geneID,uniprotAnnotation,annotation,goTerms))
            else:
                groupingID = 'NotSyntenic'
                if assemblyID not in orderedGenes:
                    orderedGenes[assemblyID] = {}
                if chromID not in orderedGenes[assemblyID]:
                    orderedGenes[assemblyID][chromID] =	[]
                orderedGenes[assemblyID][chromID].append((groupingID,transcriptStart,transcriptStop,geneID,annotation))
    for assemblyID in orderedGenes:
        for chrID in orderedGenes[assemblyID]:
            orderedGenes[assemblyID][chrID].sort(key=lambda x: x[1], reverse=False)
    OUT = open('ordered_syntenic_and_non_syntenic_genes.txt','w')
    print("ordered_syntenic_and_non_syntenic_genes.txt")
    # expressionData[geneID] = (float(femaleAvgTPM),float(maleAvgTPM),float(absDiff))
    nonSyntenicGenes = {}
    syntenicGenes = {}
    for assemblyID in orderedGenes:
        if assemblyID != 'AH3Ma' and assemblyID != 'AH3Mb':
            for chrID in orderedGenes[assemblyID]:
                for item in orderedGenes[assemblyID][chrID]:
                    if 'NotSyntenic' in item:
                        groupingID,start,end,geneID,annotation = item
                        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,groupingID,start,end,geneID,annotation))
                        if chrID not in nonSyntenicGenes:
                            nonSyntenicGenes[chrID] = {}
                        if assemblyID not in nonSyntenicGenes[chrID]:
                            nonSyntenicGenes[chrID][assemblyID] = []
                        nonSyntenicGenes[chrID][assemblyID].append(geneID)
                    else:
                        groupingID,start,end,ogID,chrID,genomeID,geneID,uniprotAnnotation,annotation,goTerms = item
                        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,groupingID,start,end,geneID,ogID,uniprotAnnotation,annotation,goTerms))
                        if chrID not in	syntenicGenes:
                            syntenicGenes[chrID] = {}
                        if assemblyID not in syntenicGenes[chrID]:
                            syntenicGenes[chrID][assemblyID] = []
                        syntenicGenes[chrID][assemblyID].append(geneID)
        else:
            for chrID in orderedGenes[assemblyID]:
                for item in orderedGenes[assemblyID][chrID]:
                    if 'NotSyntenic' in item:
                        groupingID,start,end,geneID,annotation = item
                        if geneID in expressionData:
                            femaleAvgTPM,maleAvgTPM,absDiff = expressionData[geneID]
                            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,groupingID,start,end,geneID,femaleAvgTPM,maleAvgTPM,absDiff,annotation))
                        else:
                            print(geneID,'not in expressionData')
                            sys.exit()
                        if chrID not in nonSyntenicGenes:
                            nonSyntenicGenes[chrID] = {}
                        if assemblyID not in nonSyntenicGenes[chrID]:
                            nonSyntenicGenes[chrID][assemblyID]	= []
                        nonSyntenicGenes[chrID][assemblyID].append(geneID)
                    else:
                        groupingID,start,end,ogID,chrID,genomeID,geneID,uniprotAnnotation,annotation,goTerms = item
                        if geneID in expressionData:
                            femaleAvgTPM,maleAvgTPM,absDiff = expressionData[geneID]
                            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,groupingID,start,end,geneID,ogID,femaleAvgTPM,maleAvgTPM,absDiff,uniprotAnnotation,annotation,goTerms))
                        else:
                            print(geneID,'not in expressionData')
                            sys.exit()
                        if chrID not in syntenicGenes:
                            syntenicGenes[chrID] = {}
                        if assemblyID not in syntenicGenes[chrID]:
                            syntenicGenes[chrID][assemblyID] = []
                        syntenicGenes[chrID][assemblyID].append(geneID)
    TOTAL_COUNTS = open('genome_total_synteny_counts.txt','w')
    print("genome_total_synteny_counts.txt")
    for chrID in syntenicGenes:
        for assemblyID in syntenicGenes[chrID]:
            synGeneCount = len(syntenicGenes[chrID][assemblyID])
            # print("Syntenic:",chrID,assemblyID,synGeneCount)
            TOTAL_COUNTS.write("%s\t%s\t%s\t%s\n" % ("SYNTENIC",chrID,assemblyID,synGeneCount))
    for chrID in nonSyntenicGenes:
        for assemblyID in nonSyntenicGenes[chrID]:
            nonSynGeneCount = len(nonSyntenicGenes[chrID][assemblyID])
            # print("NonSyntenic:",chrID,assemblyID,nonSynGeneCount)
            TOTAL_COUNTS.write("%s\t%s\t%s\t%s\n" % ("NON_SYNTENIC",chrID,assemblyID,nonSynGeneCount))
    return(syntenicGenes,nonSyntenicGenes)


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
        # print(exonID,dataStdDev)
        fullDataArray.append(np.asarray(dataList))
        fullStdDevList.append(dataStdDev)
    # https://stackoverflow.com/questions/19931975/sort-multiple-lists-simultaneously
    idList,fullDataArray,fullStdDevList = map(list, zip(*sorted(zip(idList, fullDataArray, fullStdDevList), reverse=False)))
    return(fullDataArray,fullStdDevList,idList)


def set_axis_style(ax, labels):
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    plt.xticks(rotation=45, ha="right")
    ax.tick_params(axis='x', labelsize=6) 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def createBoxplot(dataArray,stdDevList,idList,keywordID,labelID):
    colorDict = {'syntenic':'blue', 'nonSyntenic':'red'}
    # Create a figure instance for subs
    plt.rcParams["figure.figsize"] = [4,2]
    colors = ['#278B9A']*len(stdDevList)
    
    # Create the boxplot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    labelList = []
    exonCount = 0
    for infoID in idList:
        label = infoID
        labelList.append(label)

    bp = ax.boxplot(dataArray, labels = labelList, patch_artist=True, showfliers=False, zorder=0, showmeans=True, meanline=True)
    ax.tick_params(axis='x', labelsize=6)
    plt.xticks(rotation=45, ha="right")

    plt.ylabel('Gene counts',size=10)

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_alpha(0.6)

    for i in range(len(dataArray)):
        y1 = dataArray[i]
        x1 = np.random.normal(i+1, 0.02, len(y1))
        plt.plot(x1, y1, 'black', alpha=0.3, marker="o", markersize=3, linestyle=' ',
                 fillstyle='full',
                 markeredgecolor='black',
                 markeredgewidth=0.0)

    plt.title(labelID + ' genes')
    plt.savefig(keywordID + '_boxplot.png', bbox_inches='tight', dpi=600)
    plt.savefig(keywordID + '_boxplot.svg')
    print(keywordID + '_boxplot.png')
    plt.close()


def prepBoxPlot(syntenicGenes,nonSyntenicGenes):
    fullSyntenicBoxData = {}
    fullNonSyntenicBoxData = {}
    fullBoxData = {}
    boxData = {}
    for chrID in nonSyntenicGenes:
        for assemblyID in nonSyntenicGenes[chrID]:
            nonSynGeneCount = len(nonSyntenicGenes[chrID][assemblyID])
            if chrID not in boxData:
                boxData[chrID] = {}
            if 'nonSyntenic' not in boxData[chrID]:
                boxData[chrID]['nonSyntenic'] = []
            boxData[chrID]['nonSyntenic'].append(nonSynGeneCount)
            if chrID not in fullNonSyntenicBoxData:
                fullNonSyntenicBoxData[chrID] = []
            fullNonSyntenicBoxData[chrID].append(nonSynGeneCount)
    for chrID in syntenicGenes:
        for assemblyID in syntenicGenes[chrID]:
            synGeneCount = len(syntenicGenes[chrID][assemblyID])
            if chrID not in boxData:
                boxData[chrID] = {}
            if 'syntenic' not in boxData[chrID]:
                boxData[chrID]['syntenic'] = []
            boxData[chrID]['syntenic'].append(synGeneCount)
            if chrID not in fullSyntenicBoxData:
                fullSyntenicBoxData[chrID] = []
            fullSyntenicBoxData[chrID].append(synGeneCount)
    '''
    finalBoxData = {}
    labelsList = []
    for chrID in boxData:
        nonSyntenicCounts = boxData[chrID]['nonSyntenic']
        syntenicCounts = boxData[chrID]['syntenic']
        labelsList.append(chrID)
        if chrID not in finalBoxData:
            finalBoxData[chrID] = [syntenicCounts, nonSyntenicCounts]
    '''
    return(boxData,fullSyntenicBoxData,fullNonSyntenicBoxData)

    
def readSummaryFile(summaryFile):
    summaryAnnotDict = {}
    fullSummaryDict = {}
    goDict = {}
    fileID = summaryFile.strip()
    maskedThreshold = 1.0
    with open(summaryFile,'r') as F:
        for line in F:
            if 'PERCENT_MASKED' not in line:
                line = line.strip().split('\t')
                geneID = line[0]
                transcriptID = line[1]
                annotation = line[-1]
                fractionMasked = line[8]
                fractionMasked = float(fractionMasked)
                # print(transcriptID,fractionMasked)
                start = line[3]
                stop = line[4]
                start = int(start)
                stop = int(stop)
                # GID     TID     CTG     START   STOP
                isPrimary = line[9]
                if isPrimary == 'True':
                    if fractionMasked < maskedThreshold:
                        if transcriptID not in summaryAnnotDict:
                            summaryAnnotDict[transcriptID] = (start,stop,annotation)
    df = pd.read_csv(fileID, sep="\t", low_memory=False, on_bad_lines='warn')
    for (colname,colval) in df.items():
        if colname not in fullSummaryDict:
            fullSummaryDict[colname] = []
        fullSummaryDict[colname].append(colval.values)
    combinedList = []
    for i in range(len(fullSummaryDict['TID'])):
        for j in range(len(fullSummaryDict['TID'][i])):
            combinedList.append((fullSummaryDict['TID'][i][j],str(fullSummaryDict['GOS'][i][j])))
            if fullSummaryDict['TID'][i][j] not in goDict:
                goDict[fullSummaryDict['TID'][i][j]] = fullSummaryDict['GOS'][i][j]
    return(goDict,summaryAnnotDict)
        

############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <gffWithOgs.txt> <concatenated summary file> <genome list file> <terpene annotation file> <uniprot hit file> <male-female average expression file>\n"
if len(sys.argv) != 7:
    print(usage)
    sys.exit()

gffWithOGsFile = sys.argv[1]
summaryFile = sys.argv[2]
genomeList = sys.argv[3]
terpeneInfoFile = sys.argv[4]
uniprotHitFile = sys.argv[5]
maleFemaleAvgExpFile = sys.argv[6]

expressionData = readMaleFemaleAvgExpressionFile(maleFemaleAvgExpFile)
uniprotHits = readUniprotHitFile(uniprotHitFile)
terpeneInfoDict = readTerpeneInfoFile(terpeneInfoFile)
genomeDict,genomeCount = readGenomeList(genomeList)
goDict,summaryAnnotDict = readSummaryFile(summaryFile)
syntenic_groupings = readGFFWithOGsFile(gffWithOGsFile,goDict,summaryAnnotDict,genomeDict,genomeCount,terpeneInfoDict,uniprotHits)

syntenicGenes,nonSyntenicGenes = assessSynteny(summaryAnnotDict,syntenic_groupings,expressionData)

orderedGenesForBoxplot,syntenicGenesForPlot,nonSyntenicGenesForPlot = collectDataForBoxplot(summaryAnnotDict,syntenic_groupings,expressionData)
boxData,fullSyntenicBoxData,fullNonSyntenicBoxData = prepBoxPlot(syntenicGenesForPlot,nonSyntenicGenesForPlot)

syn_fullDataArray,syn_fullStdDevList,syn_idList = createDataArray(fullSyntenicBoxData)
non_syn_fullDataArray,non_syn_fullStdDevList,non_syn_idList = createDataArray(fullNonSyntenicBoxData)
createBoxplot(syn_fullDataArray,syn_fullStdDevList,syn_idList,'syntenic','Syntenic')
createBoxplot(non_syn_fullDataArray,non_syn_fullStdDevList,non_syn_idList,'nonSyntenic','Non-syntenic')
