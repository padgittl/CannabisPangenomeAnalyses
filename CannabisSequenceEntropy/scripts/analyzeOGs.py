import sys, re, os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from scipy.stats import entropy
from collections import Counter


###############
# SUBROUTINES #
###############


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
                        summaryAnnotDict[transcriptID] = (annotation,fractionMasked)
    return(summaryAnnotDict)


def readFastaFile(fastaFileList,notCannabis):
    recordDict = {}
    for fastaFile in fastaFileList:
        fileBase = os.path.basename(fastaFile)
        orthogroupID,fileSuffix = os.path.splitext(fileBase)
        for record in SeqIO.parse(fastaFile,"fasta"):
            if '_' in record.id:
                items = record.id.split('_')
                if len(items) == 2:
                    assemblyID,geneID = items
                else:
                    assemblyID = items[0]
                    geneID = items[1] + '_' + items[2]
                    if assemblyID in notCannabis:
                        continue
                    else:
                        print("unexpected exception",assemblyID,geneID)
                        sys.exit()
            else:
                print("underscore not in",fastaFile,record.id)
                sys.exit()
            if orthogroupID not in recordDict:
                recordDict[orthogroupID] = {}
            if geneID not in recordDict[orthogroupID]:
                recordDict[orthogroupID][geneID] = (assemblyID,record)
    return(recordDict)


'''
FragariaVesca
LotusJaponicus
MalusDomestica
PrunusPersica
RosaChinensis

orthogroupDict[orthogroupID].append(geneID)
'''
def readOrthogroupTSV(orthogroupTSV,msaDataDict,populationDict,notCannabis,summaryAnnotDict,outputLabel):
    orthogroupDict = {}
    orthogroupProteinCount = {}
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
                            i = i.strip()
                            # only cannabis gene annotations will have '.'
                            if '.' in i:
                                cols = i.split('.')
                                assemblyID = cols[0]
                                if assemblyID != 'F49v1a':
                                    if assemblyID not in notCannabis:
                                        if orthogroupID not in orthogroupDict:
                                            orthogroupDict[orthogroupID] = []
                                        orthogroupDict[orthogroupID].append(i)
                    else:
                        # item = list(filter(str.strip, item))
                        geneID = item.strip()
                        # print(geneID)
                        if '.' in geneID:
                            items = geneID.split('.')
                            assemblyID = items[0]
                            if assemblyID != 'F49v1a':
                                if assemblyID not in notCannabis:
                                    if orthogroupID not in orthogroupDict:
                                        orthogroupDict[orthogroupID] = []
                                    orthogroupDict[orthogroupID].append(geneID)
    #### start here
    compiledData = {}
    for orthogroupID in orthogroupDict:
        if orthogroupID in msaDataDict:
            for geneID in orthogroupDict[orthogroupID]:
                if geneID in msaDataDict[orthogroupID]:
                    genomeID,proteinAlignment = msaDataDict[orthogroupID][geneID]
                    if genomeID in populationDict:
                        populationID = populationDict[genomeID]
                        # print(orthogroupID,geneID,populationID)
                        
                        if orthogroupID not in orthogroupProteinCount:
                            orthogroupProteinCount[orthogroupID] = {}
                        if populationID not in orthogroupProteinCount[orthogroupID]:
                            orthogroupProteinCount[orthogroupID][populationID] = {}
                        if geneID not in orthogroupProteinCount[orthogroupID][populationID]:
                            orthogroupProteinCount[orthogroupID][populationID][geneID] = (genomeID,proteinAlignment)
                            
                        if orthogroupID not in compiledData:
                            compiledData[orthogroupID] = {}
                        if populationID not in compiledData[orthogroupID]:
                            compiledData[orthogroupID][populationID] = []
                        compiledData[orthogroupID][populationID].append(proteinAlignment)
                        # print(proteinAlignment)
    # https://stackoverflow.com/questions/59118402/how-to-create-multiple-sequence-alignments-with-fasta-files-rather-then-strings
    OUT = open('population_average_entropies_' + outputLabel + '.txt','w')
    OUT.write("### orthogroupID\togProteinCount\tpopulationID\tgenomeID\tgeneID\taverageEntropy\tannotation\n")
    OUT.write("### Max theoretical entropy value is 4.3 (no similarity between sequences)\n")
    OUT.write("### Entropy equal to zero corresponds to no difference between sequences\n")
    averageEntropies = {}
    for orthogroupID in compiledData:
        for populationID in compiledData[orthogroupID]:
            newAlignment = MultipleSeqAlignment(compiledData[orthogroupID][populationID])
            average_entropy = compute_alignment_entropy(newAlignment)
            if orthogroupID not in averageEntropies:
                averageEntropies[orthogroupID] = {}
            if populationID not in averageEntropies[orthogroupID]:
                averageEntropies[orthogroupID][populationID] = average_entropy

    for orthogroupID in orthogroupProteinCount:
        for populationID in orthogroupProteinCount[orthogroupID]:
            ogNumberProteins = len(orthogroupProteinCount[orthogroupID][populationID])
            if orthogroupID in averageEntropies and populationID in averageEntropies[orthogroupID]:
                average_entropy = averageEntropies[orthogroupID][populationID]
                # average_entropy = round(average_entropy,3)
            else:
                print(orthogroupID,populationID,'not in averageEntropies')
                sys.exit()
            for geneID in orthogroupProteinCount[orthogroupID][populationID]:
                genomeID,proteinAlignment = orthogroupProteinCount[orthogroupID][populationID][geneID]
                if geneID in summaryAnnotDict:
                    annotation,fractionMasked = summaryAnnotDict[geneID]
                else:
                    print(geneID,'not in summaryAnnotDict')
                    sys.exit()
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (orthogroupID,ogNumberProteins,populationID,genomeID,geneID,average_entropy,annotation))

            
    '''
    for indexID in notCannabisIndices:
        genomeID = notCannabisIndices[indexID]
        print(indexID,genomeID)
    '''


def compute_alignment_entropy(alignment):
    num_columns = alignment.get_alignment_length()
    entropies = []
    for i in range(num_columns):
        # print(alignment[i])
        column = alignment[:, i]
        columnSet = set(column)
        # {'-'}
        # print(column)
        if '-' in columnSet and len(columnSet) == 1:
            continue
        else:
            # print(column)
            column_entropy = compute_column_entropy(column)
            entropies.append(column_entropy)
    # print(entropies)
    # Average entropy across columns
    return np.mean(entropies)


def compute_column_entropy(column):
    # Count the frequencies of characters in the column
    counts = Counter(column)
    
    # Get the frequencies as an array
    frequencies = np.array(list(counts.values()))
    
    # Normalize the frequencies to get probabilities
    probabilities = frequencies / sum(frequencies)
    
    # Compute entropy using the base-2 logarithm
    return entropy(probabilities, base=2)


def readChemotypeFile(chemotypeFile):
    #Abv,Chemotype,Type
    #79X,1,mj
    chemotypeDict = {}
    populationDict = {}
    with open(chemotypeFile,'r') as F:
        for line in F:
            if 'Chemotype' not in line:
                assemblyID,chemotypeID,populationID = line.strip().split(',')
                newChemotypeID = 'type' + chemotypeID
                if assemblyID not in chemotypeDict:
                    chemotypeDict[assemblyID] = newChemotypeID
                if assemblyID not in populationDict:
                    populationDict[assemblyID] = populationID
    return(chemotypeDict,populationDict)

                    
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <chemotype file> <orthogroup TSV file> <msa file list> <summary file list> <output label> \n"
if len(sys.argv) != 6:
    print(usage)
    sys.exit()

chemotypeFile = sys.argv[1]
orthogroupTSV = sys.argv[2]
msaFileList = sys.argv[3]
summaryFileList = sys.argv[4]
outputLabel = sys.argv[5]

notCannabis = {'FragariaVesca':'1', 'LotusJaponicus':'1', 'MalusDomestica':'1', 'PrunusPersica':'1', 'RosaChinensis':'1'}

summaryDataList = readFileList(summaryFileList)
summaryAnnotDict = readSummaryFile(summaryDataList)
                              
chemotypeDict,populationDict = readChemotypeFile(chemotypeFile)
msaDataList = readFileList(msaFileList)
msaDataDict = readFastaFile(msaDataList,notCannabis)
readOrthogroupTSV(orthogroupTSV,msaDataDict,populationDict,notCannabis,summaryAnnotDict,outputLabel)

B
