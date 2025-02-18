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
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


###############
# SUBROUTINES #
###############

### Entropy equal to zero corresponds to no difference between sequences
# orthogroupID        ogProteinCount  populationID    genomeID        geneID  averageEntropy  annotation
# OG0000000       379     asian_hemp      CNBv1a  CNBv1a.000459.v1.g351550.t1     0.34479801927867987     Encoded by
def readEntropyGeneFile(entropyGeneFile,specificOrthogroupID):
    averageEntropyDict = {}
    with open(entropyGeneFile,'r') as F:
        for line in F:
            if '#' not in line:
                orthogroupID,ogProteinCount,populationID,genomeID,geneID,averageEntropy,annotation = line.strip().split('\t')
                averageEntropy = float(averageEntropy)
                if orthogroupID == specificOrthogroupID:
                    if orthogroupID not in averageEntropyDict:
                        averageEntropyDict[orthogroupID] = {}
                    if populationID not in averageEntropyDict[orthogroupID]:
                        averageEntropyDict[orthogroupID][populationID] = averageEntropy
    return(averageEntropyDict)
                

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


def readFastaFile(fastaFile,notCannabis):
    recordDict = {}
    fileBase = os.path.basename(fastaFile)
    orthogroupID,fileSuffix = os.path.splitext(fileBase)
    # print(orthogroupID)
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
'''
def readOrthogroupTSV(orthogroupTSV,msaDataDict,populationDict,notCannabis,summaryAnnotDict,outputLabel,specificOrthogroupID):
    orthogroupDict = {}
    orthogroupProteinCount = {}
    with open(orthogroupTSV,'r') as F:
        for line in F:
            orthoGroupInfo = line.strip().split('\t')
            if not line.startswith('Orthogroup'):
                orthogroupID = orthoGroupInfo[0]
                if orthogroupID == specificOrthogroupID:
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
                            geneID = item.strip()
                            if '.' in geneID:
                                cols = geneID.split('.')
                                assemblyID = cols[0]
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
    averageEntropies = {}
    for orthogroupID in compiledData:
        for populationID in compiledData[orthogroupID]:
            newAlignment = MultipleSeqAlignment(compiledData[orthogroupID][populationID])
            average_entropy_dict = compute_alignment_entropy(newAlignment)
            if orthogroupID not in averageEntropies:
                averageEntropies[orthogroupID] = {}
            if populationID not in averageEntropies[orthogroupID]:
                averageEntropies[orthogroupID][populationID] = average_entropy_dict

    # orthogroupProteinCount[orthogroupID][populationID][geneID] = (genomeID,proteinAlignment)
    plotDict = {}
    for orthogroupID in orthogroupProteinCount:
        if orthogroupID not in plotDict:
            plotDict[orthogroupID] = {}
        for populationID in orthogroupProteinCount[orthogroupID]:
            if populationID not in plotDict[orthogroupID]:
                plotDict[orthogroupID][populationID] = {}
            if 'xVals' not in plotDict[orthogroupID][populationID]:
                plotDict[orthogroupID][populationID]['xVals'] = []
            if 'yVals' not in plotDict[orthogroupID][populationID]:
                plotDict[orthogroupID][populationID]['yVals'] = []
            ogNumberProteins = len(orthogroupProteinCount[orthogroupID][populationID])
            if orthogroupID in averageEntropies and populationID in averageEntropies[orthogroupID]:
                average_entropy_dict = averageEntropies[orthogroupID][populationID]
                # entropies[i] = column_entropy # dict
                for pos in average_entropy_dict:
                    plotDict[orthogroupID][populationID]['xVals'].append(pos)
                    entropy_value = average_entropy_dict[pos]
                    plotDict[orthogroupID][populationID]['yVals'].append(entropy_value)
                    # print(orthogroupID,populationID,pos,entropy_value)
            else:
                print(orthogroupID,populationID,'not in averageEntropies')
                sys.exit()
            for geneID in orthogroupProteinCount[orthogroupID][populationID]:
                genomeID,proteinAlignment = orthogroupProteinCount[orthogroupID][populationID][geneID]
                # summaryAnnotDict[transcriptID] = (annotation,fractionMasked)
                if geneID in summaryAnnotDict:
                    annotation,fractionMasked = summaryAnnotDict[geneID]
                else:
                    print(geneID,'not in summaryAnnotDict')
                    sys.exit()
            
    '''
    for indexID in notCannabisIndices:
        genomeID = notCannabisIndices[indexID]
        print(indexID,genomeID)
    '''
    return(plotDict,orthogroupProteinCount)



def plotEntropies(plotDict,orthogroupProteinCount,averageEntropyDict,outputLabel):
    # swap MJ and F1
    cleanLabels = {'asian_hemp':'Asian hemp', 'feral':'feral', 'F1':'F1', 'hc_hemp':'HC hemp', 'hemp':'hemp', 'mj':'MJ'}
    colorDict = {'asian_hemp':'#d1e5f0', 'feral':'#fddbc7', 'F1':'#67001f', 'hc_hemp':'#b2abd2', 'hemp':'#4393c3', 'mj':'#d6604d'}
    for orthogroupID in plotDict:
        fig, ax1 = plt.subplots(figsize=(20,5))
        legendList = []
        for populationID in plotDict[orthogroupID]:
            xValues = plotDict[orthogroupID][populationID]['xVals']
            yValues = plotDict[orthogroupID][populationID]['yVals']
            if orthogroupID in orthogroupProteinCount and populationID in orthogroupProteinCount[orthogroupID]:
                ogNumberProteins = len(orthogroupProteinCount[orthogroupID][populationID])
            else:
                print(orthogroupID,populationID,'not in orthogroupProteinCount')
                sys.exit()
            if orthogroupID in averageEntropyDict and populationID in averageEntropyDict[orthogroupID]:
                averageEntropy = averageEntropyDict[orthogroupID][populationID]
                if averageEntropy > 0.01:
                    averageEntropy = round(averageEntropy,2)
                else:
                    averageEntropy = round(averageEntropy,5)
            else:
                print(orthogroupID,populationID,'not in averageEntropyDict')
                sys.exit()
            ax1.plot(xValues, yValues, color=colorDict[populationID], label=cleanLabels[populationID], linewidth=4)
            # ax2 = ax1.twinx()
            # ax2.plot(xValues, yValues2, color='blue', label='')
            ax1.set_ylabel('Entropy')
            #ax2.set_ylabel('Occupancy (Position count)')
            ax1.set_xlabel('Position in amino acid alignment')
            legendInfo = mlines.Line2D([], [], color=colorDict[populationID], marker='s',linestyle="None",markersize=12, label=cleanLabels[populationID] + '\n' + str(ogNumberProteins) + ' proteins\nAverage entropy ' + str(averageEntropy), alpha=1.0, fillstyle='full', markeredgecolor=colorDict[populationID], markeredgewidth=0.0)
            # legendInfo2 = mlines.Line2D([], [], color='blue', marker='s',linestyle="None",markersize=6, label='Occupancy (Position count)', alpha=1.0, fillstyle='full', markeredgecolor='blue', markeredgewidth=0.0)
            legendList.append(legendInfo)
        plt.legend(handles=legendList, frameon=False, ncols=6, loc='upper center')
        plt.savefig(orthogroupID + '_' + outputLabel + '.png', bbox_inches='tight', dpi=600)
        plt.savefig(orthogroupID + '_' + outputLabel + '.svg', bbox_inches='tight')
        print(orthogroupID + '_' + outputLabel + '.png')
        plt.close()


def compute_alignment_entropy(alignment):
    num_columns = alignment.get_alignment_length()
    entropies = {}
    for i in range(num_columns):
        column = alignment[:, i]
        columnSet = set(column)
        column_entropy = compute_column_entropy(column)
        if i not in entropies:
            entropies[i] = column_entropy
    return entropies


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

usage = "Usage: " + sys.argv[0] + " <chemotype file> <orthogroup file> <msa file> <summary file list> <entropy-gene file> <specific orthogroup ID> <output label> \n"
if len(sys.argv) != 8:
    print(usage)
    sys.exit()

chemotypeFile = sys.argv[1]
orthogroupTSV = sys.argv[2]
msaFile = sys.argv[3]
summaryFileList = sys.argv[4]
entropyGeneFile = sys.argv[5]
specificOrthogroupID = sys.argv[6]
outputLabel = sys.argv[7]

notCannabis = {'FragariaVesca':'1', 'LotusJaponicus':'1', 'MalusDomestica':'1', 'PrunusPersica':'1', 'RosaChinensis':'1'}

averageEntropyDict = readEntropyGeneFile(entropyGeneFile,specificOrthogroupID)

summaryDataList = readFileList(summaryFileList)
summaryAnnotDict = readSummaryFile(summaryDataList)
                              
chemotypeDict,populationDict = readChemotypeFile(chemotypeFile)
# msaDataList = readFileList(msaFileList)
msaDataDict = readFastaFile(msaFile,notCannabis)
plotDict,orthogroupProteinCount = readOrthogroupTSV(orthogroupTSV,msaDataDict,populationDict,notCannabis,summaryAnnotDict,outputLabel,specificOrthogroupID)

plotEntropies(plotDict,orthogroupProteinCount,averageEntropyDict,outputLabel)
