import os,re,sys,math
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches


def createHistogram(histTextFile,outputLabel,sexID):
    divisor = 1000000
    values = {}
    syntenyColorList = {}
    colorMap = {'NotSyntenic':'black', 'bothXY':'blue', 'onlyChrX':'red', 'onlyChrY':'orange',
            'syntenicAutosomeAllHaplotypes':'violet', 'syntenicAutosomeFemaleHaplotypes':'green',
            'syntenicAutosomeMaleHaplotypes':'cyan'}
    syntenyLabels = {}
    with open(histTextFile,'r') as F:
        for line in F:
            if 'femaleAvgTPM' not in line:
                line = line.strip().split('\t')
                geneID = line[1]
                items = geneID.split('.')
                assemblyID = items[0]
                chrID = items[1]
                pos = line[2]
                pos = int(pos)
                pos = float(pos) / divisor
                syntenyID = line[4]
                femaleAvgTPM = line[6]
                femaleAvgTPM = float(femaleAvgTPM)
                # math.log(expValue+1,2)
                femaleAvgLogTPM = math.log(femaleAvgTPM+1,2)
                maleAvgTPM = line[8]
                maleAvgTPM = float(maleAvgTPM)
                maleAvgLogTPM = math.log(maleAvgTPM+1,2)
                tpmLogDiff = maleAvgLogTPM - femaleAvgLogTPM
                tpmDiff = line[9]
                tpmDiff = float(tpmDiff)
                if assemblyID not in syntenyColorList:
                    syntenyColorList[assemblyID] = {}
                if chrID not in syntenyColorList[assemblyID]:
                    syntenyColorList[assemblyID][chrID] = []
                syntenyColorList[assemblyID][chrID].append(colorMap[syntenyID])
                if assemblyID not in syntenyLabels:
                    syntenyLabels[assemblyID] = {}
                if chrID not in syntenyLabels[assemblyID]:
                    syntenyLabels[assemblyID][chrID] = {}
                if syntenyID not in syntenyLabels[assemblyID][chrID]:
                    syntenyLabels[assemblyID][chrID][syntenyID] = 1
                
                if assemblyID not in values:
                    values[assemblyID] = {}
                if chrID not in values[assemblyID]:
                    values[assemblyID][chrID] = {}
                if syntenyID not in values[assemblyID][chrID]:
                    values[assemblyID][chrID][syntenyID] = []
                values[assemblyID][chrID][syntenyID].append((pos,maleAvgTPM,maleAvgLogTPM,tpmDiff,tpmLogDiff))

    allValues = {}
    for assemblyID in values:
        for chrID in values[assemblyID]:
            for syntenyID in values[assemblyID][chrID]:
                posList,maleAvgTPM,maleAvgLogTPM,tpmDiffList,tpmLogDiff = zip(*values[assemblyID][chrID][syntenyID])
                if assemblyID not in allValues:
                    allValues[assemblyID] = {}
                if chrID not in allValues[assemblyID]:
                    allValues[assemblyID][chrID] = []
                for i in tpmLogDiff:
                    allValues[assemblyID][chrID].append(i)

    allPos = {}
    for assemblyID in values:
        for chrID in values[assemblyID]:
            for syntenyID in values[assemblyID][chrID]:
                posList,maleAvgTPM,maleAvgLogTPM,tpmDiffList,tpmLogDiff = zip(*values[assemblyID][chrID][syntenyID])
                if assemblyID not in allPos:
                    allPos[assemblyID] = {}
                if chrID not in allPos[assemblyID]:
                    allPos[assemblyID][chrID] = []
                for i in posList:
                    allPos[assemblyID][chrID].append(i)

    newLabels = {'onlyChrX':'X-specific region', 'onlyChrY':'Sex-determining region (SDR)', 'bothXY':'Pseudo-autosomal region (PAR)', 'NotSyntenic':'Not syntenic'}
    binWidth = 0.5
    posBinWidth = 1
    for assemblyID in values:
        for chrID in values[assemblyID]:
            posBinList = allPos[assemblyID][chrID]
            posBins = np.arange(min(posBinList),max(posBinList),posBinWidth)
            dataList = allValues[assemblyID][chrID]
            bins = np.arange(min(dataList),max(dataList),binWidth)
            handleIDs = []
            fig = plt.figure()
            gs = GridSpec(4,4)
            ax_joint = fig.add_subplot(gs[1:4,0:3])
            ax_marg_x = fig.add_subplot(gs[0,0:3])
            ax_marg_y = fig.add_subplot(gs[1:4,3])
            for syntenyID in values[assemblyID][chrID]:
                posList,maleAvgTPM,maleAvgLogTPM,tpmDiffList,tpmLogDiff = zip(*values[assemblyID][chrID][syntenyID])
                if syntenyID in newLabels:
                    newLabel = newLabels[syntenyID]
                    patchID = mpatches.Patch(color=colorMap[syntenyID], label=newLabel)
                else:
                    patchID = mpatches.Patch(color=colorMap[syntenyID], label=syntenyID)
                handleIDs.append(patchID)
                ax_joint.scatter(posList,tpmLogDiff, color=colorMap[syntenyID], label=syntenyID, alpha=0.5, edgecolor='none')
                ax_marg_x.hist(posList, color=colorMap[syntenyID], alpha=1.0, histtype='step', bins=posBins, linewidth=2)
                ax_marg_y.hist(tpmLogDiff, orientation="horizontal", color=colorMap[syntenyID], density=True, cumulative=-1, alpha=1, histtype='step', bins = bins, linewidth=2)
            plt.setp(ax_marg_x.get_xticklabels(), visible=False)
            plt.setp(ax_marg_y.get_yticklabels(), visible=False)

            ax_joint.set_xlabel('Gene start positions (Mb) on ' + chrID)
            ax_joint.set_ylabel('Difference of log' + r'$_{2}$' + ' average TPM\nin genes with ' + sexID + '-biased expression')
            
            ax_marg_y.set_xlabel('Reversed\ncumulative\nhistogram')
            ax_marg_x.set_ylabel('Frequency\nhistogram')
            ax_joint.legend(handles=handleIDs)
            plt.tight_layout()
            plt.title(chrID)
            plt.savefig(assemblyID + '_' + chrID + '_' + outputLabel + '_pos_hist.png', dpi=600)
            print(assemblyID + '_' + chrID + '_' + outputLabel + '_pos_hist.png')
            plt.close()


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hist txt file> <output label> <sex ID>\n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

histTextFile = sys.argv[1]
outputLabel = sys.argv[2]
sexID = sys.argv[3]

createHistogram(histTextFile,outputLabel,sexID)
