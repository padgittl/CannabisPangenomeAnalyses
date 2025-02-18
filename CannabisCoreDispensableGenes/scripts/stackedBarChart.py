import sys, re, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


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


def readOrderedGenomeIDFile(orderedGenomeIDFile):
    orderedGenomes = {}
    with open(orderedGenomeIDFile,'r') as F:
        for line in F:
            genomeID = line.strip()
            if genomeID not in orderedGenomes:
                orderedGenomes[genomeID] = 1
    return(orderedGenomes)

'''
ACBD private 762 765 42681
ACBD core 3360 10300 42681
ACBD softcore 13039 23352 42681
ACBD shell 4930 8070 42681
ACBD cloud 185 194 42681
'''
def readDataFileForPie(dataFile,orderedGenomes,populationDict):
    initData = {}
    with open(dataFile,'r') as F:
        for line in F:
            # 79X private 2965 2965 40864
            genomeID,groupID,ogCount,geneCount,totalGeneCount = line.strip().split()
            if groupID == 'private':
                groupID = 'Unique'
            elif groupID == 'softcore':
                groupID = 'Nearly core'
            else:
                groupID = groupID.capitalize()
            ogCount = int(ogCount)
            geneCount = int(geneCount)
            totalGeneCount = int(totalGeneCount)
            if genomeID not in initData:
                initData[genomeID] = []
            initData[genomeID].append((groupID,ogCount,geneCount,totalGeneCount))
    popData = {}
    for genomeID in orderedGenomes:
        if genomeID in initData:
            if genomeID in populationDict:
                populationID = populationDict[genomeID]
                if '?' not in populationID:
                    if populationID not in popData:
                        popData[populationID] = {}
                    for groupID,ogCount,geneCount,totalGeneCount in initData[genomeID]:
                        if groupID not in popData[populationID]:
                            popData[populationID][groupID] = 0
                        popData[populationID][groupID] += geneCount
        else:
            print(genomeID,'not in initData')
            sys.exit()
    return(popData)
    
'''
79X private 2965 2965 40864
79X core 3360 9170 40864
79X softcore 13054 20855 40864
79X shell 5025 7643 40864
79X cloud 224 231 40864

privateThreshold = 2
pangenomeCount = len(pangenomeIDs)
coreThreshold = pangenomeCount
lowerSoftCore = 183 # 183/193=0.95
upperSoftCore = 192 # 192/193=0.99
lowerShell = 10 # 10/193=0.05
upperShell = 182 # 182/193=0.94
lowerCloud = 3 # 3/193=0.02
upperCloud = 9 # 9/193=0.05
'''
def readDataFile(dataFile,orderedGenomes):
    fullData = {}
    labels = []
    privateList = []
    coreList = []
    softcoreList = []
    shellList = []
    cloudList = []
    # RGBA
    # 7d2e68ff magenta
    # f39237ff orange
    # 52ad9cff light green
    # 005e7cff teal
    # e181b3ff lighter pink
    initData = {}
    initGenomeIDs = []
    totalPrivate = 0
    totalShell = 0
    totalCloud = 0
    totalSoftCore = 0
    totalCore = 0
    with open(dataFile,'r') as F:
        for line in F:
            # 79X private 2965 2965 40864
            genomeID,groupID,ogCount,geneCount,totalGeneCount = line.strip().split()
            ogCount = int(ogCount)
            geneCount = int(geneCount)
            totalGeneCount = int(totalGeneCount)
            if genomeID not in initData:
                initData[genomeID] = []
            initData[genomeID].append((groupID,ogCount,geneCount,totalGeneCount))
            initGenomeIDs.append(genomeID)
    for genomeID in orderedGenomes:
        if genomeID in initData:
            for groupID,ogCount,geneCount,totalGeneCount in initData[genomeID]:
                if genomeID not in fullData:
                    fullData[genomeID] = {}
                labels.append(genomeID)
                if groupID == 'private':
                    totalPrivate += geneCount
                    privateList.append(geneCount)
                    # print(genomeID,groupID,ogCount,geneCount,totalGeneCount)
                    if 'PRIVATE' not in fullData[genomeID]:
                        fullData[genomeID]['PRIVATE'] = geneCount
                elif groupID == 'core' and 'soft' not in groupID:
                    totalCore += geneCount
                    coreList.append(geneCount)
                    # print(genomeID,groupID,ogCount,geneCount,totalGeneCount)
                    if 'CORE' not in fullData[genomeID]:
                        fullData[genomeID]['CORE'] = geneCount
                elif groupID == 'softcore':
                    totalSoftCore += geneCount
                    softcoreList.append(geneCount)
                    if 'SOFTCORE' not in fullData[genomeID]:
                        fullData[genomeID]['SOFTCORE'] = geneCount
                        # print(genomeID,groupID,ogCount,geneCount,totalGeneCount)
                elif groupID == 'shell':
                    totalShell += geneCount
                    shellList.append(geneCount)
                    # print(genomeID,groupID,ogCount,geneCount,totalGeneCount)
                    if 'SHELL' not in fullData[genomeID]:
                        fullData[genomeID]['SHELL'] = geneCount
                elif groupID == 'cloud':
                    totalCloud += geneCount
                    cloudList.append(geneCount)
                    #print(genomeID,groupID,ogCount,geneCount,totalGeneCount)
                    if 'CLOUD' not in fullData[genomeID]:
                        fullData[genomeID]['CLOUD'] = geneCount
                else:
                    print(genomeID,groupID,ogCount,geneCount,totalGeneCount)
                    sys.exit()
        else:
            print(genomeID,'not in initData')
            sys.exit()
    if len(labels) != len(privateList) != len(coreList) != len(softcoreList) != len(shellList) != len(cloudList):
        print("data lists not equal")
        sys.exit()
    return(labels,privateList,coreList,softcoreList,shellList,cloudList,fullData,
           totalPrivate,totalCore,totalSoftCore,totalShell,totalCloud)


def createPopulationPiechart(popData,keyword):
    coreColor = '#d1e5f0'
    nearlyCoreColor = '#053061'
    shellColor = '#f4a582'
    cloudColor = '#4393c3'
    uniqueColor = '#b2182b'
    overallPopTotals = {}
    for popID in popData:
        if popID not in overallPopTotals:
            overallPopTotals[popID] = 0
        for groupID in popData[popID]:
            popTotal = popData[popID][groupID]
            overallPopTotals[popID] += popTotal
    for popID in popData:
        sizes = []
        pie_labels = []
        plt.rcParams['font.size'] = 64
        overallPopTotal = overallPopTotals[popID]
        for groupID in popData[popID]:
            popTotal = popData[popID][groupID]
            percentValue = float(popTotal) / overallPopTotal * 100
            sizes.append(percentValue)
            pie_labels.append(groupID)
        plt.pie(sizes, startangle=90, radius=1.2, autopct='%1.1f%%', colors=[uniqueColor, coreColor, nearlyCoreColor, shellColor, cloudColor])
        plt.savefig(keyword + '_' + popID + '_pieChart.png', dpi=600)
        plt.savefig(keyword + '_' + popID + '_pieChart.svg')
        print(keyword + '_' + popID + '_pieChart.png')
        #plt.savefig('allPops_pieChart.png', dpi=600)
        plt.close()


def createPopulationBarChart(popData,keyword):
    coreColor = '#d1e5f0'
    nearlyCoreColor = '#053061'
    shellColor = '#f4a582'
    cloudColor = '#4393c3'
    uniqueColor = '#b2182b'
    overallPopTotals = {}
    groupData = {}
    popIDs = {}
    colorDict = {'Unique':'#b2182b', 'Core':'#d1e5f0', 'Nearly core':'#053061', 'Shell':'#f4a582', 'Cloud':'#4393c3'}
    for popID in popData:
        if popID not in overallPopTotals:
            overallPopTotals[popID] = 0
        for groupID in popData[popID]:
            popTotal = popData[popID][groupID]
            overallPopTotals[popID] += popTotal
    for popID in popData:
        overallPopTotal = overallPopTotals[popID]
        if popID not in popIDs:
            popIDs[popID] = 1
        for groupID in popData[popID]:
            popTotal = popData[popID][groupID]
            percentValue = float(popTotal) / overallPopTotal * 100
            percentValue = round(percentValue,2)
            if groupID not in groupData:
                groupData[groupID] = []
            groupData[groupID].append(percentValue)
    #print(popIDs)
    #print(groupData)
    x = np.arange(len(popIDs))
    # width = 0.25
    width = 0.1
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained')

    # print(groupData.items())
    for attribute, measurement in groupData.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute, color=colorDict[attribute])
        ax.bar_label(rects, padding=3)
        multiplier += 1

    ax.set_ylabel('Percent (%)')
    ax.set_title('Orthogroup membership by population')
    ax.set_xticks(x + width, popIDs)
    ax.legend(loc='center left', ncols=5)
    # ax.set_ylim(0, 250)
    plt.savefig(keyword + '_barchart.png',dpi=600)
    plt.savefig(keyword + '_barchart.svg',dpi=600)
    print(keyword + '_barchart.png')
    
        
def createPiechart(totalPrivate,totalCore,totalSoftCore,totalShell,totalCloud,keyword):
    # overlapping percent labels --> 
    # https://stackoverflow.com/questions/23577505/how-to-avoid-overlapping-of-labels-autopct-in-a-pie-chart

    coreColor = '#d1e5f0'
    nearlyCoreColor = '#053061'
    shellColor = '#f4a582'
    cloudColor = '#4393c3'
    uniqueColor = '#b2182b'
    
    overallTotal = totalPrivate + totalCore + totalSoftCore + totalShell + totalCloud
    
    totalPrivatePercent = float(totalPrivate) / overallTotal * 100
    totalCorePercent = float(totalCore) / overallTotal * 100
    totalSoftCorePercent = float(totalSoftCore) / overallTotal * 100
    totalShellPercent = float(totalShell) / overallTotal * 100
    totalCloudPercent = float(totalCloud) / overallTotal * 100
    
    sizes = [totalPrivatePercent, totalCorePercent, totalSoftCorePercent, totalShellPercent, totalCloudPercent]
    pie_labels = ['Unique', 'Core', 'Nearly core', 'Shell', 'Cloud']
    plt.rcParams['font.size'] = 64
    plt.pie(sizes, startangle=90, radius=1.2, autopct='%1.1f%%', colors=[uniqueColor, coreColor, nearlyCoreColor, shellColor, cloudColor])
    plt.savefig(keyword + '_pieChart.png', dpi=600)
    plt.savefig(keyword + '_pieChart.svg')
    print(keyword + '_pieChart.png')
    plt.close()



def plot(labels,privateList,coreList,softcoreList,shellList,cloudList,fullData):
    coreColor = '#d1e5f0'
    nearlyCoreColor = '#053061'
    shellColor = '#f4a582'
    cloudColor = '#4393c3'
    uniqueColor = '#b2182b'

    # x = labels
    x = []
    y1 = np.array(coreList)
    y2 = np.array(softcoreList)
    y3 = np.array(shellList)
    y4 = np.array(cloudList)
    y5 = np.array(privateList)
    barData = []
    for genomeID in fullData:
        x.append(genomeID)
        barData.append((genomeID,
                       fullData[genomeID]['CORE'],
                       fullData[genomeID]['SOFTCORE'],
                       fullData[genomeID]['SHELL'],
                       fullData[genomeID]['CLOUD'],
                       fullData[genomeID]['PRIVATE']))
    df = pd.DataFrame(barData, columns=['Genomes', 'Core', 'Nearly core', 'Shell', 'Cloud', 'Unique'])
    # print(df)
    plt.rcParams["figure.figsize"] = [40,20]
    plt.rc(('xtick.major', 'ytick.major'), width=2.5, size=24)
    plt.rcParams['font.size'] = 24

    df.plot(x='Genomes', kind='bar', stacked=True,
            title='Pangenome orthogroup membership', color=[coreColor, nearlyCoreColor, shellColor, cloudColor, uniqueColor])
    main_axes = plt.gca()
    main_axes.set_xticklabels(main_axes.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=12)
    '''
    plt.bar(x, y1, color='r')
    plt.bar(x, y2, bottom=y1, color='b')
    plt.bar(x, y3, bottom=y1+y2, color='y')
    plt.bar(x, y4, bottom=y1+y2+y3, color='g')
    plt.bar(x, y5, bottom=y1+y2+y3+y4, color='black')

    lowerSoftCore = 183 # 183/193=0.95
    upperSoftCore = 192 # 192/193=0.99
    lowerShell = 10 # 10/193=0.05
    upperShell = 182 # 182/193=0.94
    lowerCloud = 3 # 3/193=0.02
    upperCloud = 9 # 9/193=0.05
    '''
    
    plt.xlabel("Genome ID", fontsize=32)
    plt.ylabel("Gene counts", fontsize=32)
    plt.legend(['Core (193 genomes [100%])', 'Nearly core (183-192 [95-99%])', 'Shell (10-182 [5-94%])', 'Cloud (3-9 [2-5%])', 'Unique (1 or 2 genomes [0.5-1%])'], fontsize="24")
    plt.tight_layout()
    plt.savefig('stackedBarChart.png',dpi=600)
    plt.savefig('stackedBarChart.svg')
    print('stackedBarChart.png')
    plt.close()

                
usage = "Usage: " + sys.argv[0] + " <data> <ordered genome ID file> <chemotype/population file> <scaffolded genomes> <unscaffolded genomes> <public genomes> \n"
if len(sys.argv) != 7:
    print(usage)
    sys.exit()

dataFile = sys.argv[1]
orderedGenomeIDFile = sys.argv[2]
chemotypeFile = sys.argv[3]
scaffoldedGenomeFile = sys.argv[4]
unscaffoldedGenomeFile = sys.argv[5]
publicGenomeFile = sys.argv[6]

chemotypeDict,populationDict = readChemotypeFile(chemotypeFile)
orderedGenomes = readOrderedGenomeIDFile(orderedGenomeIDFile)

scaffoldedGenomes = readOrderedGenomeIDFile(scaffoldedGenomeFile)
unscaffoldedGenomes = readOrderedGenomeIDFile(unscaffoldedGenomeFile)
publicGenomes = readOrderedGenomeIDFile(publicGenomeFile)

labels,privateList,coreList,softcoreList,shellList,cloudList,fullData,totalPrivate,totalCore,totalSoftCore,totalShell,totalCloud = readDataFile(dataFile,orderedGenomes)
plot(labels,privateList,coreList,softcoreList,shellList,cloudList,fullData)
createPiechart(totalPrivate,totalCore,totalSoftCore,totalShell,totalCloud,'fullPangenome')

s_labels,s_privateList,s_coreList,s_softcoreList,s_shellList,s_cloudList,s_fullData,s_totalPrivate,s_totalCore,s_totalSoftCore,s_totalShell,s_totalCloud = readDataFile(dataFile,scaffoldedGenomes)
uns_labels,uns_privateList,uns_coreList,uns_softcoreList,uns_shellList,uns_cloudList,uns_fullData,uns_totalPrivate,uns_totalCore,uns_totalSoftCore,uns_totalShell,uns_totalCloud = readDataFile(dataFile,unscaffoldedGenomes)
pub_labels,pub_privateList,pub_coreList,pub_softcoreList,pub_shellList,pub_cloudList,pub_fullData,pub_totalPrivate,pub_totalCore,pub_totalSoftCore,pub_totalShell,pub_totalCloud = readDataFile(dataFile,publicGenomes)

# create pie chart for each of the assembly groupings that includes all populations -- this is to demonstrate similarity overall across methods
createPiechart(s_totalPrivate,s_totalCore,s_totalSoftCore,s_totalShell,s_totalCloud,'scaffolded')
createPiechart(uns_totalPrivate,uns_totalCore,uns_totalSoftCore,uns_totalShell,uns_totalCloud,'unscaffolded')
createPiechart(pub_totalPrivate,pub_totalCore,pub_totalSoftCore,pub_totalShell,pub_totalCloud,'publics')

allPopData = readDataFileForPie(dataFile,orderedGenomes,populationDict)

scaffoldedPieData = readDataFileForPie(dataFile,scaffoldedGenomes,populationDict)
unscaffoldedPieData = readDataFileForPie(dataFile,unscaffoldedGenomes,populationDict)
publicPieData = readDataFileForPie(dataFile,publicGenomes,populationDict)

createPopulationPiechart(allPopData,'allGenomes')
createPopulationPiechart(scaffoldedPieData,'chromLevelHapResolved')
createPopulationPiechart(unscaffoldedPieData,'contigScaffoldLevel')
createPopulationPiechart(publicPieData,'publics')

createPopulationBarChart(allPopData,'allGenomes')
createPopulationBarChart(scaffoldedPieData,'chromLevelHapResolved')
createPopulationBarChart(unscaffoldedPieData,'contigScaffoldLevel')
createPopulationBarChart(publicPieData,'publics')
