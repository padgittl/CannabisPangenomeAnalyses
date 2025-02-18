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

def ogAnnotationFile(ogAnnotFile):
    annotDict = {}
    with open(ogAnnotFile,'r') as F:
        for line in F:
            orthogroupID,annotation = line.strip().split('\t')
            if orthogroupID not in annotDict:
                annotDict[orthogroupID] = annotation
    return(annotDict)


def readGenomeIDFile(genomeIDFile):
    pangenomeIDs = {}
    with open(genomeIDFile,'r') as F:
        for line in F:
            genomeID = line.strip()
            if genomeID not in pangenomeIDs:
                pangenomeIDs[genomeID] = 1
    return(pangenomeIDs)

'''
FragariaVesca
LotusJaponicus
MalusDomestica
PrunusPersica
RosaChinensis
'''
'''
soft core, shell, and cloud genes which are present in 95–98%, 5–94% and 1–5% of individuals

pangenomeCount = len(pangenomeIDs)
coreThreshold = pangenomeCount * 0.99
lowerSoftCore = pangenomeCount * 0.95
upperSoftCore = pangenomeCount * 0.98
lowerShell = pangenomeCount * 0.05
upperShell = pangenomeCount * 0.94
lowerCloud = pangenomeCount * 0.01
upperCloud = pangenomeCount * 0.05

pangenomeCount 193
coreThreshold 191.07
lowerSoftCore 183.35
upperSoftCore 189.14
lowerShell 9.65
upperShell 181.42
lowerCloud 1.93
upperCloud 9.65
'''
def countOGs(orthogroupDict,populationDict,pangenomeIDs,annotDict,geneDict,unassignedGroups,unassignedGeneDict,scaffoldedGenomeIDs,publicScaffoldedGenomeIDFile):
    # https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-024-10931-w
    groupOGs = {}
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
    print("pangenomeCount",pangenomeCount)
    print("coreThreshold",int(coreThreshold))
    print("lowerSoftCore",int(lowerSoftCore))
    print("upperSoftCore",int(upperSoftCore))
    print("lowerShell",lowerShell)
    print("upperShell",upperShell)
    print("lowerCloud",lowerCloud)
    print("upperCloud",upperCloud)
    #  orthogroupDict[orthogroupID][assemblyID].append(geneSet)
    '''
    OUT = open('core_shell_variable_designations.txt','w')
    CORE = open('core_genes.txt','w')
    SOFTCORE = open('softcore_genes.txt','w')
    SHELL = open('shell_genes.txt','w')
    CLOUD = open('cloud_genes.txt','w')
    PRIVATE = open('private_genes.txt','w')
    '''
    print("OUTPUT FILES")
    print("core_shell_variable_designations.txt")
    print("core_genes.txt")
    print("softcore_genes.txt")
    print("shell_genes.txt")
    print("cloud_genes.txt")
    print("private_genes.txt")
    '''
    # annotDict[orthogroupID] = annotation
    # geneDict[orthogroupID].append(geneSet)
    assemblyDict = {}
    '''
    print("OUTPUT FILES")
    print("core_shell_variable_designations.txt")
    print("core_genes.txt")
    print("softcore_genes.txt")
    print("shell_genes.txt")
    print("cloud_genes.txt")
    print("private_genes.txt")
    # confirm these occur once, yes
    '''
    # assign unassigned genes to private group
    for orthogroupID in unassignedGroups:
        genomeCount = len(unassignedGroups[orthogroupID])
        # there are examples with zero genes, these are non-cannabis and would get filtered out
        if genomeCount > 1:
            print(orthogroupID,genomeCount,'more than one gene in unassigned OG')
            sys.exit()
        if orthogroupID in annotDict:
            annotation = annotDict[orthogroupID]
        else:
            annotation = 'NoAnnotationReported'
        if genomeCount <= privateThreshold:
            if 'private' not in groupOGs:
                groupOGs['private'] = {}
            if orthogroupID not in groupOGs['private']:
                groupOGs['private'][orthogroupID] = genomeCount
            OUT.write("%s\t%s\t%s\t%s\n" % ('PRIVATE',orthogroupID,genomeCount,annotation))
            if orthogroupID in unassignedGeneDict:
                for geneID in unassignedGeneDict[orthogroupID]:
                    PRIVATE.write("%s%s%s\n" % ('"',geneID,'"'))
                    geneID = geneID.strip()
                    items = geneID.split('.')
                    assemblyID = items[0]
                    if assemblyID not in assemblyDict:
                        assemblyDict[assemblyID] = {}
                    if 'private' not in assemblyDict[assemblyID]:
                        assemblyDict[assemblyID]['private'] = {}
                    if orthogroupID not in assemblyDict[assemblyID]['private']:
                        assemblyDict[assemblyID]['private'][orthogroupID] = []
                    assemblyDict[assemblyID]['private'][orthogroupID].append(geneID)
            
    for orthogroupID in orthogroupDict:
        if orthogroupID in annotDict:
            annotation = annotDict[orthogroupID]
        else:
            annotation = 'NoAnnotationReported'
        genomeCount = len(orthogroupDict[orthogroupID])
        # print(orthogroupID,genomeCount)
        if genomeCount <= privateThreshold:
            # print(orthogroupID,genomeCount)
            if 'private' not in groupOGs:
                groupOGs['private'] = {}
            if orthogroupID not in groupOGs['private']:
                groupOGs['private'][orthogroupID] = genomeCount
            OUT.write("%s\t%s\t%s\t%s\n" % ('PRIVATE',orthogroupID,genomeCount,annotation))
            if orthogroupID in geneDict:
                for geneID in geneDict[orthogroupID]:
                    PRIVATE.write("%s%s%s\n" % ('"',geneID,'"'))
                    geneID = geneID.strip()
                    items = geneID.split('.')
                    assemblyID = items[0]
                    if assemblyID not in assemblyDict:
                        assemblyDict[assemblyID] = {}
                    if 'private' not in assemblyDict[assemblyID]:
                        assemblyDict[assemblyID]['private'] = {}
                    if orthogroupID not in assemblyDict[assemblyID]['private']:
                        assemblyDict[assemblyID]['private'][orthogroupID] = []
                    assemblyDict[assemblyID]['private'][orthogroupID].append(geneID)
        elif genomeCount >= lowerCloud and genomeCount <= upperCloud:
            if 'cloud' not in groupOGs:
                groupOGs['cloud'] = {}
            if orthogroupID not in groupOGs['cloud']:
                groupOGs['cloud'][orthogroupID] = genomeCount
            OUT.write("%s\t%s\t%s\t%s\n" % ('CLOUD',orthogroupID,genomeCount,annotation))
            if orthogroupID in geneDict:
                for geneID in geneDict[orthogroupID]:
                    CLOUD.write("%s%s%s\n" % ('"',geneID,'"'))
                    geneID = geneID.strip()
                    items = geneID.split('.')
                    assemblyID = items[0]
                    if assemblyID not in assemblyDict:
                        assemblyDict[assemblyID] = {}
                    if 'cloud' not in assemblyDict[assemblyID]:
                        assemblyDict[assemblyID]['cloud'] = {}
                    if orthogroupID not in assemblyDict[assemblyID]['cloud']:
                        assemblyDict[assemblyID]['cloud'][orthogroupID] = []
                    assemblyDict[assemblyID]['cloud'][orthogroupID].append(geneID)
        elif genomeCount >= lowerShell and genomeCount <= upperShell:
            if 'shell' not in groupOGs:
                groupOGs['shell'] = {}
            if orthogroupID not in groupOGs['shell']:
                groupOGs['shell'][orthogroupID]	= genomeCount
            OUT.write("%s\t%s\t%s\t%s\n" % ('SHELL',orthogroupID,genomeCount,annotation))
            if orthogroupID in geneDict:
                for geneID in geneDict[orthogroupID]:
                    SHELL.write("%s%s%s\n" % ('"',geneID,'"'))
                    geneID = geneID.strip()
                    items = geneID.split('.')
                    assemblyID = items[0]
                    if assemblyID not in assemblyDict:
                        assemblyDict[assemblyID] = {}
                    if 'shell' not in   assemblyDict[assemblyID]:
                        assemblyDict[assemblyID]['shell'] = {}
                    if orthogroupID not in assemblyDict[assemblyID]['shell']:
                        assemblyDict[assemblyID]['shell'][orthogroupID] = []
                    assemblyDict[assemblyID]['shell'][orthogroupID].append(geneID)
        elif genomeCount >= lowerSoftCore and genomeCount <= upperSoftCore:
            if 'softCore' not in groupOGs:
                groupOGs['softCore'] = {}
            if orthogroupID not in groupOGs['softCore']:
                groupOGs['softCore'][orthogroupID] = genomeCount
            OUT.write("%s\t%s\t%s\t%s\n" % ('SOFTCORE',orthogroupID,genomeCount,annotation))
            if orthogroupID in geneDict:
                for geneID in geneDict[orthogroupID]:
                    SOFTCORE.write("%s%s%s\n" % ('"',geneID,'"'))
                    geneID = geneID.strip()
                    items = geneID.split('.')
                    assemblyID = items[0]
                    if assemblyID not in assemblyDict:
                        assemblyDict[assemblyID] = {}
                    if 'softcore' not in   assemblyDict[assemblyID]:
                        assemblyDict[assemblyID]['softcore'] = {}
                    if orthogroupID not in assemblyDict[assemblyID]['softcore']:
                        assemblyDict[assemblyID]['softcore'][orthogroupID] = []
                    assemblyDict[assemblyID]['softcore'][orthogroupID].append(geneID)
        elif genomeCount == coreThreshold:
            if 'core' not in groupOGs:
                groupOGs['core'] = {}
            if orthogroupID not in groupOGs['core']:
                groupOGs['core'][orthogroupID] = genomeCount
            OUT.write("%s\t%s\t%s\t%s\n" % ('CORE',orthogroupID,genomeCount,annotation))
            # print('CORE',orthogroupID,genomeCount,annotation)
            if orthogroupID in geneDict:
                for geneID in geneDict[orthogroupID]:
                    CORE.write("%s%s%s\n" % ('"',geneID,'"'))
                    geneID = geneID.strip()
                    items = geneID.split('.')
                    assemblyID = items[0]
                    if assemblyID not in assemblyDict:
                        assemblyDict[assemblyID] = {}
                    if 'core' not in   assemblyDict[assemblyID]:
                        assemblyDict[assemblyID]['core'] = {}
                    if orthogroupID not in assemblyDict[assemblyID]['core']:
                        assemblyDict[assemblyID]['core'][orthogroupID] = []
                    assemblyDict[assemblyID]['core'][orthogroupID].append(geneID)
        else:
            print("check count",orthogroupID,genomeCount)
            sys.exit()
    geneCounts = {}
    for assemblyID in assemblyDict:
        for geneGroup in assemblyDict[assemblyID]:
            for orthogroupID in assemblyDict[assemblyID][geneGroup]:
                geneCount = len(assemblyDict[assemblyID][geneGroup][orthogroupID])
                # print(assemblyID,geneGroup,ogCount)
                if assemblyID not in geneCounts:
                    geneCounts[assemblyID] = {}
                if geneGroup not in geneCounts[assemblyID]:
                    geneCounts[assemblyID][geneGroup] = 0
                geneCounts[assemblyID][geneGroup] += geneCount
    totalGeneCountDict = {}
    for assemblyID in geneCounts:
        if assemblyID not in totalGeneCountDict:
            totalGeneCountDict[assemblyID] = 0
        for geneGroup in geneCounts[assemblyID]:
            geneCount = geneCounts[assemblyID][geneGroup]
            totalGeneCountDict[assemblyID] += geneCount
    for assemblyID in assemblyDict:
        for geneGroup in assemblyDict[assemblyID]:
            ogCount = len(assemblyDict[assemblyID][geneGroup])
            if assemblyID in geneCounts and geneGroup in geneCounts[assemblyID]:
                geneCount = geneCounts[assemblyID][geneGroup]
                if assemblyID in totalGeneCountDict:
                    totalGeneCount = totalGeneCountDict[assemblyID]
                else:
                    print(assemblyID,'not in totalGeneCountDict')
                    sys.exit()
                print(assemblyID,geneGroup,ogCount,geneCount,totalGeneCount)

                
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


def readUnassignedGenesFile(unassignedGenesFile,pangenomeIDs,scaffoldedGenomeIDs,publicScaffoldedGenomeIDFile,contigIDsFromEH23a):
    unassignedGroups = {}
    unassignedGeneDict = {}
    with open(unassignedGenesFile,'r') as F:
        for line in F:
            if 'Orthogroup' not in line:
                orthoGroupInfo = line.strip().split('\t')
                orthogroupID = orthoGroupInfo[0]
                if orthogroupID not in unassignedGroups:
                    unassignedGroups[orthogroupID] = {}
                for item in orthoGroupInfo[1:]:
                    item = item.strip().split(',')
                    item = list(filter(str.strip, item))
                    for i in item:
                        i = i.strip()
                        if '.' in i:
                            cols = i.split('.')
                            if len(cols) == 5:
                                assemblyID = cols[0]
                                geneID = i
                                # print(geneID)
                                # 79X.00000017
                                # 79X.00000043
                                if assemblyID != 'F49v1a':
                                    if assemblyID in pangenomeIDs:
                                        items = geneID.split('.')
                                        genomeID = items[0]
                                        chrID = items[1]
                                        fullChromID = genomeID + '.' + chrID
                                        # contigIDsFromEH23a
                                        if assemblyID in scaffoldedGenomeIDs or assemblyID in publicScaffoldedGenomeIDFile:
                                            if 'chr' in geneID:
                                                if assemblyID not in unassignedGroups[orthogroupID]:
                                                    unassignedGroups[orthogroupID][assemblyID] = []
                                                unassignedGroups[orthogroupID][assemblyID].append(geneID)
                                                if orthogroupID not in unassignedGeneDict:
                                                    unassignedGeneDict[orthogroupID] = []
                                                unassignedGeneDict[orthogroupID].append(geneID)
                                        else:
                                            # these are examples where we need to filter the contig genomes
                                            if fullChromID in contigIDsFromEH23a:
                                                # print(fullChromID)
                                                if assemblyID not in unassignedGroups[orthogroupID]:
                                                    unassignedGroups[orthogroupID][assemblyID] = []
                                                unassignedGroups[orthogroupID][assemblyID].append(geneID)
                                                if orthogroupID not in unassignedGeneDict:
                                                    unassignedGeneDict[orthogroupID] = []
                                                unassignedGeneDict[orthogroupID].append(geneID)
    return(unassignedGroups,unassignedGeneDict)

                                        
def readOrthogroupTSV(orthogroupTSV,pangenomeIDs,scaffoldedGenomeIDs,publicScaffoldedGenomeIDFile,contigIDsFromEH23a):
    geneDict = {}
    orthogroupDict = {}
    orthogroupAssemblyCounts = {}
    with open(orthogroupTSV,'r') as F:
        for line in F:
            if not line.startswith('Orthogroup'):
                orthoGroupInfo = line.strip().split('\t')
                orthogroupID = orthoGroupInfo[0]
                if orthogroupID not in orthogroupDict:
                    orthogroupDict[orthogroupID] = {}
                if orthogroupID not in orthogroupAssemblyCounts:
                    orthogroupAssemblyCounts[orthogroupID] = {}
                for item in orthoGroupInfo[1:]:
                    if ',' in item:
                        item = item.strip().split(',')
                        # filter empty lines
                        item = list(filter(str.strip, item))
                        for i in item:
                            i = i.strip()
                            # 79X.00000625.v1.g046230.t1
                            if '.' in i:
                                cols = i.split('.')
                                if len(cols) == 5:
                                    assemblyID = cols[0]
                                    geneID = i
                                    if assemblyID in pangenomeIDs:
                                        items = geneID.split('.')
                                        genomeID = items[0]
                                        chrID = items[1]
                                        fullChromID = genomeID + '.' + chrID
                                        # print(pangenomeIDs)
                                        # these are scaffolded examples and so are filtered on basis of 10 chromosomes
                                        if assemblyID in scaffoldedGenomeIDs or assemblyID in publicScaffoldedGenomeIDFile:
                                            if 'chr' in geneID:
                                                if assemblyID not in orthogroupDict[orthogroupID]:
                                                    orthogroupDict[orthogroupID][assemblyID] = []
                                                orthogroupDict[orthogroupID][assemblyID].append(geneID)
                                                if orthogroupID not in geneDict:
                                                    geneDict[orthogroupID] = []
                                                geneDict[orthogroupID].append(geneID)
                                                if assemblyID not in orthogroupAssemblyCounts[orthogroupID]:
                                                    orthogroupAssemblyCounts[orthogroupID][assemblyID] = 0
                                                orthogroupAssemblyCounts[orthogroupID][assemblyID] += 1
                                        # these are contig-level and so will not have 10 chromosomes included
                                        else:
                                            if fullChromID in contigIDsFromEH23a:
                                                if assemblyID not in orthogroupDict[orthogroupID]:
                                                    orthogroupDict[orthogroupID][assemblyID] = []
                                                orthogroupDict[orthogroupID][assemblyID].append(geneID)
                                                if orthogroupID not in geneDict:
                                                    geneDict[orthogroupID] = []
                                                geneDict[orthogroupID].append(geneID)
                                                if assemblyID not in orthogroupAssemblyCounts[orthogroupID]:
                                                    orthogroupAssemblyCounts[orthogroupID][assemblyID] = 0
                                                orthogroupAssemblyCounts[orthogroupID][assemblyID] += 1
                    else:
                        geneID = item.strip()
                        if '.' in geneID:
                            items = geneID.split('.')
                            assemblyID = items[0]
                            chrID = items[1]
                            fullChromID = assemblyID + '.' + chrID
                            if assemblyID in pangenomeIDs:
                                if assemblyID in scaffoldedGenomeIDs or assemblyID in publicScaffoldedGenomeIDFile:
                                    if 'chr' in geneID:
                                        if assemblyID not in orthogroupDict[orthogroupID]:
                                            orthogroupDict[orthogroupID][assemblyID] = []
                                        orthogroupDict[orthogroupID][assemblyID].append(geneID)
                                        if orthogroupID not in geneDict:
                                            geneDict[orthogroupID] = []
                                        geneDict[orthogroupID].append(geneID)
                                        if assemblyID not in orthogroupAssemblyCounts[orthogroupID]:
                                            orthogroupAssemblyCounts[orthogroupID][assemblyID] = 0
                                        orthogroupAssemblyCounts[orthogroupID][assemblyID] += 1
                                else:
                                    if fullChromID in contigIDsFromEH23a:
                                        if assemblyID not in orthogroupDict[orthogroupID]:
                                            orthogroupDict[orthogroupID][assemblyID] = []
                                        orthogroupDict[orthogroupID][assemblyID].append(geneID)
                                        if orthogroupID not in geneDict:
                                            geneDict[orthogroupID] = []
                                        geneDict[orthogroupID].append(geneID)
                                        if assemblyID not in orthogroupAssemblyCounts[orthogroupID]:
                                            orthogroupAssemblyCounts[orthogroupID][assemblyID] = 0
                                        orthogroupAssemblyCounts[orthogroupID][assemblyID] += 1
                        
    '''
    for orthogroupID in orthogroupAssemblyCounts:
        for assemblyID in orthogroupAssemblyCounts[orthogroupID]:
            print(orthogroupID,assemblyID,orthogroupAssemblyCounts[orthogroupID][assemblyID])
    '''
    return(orthogroupDict,orthogroupAssemblyCounts,geneDict)

                    
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <chemotype file> <genome ID file> <og annotation file> <orthogroups tsv file> <unassigned gene/orthogroup file> <scaffolded (78) genome IDs> <public scaffolded genome ID file> <contig IDs from EH23a-jbrowse-paf> \n"
if len(sys.argv) != 9:
    print(usage)
    sys.exit()

chemotypeFile = sys.argv[1]
genomeIDFile = sys.argv[2]
ogAnnotFile = sys.argv[3]
orthogroupTSV = sys.argv[4]
unassignedGenesFile = sys.argv[5]
scaffoldedGenomeIDFile = sys.argv[6]
publicScaffoldedGenomeIDFile = sys.argv[7]
contigIDsFromEH23aPAF = sys.argv[8]

notCannabis = {'FragariaVesca':'1', 'LotusJaponicus':'1', 'MalusDomestica':'1', 'PrunusPersica':'1', 'RosaChinensis':'1'}
                              
annotDict = ogAnnotationFile(ogAnnotFile)
chemotypeDict,populationDict = readChemotypeFile(chemotypeFile)
pangenomeIDs = readGenomeIDFile(genomeIDFile)
scaffoldedGenomeIDs = readGenomeIDFile(scaffoldedGenomeIDFile)
contigIDsFromEH23a = readGenomeIDFile(contigIDsFromEH23aPAF)

orthogroupDict,orthogroupAssemblyCounts,geneDict = readOrthogroupTSV(orthogroupTSV,pangenomeIDs,scaffoldedGenomeIDs,publicScaffoldedGenomeIDFile,contigIDsFromEH23a)
unassignedGroups,unassignedGeneDict = readUnassignedGenesFile(unassignedGenesFile,pangenomeIDs,scaffoldedGenomeIDs,publicScaffoldedGenomeIDFile,contigIDsFromEH23a)

countOGs(orthogroupDict,populationDict,pangenomeIDs,annotDict,geneDict,unassignedGroups,unassignedGeneDict,scaffoldedGenomeIDs,publicScaffoldedGenomeIDFile)
