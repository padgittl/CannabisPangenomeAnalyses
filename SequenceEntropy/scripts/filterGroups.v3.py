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

'''
OG0000256
OG0000395
OG0000371
OG0000235
OG0000196
OG0000040
OG0000487
OG0000333
OG0000137
OG0000242
OG0000706
OG0000026

N3.HOG0003585 {'Reverse transcriptase (RNA-dependent DNA polymerase)': 1, 'DNA polymerase': 1}
N3.HOG0003590 {'Belongs to the ABC transporter superfamily. ABCG family. PDR (TC 3.A.1.205) subfamily': 1, 'Reverse transcriptase (RNA-dependent DNA polymerase)': 1}
N3.HOG0003720 {'ribonuclease H protein': 1, 'Ribonuclease H protein': 1}
N3.HOG0003720 {'ribonuclease H protein': 1, 'Ribonuclease H protein': 1}
N3.HOG0004057 {'Domain of unknown function (DUF4283)': 1, 'ribonuclease H protein': 1}
N3.HOG0004060 {'Domain of unknown function (DUF4283)': 1, 'Ribonuclease H protein': 1}
N3.HOG0004162 {'ribonuclease H protein At1g65750-like': 1, 'ribonuclease H protein': 1}
N3.HOG0004162 {'ribonuclease H protein At1g65750-like': 1, 'ribonuclease H protein': 1}
N3.HOG0004167 {'ribonuclease H protein At1g65750-like': 1, 'ribonuclease H protein': 1}
N3.HOG0004167 {'ribonuclease H protein At1g65750-like': 1, 'ribonuclease H protein': 1}
N3.HOG0004169 {'zinc-binding in reverse transcriptase': 1, 'ribonuclease H protein': 1}
N3.HOG0004174 {'zinc-binding in reverse transcriptase': 1, 'Ribonuclease H protein': 1}
'''

# OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (orthogroupID,hogID,hogNumberProteins,populationID,genomeID,geneID,average_entropy,annotation))
# orthogroupID        ogProteinCount  populationID    genomeID        geneID  averageEntropy  annotation
def readEntropyFile(entropyFile,confirmedTEs,repeatDescriptions,orthogroupSetID):
    dataDict = {}
    lines = {}
    with open(entropyFile,'r') as F:
        for line in F:
            if '#' not in line:
                orthogroupID,ogProteinCount,populationID,genomeID,geneID,average_entropy,annotation = line.strip().split('\t')
                if orthogroupID not in confirmedTEs:
                    if orthogroupID not in dataDict:
                        dataDict[orthogroupID] = {}
                    if annotation not in dataDict[orthogroupID]:
                        dataDict[orthogroupID][annotation] = 1
                    if orthogroupID not in lines:
                        lines[orthogroupID] = []
                    lines[orthogroupID].append((orthogroupID,ogProteinCount,populationID,genomeID,geneID,average_entropy,annotation))
                    
    notRepeat = {}
    repeat = {}
    for orthogroupID in dataDict:
        annotSet = set(dataDict[orthogroupID])
        for annotation in dataDict[orthogroupID]:
            # check if annotation is associated with TEs
            if annotation in repeatDescriptions:
                if orthogroupID not in repeat:
                    repeat[orthogroupID] = {}
                if annotation not in repeat[orthogroupID]:
                    repeat[orthogroupID][annotation] = 1
            # annotation is not associated with TEs
            else:
                if orthogroupID not in notRepeat:
                    notRepeat[orthogroupID] = {}
                if annotation not in notRepeat[orthogroupID]:
                    notRepeat[orthogroupID][annotation] = 1

    ogRemoval = {}
    for orthogroupID in dataDict:
        if orthogroupID in repeat:
            # cases where orthogroupID has TE and non-TE annotations -- BOTH
            if orthogroupID in notRepeat:
                notRepeatSet = set(notRepeat[orthogroupID].keys())
                # these are TE-related, plus annotations denoted as 'True,' meaning no similarity to known gene -- not worth keeping, probably
                if len(notRepeatSet) == 1 and 'True' in notRepeatSet:
                    #print(orthogroupID,notRepeatSet,dataDict[orthogroupID])
                    if orthogroupID not in ogRemoval:
                        ogRemoval[orthogroupID] = dataDict[orthogroupID]
                if len(notRepeatSet) == 1 and 'NA|NA|NA' in notRepeatSet:
                    if orthogroupID not in ogRemoval:
                        ogRemoval[orthogroupID] = dataDict[orthogroupID]
            # these are straight-up TEs only
            else:
                if orthogroupID not in ogRemoval:
                    ogRemoval[orthogroupID] = dataDict[orthogroupID]
                    
    #for orthogroupID in ogRemoval:
    #    print(orthogroupID,ogRemoval[orthogroupID])
    # lines[orthogroupID].append((orthogroupID,orthogroupID,ogProteinCount,populationID,genomeID,geneID,average_entropy,annotation))
    REPEATS = open('repeats_identified_' + orthogroupSetID + '.txt','w')
    REPEATS.write("orthogroupID\trepeatAnnotation\n")
    for orthogroupID in ogRemoval:
        repeatAnnotation = ogRemoval[orthogroupID]
        REPEATS.write("%s\t%s\n" % (orthogroupID,repeatAnnotation))
    
    ### orthogroupID        orthogroupID   hogProteinCount populationID    genomeID        geneID  averageEntropy  annotation
    ### Max theoretical entropy value is 4.3 (no similarity between sequences)
    ### Entropy equal to zero corresponds to no difference between sequences
    OUT = open('all_population_average_entropies_' + orthogroupSetID + '.repeat_filtered.txt','w')
    OUT.write("### orthogroupID\togProteinCount\tpopulationID\tgenomeID\tgeneID\taverageEntropy\tannotation\n")
    OUT.write("### Max theoretical entropy value is 4.3 (no similarity between sequences)\n")
    OUT.write("### Entropy equal to zero corresponds to no difference between sequences\n")
    for orthogroupID in lines:
        if orthogroupID not in ogRemoval:
            # transcriptCoordDict[strand][chromID][transcriptID].sort(key=lambda x: x[1], reverse=True)
            lines[orthogroupID].sort(key=lambda x: x[5], reverse=True)
            for orthogroupID,ogProteinCount,populationID,genomeID,geneID,average_entropy,annotation in lines[orthogroupID]:
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (orthogroupID,ogProteinCount,populationID,genomeID,geneID,average_entropy,annotation))
    

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

usage = "Usage: " + sys.argv[0] + " <chemotype file> <full entropy table> <orthogroup set ID>\n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

chemotypeFile = sys.argv[1]
entropyFile = sys.argv[2]
orthogroupSetID = sys.argv[3]

confirmedTEs = {'OG0000256':'1', 'OG0000395':'1', 'OG0000371':'1', 'OG0000235':'1', 'OG0000196':'1',
                'OG0000040':'1', 'OG0000487':'1', 'OG0000333':'1', 'OG0000137':'1', 'OG0000242':'1',
                'OG0000706':'1', 'OG0000026':'1', 'OG0000336':'1'}

# TE associated descriptions
# Uncharacterized protein K02A2.6-like ## https://www.ebi.ac.uk/interpro/entry/InterPro/IPR034128/ TE-related
repeatDescriptions = {'Reverse transcriptase (RNA-dependent DNA polymerase)':'1', 'ribonuclease H protein':'1', 'Ribonuclease H protein':'1',
                      'zinc-binding in reverse transcriptase':'1', 'GAG-pre-integrase domain':'1',
                      'gag-polypeptide of LTR copia-type':'1', 'transposition, RNA-mediated':'1', 'Retrotransposon gag protein':'1',
                      'ribonuclease H protein At1g65750-like':'1', 'ribonuclease H protein At1g65750':'1', 'reverse transcriptase':'1',
                      'Reverse transcriptase-like':'1', 'Integrase core domain':'1', 'Aspartyl protease':'1', 'Retroviral aspartyl protease':'1',
                      'Integrase core domain containing protein':'1', 'COG2801 Transposase and inactivated derivatives':'1',
                      'PFAM Integrase catalytic region':'1', 'transposon protein':'1', 'Transposase family tnp2':'1', 'transposition':'1',
                      'Retrovirus-related Pol polyprotein from transposon TNT 1-94':'1', 'Transposase-associated domain':'1', 'TNP1/EN/SPM transposase':'1',
                      'Plant transposase (Ptta/En/Spm family)':'1', 'MuDR family transposase':'1', 'MULE transposase domain':'1',
                      'Pfam:UBN2_2':'1', 'Uncharacterized protein K02A2.6-like':'1'}

readEntropyFile(entropyFile,confirmedTEs,repeatDescriptions,orthogroupSetID)
