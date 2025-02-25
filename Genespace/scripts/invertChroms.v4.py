import sys, re, os

###############
# SUBROUTINES #
###############

def readAssemblyIDFile(assemblyIDFile):
    assemblyDict = {}
    with open(assemblyIDFile,'r') as F:
        for line in F:
            assemblyID = line.strip()
            if assemblyID not in assemblyDict:
                assemblyDict[assemblyID] = 1
    return(assemblyDict)
            

def readOrientationFile(assemblyDict,orientationFile,anchorAssemblyID):
    orientationDict = {}
    anchorDict = {}
    with open(orientationFile,'r') as F:
        for line in F:
            if 'Sample' not in line:
                firstCol,assemblyID,chromID,flipBool = line.strip().split('\t')
                if assemblyID in assemblyDict:
                    if chromID not in orientationDict:
                        orientationDict[chromID] = {}
                    if assemblyID not in orientationDict[chromID]:
                        orientationDict[chromID][assemblyID] = flipBool
                    if assemblyID == anchorAssemblyID:
                        if chromID not in anchorDict:
                            anchorDict[chromID] = {}
                        if assemblyID not in anchorDict[chromID]:
                            anchorDict[chromID][assemblyID] = flipBool
    return(orientationDict,anchorDict)


def doInversion(anchorDict,orientationDict):
    assemblyInvertList = []
    chromInvertList = []
    for chromID in anchorDict:
        #print(chromID)
        for anchorAssembly in anchorDict[chromID]:
            flipBool = anchorDict[chromID][anchorAssembly]
            #print(anchorAssembly,chromID,flipBool)
            if chromID in orientationDict:
                #if 'chrY' in chromID:
                for otherAssemblyID in orientationDict[chromID]:
                    otherFlipBool = orientationDict[chromID][otherAssemblyID]
                    if anchorAssembly != otherAssemblyID:
                        if flipBool != otherFlipBool:
                            # print("\t",otherAssemblyID,otherFlipBool)
                            # command = "plot_riparianHits(gpar, genomeIDs = c(" + + ")"
                            quoteOtherAssemblyID = '\"' + otherAssemblyID + '\"'
                            quoteChromID = '\"' + otherAssemblyID + "." + chromID + '\"'
                            # print(quoteOtherAssemblyID)
                            # print(quoteChromID)
                            assemblyInvertList.append(quoteOtherAssemblyID)
                            chromInvertList.append(quoteChromID)
    joinedAssemblyIDs = ','.join(assemblyInvertList)
    joinedChromIDs = ','.join(chromInvertList)
    #print(joinedAssemblyIDs)
    #print(joinedChromIDs)
    command = "invertTheseChrs = data.frame(genome = c(" + joinedAssemblyIDs + "), chr = c(" + joinedChromIDs + "))"  
    print(command)

                
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <assembly ID file> <orientation file> <anchor assembly ID> \n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

assemblyIDFile = sys.argv[1]
orientationFile = sys.argv[2]
anchorAssemblyID = sys.argv[3]

assemblyDict = readAssemblyIDFile(assemblyIDFile)

orientationDict,anchorDict = readOrientationFile(assemblyDict,orientationFile,anchorAssemblyID)
doInversion(anchorDict,orientationDict)
