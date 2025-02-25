import sys, re, os

###############
# SUBROUTINES #
###############

'''
mkdir rawGenomes
for f in peptide/*fa; do echo mkdir rawGenomes/`basename $f .fa`/; done > makeDirFirstLayer.sh 
for f in peptide/*fa; do echo mkdir rawGenomes/`basename $f .fa`/`basename $f .fa`/; done > makeDirSecondLayer.sh 
for f in peptide/*fa; do echo mkdir rawGenomes/`basename $f .fa`/`basename $f .fa`/annotation/; done > makeDirThirdLayer.sh
'''
def readAssemblyIDFile(assemblyIDFile):
    assemblyDict = {}
    with open(assemblyIDFile,'r') as F:
        for line in F:
            assemblyID = line.strip()
            if assemblyID not in assemblyDict:
                assemblyDict[assemblyID] = 1
    return(assemblyDict)
            

def prepareCommands(assemblyDict):
    # protein and gff3 files are expected to be named by https://gitlab.com/salk-tm/snake_tsebra pipeline
    pathToFile = "/path/to/location/of/protein/and/gff3/files/" # expected to be in same dir
    outFileID = 'genespace_dir_commands.sh'
    OUT = open(outFileID,'w')
    OUT.write("mkdir rawGenomes\n")
    OUT_ASSEMBLIES = open("genespace_assemblyIDs.txt",'w')
    for assemblyID in assemblyDict:
        OUT_ASSEMBLIES.write("%s\n" % (assemblyID))
        command1 = "mkdir rawGenomes/"  + assemblyID
        command2 = "mkdir rawGenomes/"  + assemblyID + "/" + assemblyID
        command3 = "mkdir rawGenomes/"  + assemblyID + "/" + assemblyID + "/annotation"
        command4 = "ln -s " + pathToFile + assemblyID + ".primary_high_confidence.gff3 rawGenomes/"  + assemblyID + "/" + assemblyID + "/annotation/" + assemblyID + ".gff3"
        command5 = "ln -s " + pathToFile + assemblyID + ".primary_high_confidence.protein.fasta rawGenomes/" + assemblyID + "/" + assemblyID + "/annotation/" + assemblyID + ".fa"
        OUT.write("%s\n" % (command1))
        OUT.write("%s\n" % (command2))
        OUT.write("%s\n" % (command3))
        OUT.write("%s\n" % (command4))
        OUT.write("%s\n" % (command5))
                    
                
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <assembly ID file> \n"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

assemblyIDFile = sys.argv[1]

assemblyDict = readAssemblyIDFile(assemblyIDFile)
prepareCommands(assemblyDict)
