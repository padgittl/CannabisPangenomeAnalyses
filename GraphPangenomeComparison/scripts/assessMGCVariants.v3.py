#!/bin/python
import sys, os, re, math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

###############
# SUBROUTINES #
###############

# coordinate system for gaf file
# https://d-nb.info/1222502348/34
# https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/, closed and open
# gaf file is "0 start, half open" according to UCSC, which is format for tables -- "start included, end excluded"
# included means CLOSED; excluded means OPEN, this is from interval terminology
# this is the same as gaf specs, which say, "Query start coordinate (0-based; closed)"; "Query end coordinate (0-based; open)"
# size = stop - start
def readGAF(gafFile):
    graphDict = {}
    with open(gafFile,'r') as F:
        for line in F:
            line = line.strip().split('\t')
            queryID = line[0]
            # id=AH3Ma|AH3Ma.chr1
            getRefGenomeID = re.search(r'id=(.+)\|(.+)',queryID)
            refGenomeID = getRefGenomeID.group(1)
            # print(refGenomeID)
            chromosomeInfo = getRefGenomeID.group(2)
            # AH3Ma.chr1
            gID,chrID = chromosomeInfo.split('.')
            # print(chromosomeInfo)
            # print(chrID)
            queryLen = line[1]
            queryStart = line[2]
            queryStop = line[3]
            varLabel = queryStart + '-' + queryStop
            strand = line[4]
            if 'chr' in queryID:
                if refGenomeID not in graphDict:
                    graphDict[refGenomeID] = {}
                if chrID not in graphDict[refGenomeID]:
                    graphDict[refGenomeID][chrID] = {}
                if varLabel not in graphDict[refGenomeID][chrID]:
                    graphDict[refGenomeID][chrID][varLabel] = []
                delimiters = "<id=", ">id="
                results = re.split(f"[{delimiters}]", line[5])
                for i in results:
                    if i:
                        items = i.split(':')
                        if len(items) == 2:
                            #print(refGenomeID,varLabel,items,results)
                            getInfo = re.search(r'(.+)\|(.+):(\d+)-(\d+)',i)
                            #print(i)
                            genomeID = getInfo.group(1)
                            fullChromID = getInfo.group(2)
                            varStart = getInfo.group(3)
                            varStop = getInfo.group(4)
                            varStart = int(varStart)
                            varStop = int(varStop)
                            varLen = varStop - varStart
                            graphDict[refGenomeID][chrID][varLabel].append((genomeID,fullChromID,varStart,varStop,varLen))
                            # print(genomeID,fullChromID,varStart,varStop)
                        else:
                            if len(items) == 1:
                                if len(line) == 19:
                                    genomeID,fullChromID = items[0].split('|')
                                    varStart = line[7]
                                    varStop = line[8]
                                    varStart = int(varStart)
                                    varStop = int(varStop)
                                    varLen = varStop - varStart
                                    graphDict[refGenomeID][chrID][varLabel].append((genomeID,fullChromID,varStart,varStop,varLen))
                                else:
                                    print('other parsing problem',refGenomeID,varLabel,len(line),line)
                                    sys.exit()
                            else:
                                print(refGenomeID,varLabel,len(items),items)
                                sys.exit()

    OUT = open('basic_variant_length_info.txt','w')
    print("basic_variant_length_info.txt")
    OUT.write("refGenomeID\tchrID\tqueryStart\tqueryStop\tvarLen\tlogLen\n")
    for refGenomeID in graphDict:
        for chrID in graphDict[refGenomeID]:
            # queryStart + '-' + queryStop
            for varLabel in graphDict[refGenomeID][chrID]:
                queryStart,queryStop = varLabel.split('-')
                queryStart = int(queryStart)
                queryStop = int(queryStop)
                varLen = queryStop - queryStart
                logLen = math.log(varLen)
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (refGenomeID,chrID,queryStart,queryStop,varLen,logLen))

                    
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <gaf file> \n"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

gafFile = sys.argv[1]

readGAF(gafFile)
