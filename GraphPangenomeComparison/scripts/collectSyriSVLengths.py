#!/bin/python
import sys, os, re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean, median, stdev
import math

###############
# SUBROUTINES #
###############

#Chrom   Start   End     Name
#chr1    13839866        13852176        AH3Ma
#chr6    38673417        38685697        DUP     AH3Ma
#refGenomeID     chrID   queryStart      queryStop       varLen  logLen
#bed files are 0-based
def readSVFile(svFile,svType):
    OUT = open('syri_' + svType + '_lengths.txt','w')
    OUT.write("genome\tchrID\tstart\tend\tsvLength\tlogLen\n")
    print('syri_' + svType + '_lengths.txt')
    with open(svFile,'r') as F:
        for line in F:
            if 'Chrom' not in line:
                chrID,start,end,svID,genome = line.strip().split('\t')
                start = int(start)
                end = int(end)
                svLength = end - start
                logLen = math.log(svLength)
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (genome,chrID,start,end,svLength,logLen))

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <SV file> <SV type> \n"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

svFile = sys.argv[1]
svType = sys.argv[2]

readSVFile(svFile,svType)
