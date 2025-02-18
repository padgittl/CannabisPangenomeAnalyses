#!/bin/python
import sys, os, re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from statistics import mean, median, stdev
import math
import scipy.stats 

###############
# SUBROUTINES #
###############

# https://seqan.readthedocs.io/en/seqan-v2.0.2/Tutorial/VcfIO.html#:~:text=In%20text%20representations%2C%20such%20as%20VCF%2C%201%2Dbased,data%20structures%2C%20SeqAn%20uses%200%2Dbased%20half%2Dopen%20intervals.
# https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
# In text representations, such as VCF, 1-BASED CLOSED intervals are used
# included means CLOSED; excluded means OPEN, this is from interval terminology
# Size = stop - start + 1
# https://github.com/samtools/bcftools/issues/661
# https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants/
# abs(strlen(ALT)-strlen(REF))<=50
# ##INFO=<ID=SVLEN,Number=-1,Type=Integer,Description="Difference in length between REF and ALT alleles">, take longer of two variants if comma separated, but first investigate

# ['AH3Ma.chrX', '10761', '>733598>733601', 'G', 'A', '60.0', '.', 'AC=1;AF=0.5;AN=2;AT=>733598>733599>733601,>733598>733600>733601;NS=2;LV=0', 'GT', '.', '.', '.', '.', '.', '.', '1', '.', '.', '.', '0', '.']
def readVCF(vcfFile,labelGenomeID1,labelGenomeID2):
    diffThreshold = 1
    OUT1 = open(labelGenomeID1, 'w')
    OUT1.write("assemblyID\tchrID\tpos\trefAllele\taltAllele\tsvLen\tlogSVLen\n")
    # print(OUT1)
    # print(labelGenomeID1 + '_filtered.vcf')
    OUT2 = open(labelGenomeID2,'w')
    OUT2.write("assemblyID\tchrID\tstartPos\tstopPos\n")
    with open(vcfFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                # print(line)
                fullChromID = line[0]
                if 'chr' in fullChromID:
                    assemblyID,chrID = fullChromID.split('.')
                    # print(assemblyID,chrID)
                    pos = line[1]
                    startPos = int(pos)
                    refAllele = line[3]
                    altAllele = line[4]
                    if ',' in refAllele:
                        print('unexpected comman in refAllele',refAllele)
                        sys.exit()
                    if ',' in altAllele:
                        altItems = altAllele.split(',')
                        # print(len(altItems),altItems)
                        itemList = []
                        for i in altItems:
                            i = i.strip()
                            # print(i,len(i))
                            itemList.append(len(i))
                            # print(assemblyID,chrID,pos,refAllele,altAllele)
                        # print(max(itemList),itemList)
                        maxValue = max(itemList)
                        diff = abs(len(refAllele)-maxValue)
                        stopPos = startPos + diff
                        if diff >= diffThreshold:
                            logLen = math.log(diff)
                            OUT1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chrID,pos,refAllele,altAllele,diff,logLen))
                            OUT2.write("%s\t%s\t%s\t%s\n" % (assemblyID,chrID,startPos,stopPos))
                    else:
                        diff = abs(len(refAllele)-len(altAllele))
                        stopPos = startPos + diff
                        if diff >= diffThreshold:
                            logLen = math.log(diff)
                            OUT1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chrID,pos,refAllele,altAllele,diff,logLen))
                            OUT2.write("%s\t%s\t%s\t%s\n" % (assemblyID,chrID,startPos,stopPos))
                            
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <pggb vcf file> <genome ID 1> <genome ID 2> \n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

vcfFile = sys.argv[1]
labelGenomeID1 = sys.argv[2]
labelGenomeID2 = sys.argv[3]

readVCF(vcfFile,labelGenomeID1,labelGenomeID2)
