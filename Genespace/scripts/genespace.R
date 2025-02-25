if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE", upgrade = F)
library(GENESPACE)
sessionInfo() 

# 1
runwd <- file.path("/current/working/directory/")
setwd("/current/working/directory/")

file.path(runwd, "rawGenomes")
runwd

assemblyIDs <- read.delim("assemblyIDs.txt", header=FALSE)
csIDs <- paste(assemblyIDs, collapse = ",")
length(list(csIDs))
cat(csIDs)
# c("AH3Ma", "BCMa", "GRMb", "KOMPa", "AH3Mb", "BCMb", "GRMa", "KOMPb"),

# 2
?init_genespace
gpar <- init_genespace(
  genomeIDs = c("AH3Ma", "BCMa", "GRMb", "KOMPa", "AH3Mb", "BCMb", "GRMa", "KOMPb"),
  speciesIDs = c("AH3Ma", "BCMa", "GRMb", "KOMPa", "AH3Mb", "BCMb", "GRMa", "KOMPb"),
  versionIDs = c("AH3Ma", "BCMa", "GRMb", "KOMPa", "AH3Mb", "BCMb", "GRMa", "KOMPb"),
  outgroup = NULL,
  ploidy = 1,
  wd = runwd,
  overwrite = F, 
  verbose = T,
  nCores = 1,
  minPepLen = 50,
  orthofinderMethod = "fast",
  path2orthofinder = "~/path/to/orthofinder/bin/orthofinder",
  path2diamond = "~/path/to/orthofinder/bin/diamond",
  diamondMode = "fast",
  orthofinderInBlk = FALSE, 
  gffString = "gff3",
  pepString = "fa",
  path2mcscanx = "~/path/to/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))
gpar

# 3 
parse_annotations(
  gsParam = gpar,
  gffEntryType = "mRNA",
  gffIdColumn = "ID",
  gffStripText = "ID=",
  headerEntryIndex = 1,
  headerSep = "\t",
  headerStripText = "ID=")

# check to make sure genespace can find orthofinder results
find_orthofinderResults(gsParam = gpar, onlyCheckRun = F)

# run synteny 
gpar <- synteny(gsParam = gpar)

# 4 v1
plot_riparianHits(gpar)
plot_riparianHits(gpar, blackBg = FALSE)

# 4 v2
# get the chromosome lengths for AH3Ma (from gene gff3 file [tsebra annotation])
##sequence-region AH3Ma.chr1 1 65671371
##sequence-region AH3Ma.chr2 1 76481500
##sequence-region AH3Ma.chr3 1 81995129
##sequence-region AH3Ma.chr4 1 81631532
##sequence-region AH3Ma.chr5 1 76536173
##sequence-region AH3Ma.chr6 1 76260137
##sequence-region AH3Ma.chr7 1 66405468
##sequence-region AH3Ma.chr8 1 54135332
##sequence-region AH3Ma.chr9 1 66045532
##sequence-region AH3Ma.chrX 1 84231629

#### This reproduces Supplemental Figure 17a (bioRxiv Supp. Fig. 13 https://www.biorxiv.org/content/10.1101/2024.05.21.595196v1.full.pdf)
regs <- data.frame(
  genome = c("AH3Ma", "AH3Ma", "AH3Ma", "AH3Ma", "AH3Ma", 
             "AH3Ma", "AH3Ma", "AH3Ma", "AH3Ma", "AH3Ma"),
  chr = c("AH3Ma.chr1", "AH3Ma.chr2", "AH3Ma.chr3", "AH3Ma.chr4", "AH3Ma.chr5", 
          "AH3Ma.chr6", "AH3Ma.chr7", "AH3Ma.chr8", "AH3Ma.chr9", "AH3Ma.chrX"),
  start = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  end = c(65671371,76481500,81995129,81631532,76536173,76260137,66405468,54135332,66045532,84231629),
  cols = c("blue", "violet", "darkmagenta", "black", "yellow", "orange",
                 "yellowgreen", "aquamarine", "royalblue2",
                 "mediumseagreen"))
plot_riparianHits(gpar, blackBg = FALSE, refGenome = "AH3Ma",
                  onlyTheseRegions = regs, useOrder = TRUE, 
                  invertTheseChrs = data.frame(genome = c("KOMPb","BCMb","AH3Mb","BCMa","KOMPb","BCMa","GRMa","KOMPa","AH3Mb","BCMa","BCMb","GRMb","KOMPb","BCMa","GRMa","AH3Mb","BCMa","BCMb","KOMPa","KOMPb","AH3Mb","BCMb","GRMb","KOMPa","BCMa","BCMb","GRMb","KOMPa","KOMPb","AH3Mb","BCMa","KOMPb","KOMPa","BCMa","GRMb"), 
                                               chr = c("KOMPb.chrY","BCMb.chrY","AH3Mb.chrY","BCMa.chr1","KOMPb.chr1","BCMa.chr2","GRMa.chr2","KOMPa.chr2","AH3Mb.chr3","BCMa.chr3","BCMb.chr3","GRMb.chr3","KOMPb.chr3","BCMa.chr4","GRMa.chr4","AH3Mb.chr5","BCMa.chr5","BCMb.chr5","KOMPa.chr5","KOMPb.chr5","AH3Mb.chr6","BCMb.chr6","GRMb.chr6","KOMPa.chr6","BCMa.chr7","BCMb.chr7","GRMb.chr7","KOMPa.chr7","KOMPb.chr7","AH3Mb.chr8","BCMa.chr8","KOMPb.chr8","KOMPa.chr9","BCMa.chrX","GRMb.chrX")))


# 5
pg <- pangenome(gpar)
pg




