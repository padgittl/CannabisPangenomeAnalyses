# Genespace commands

### prep orthofinder files

```
# Important -- genespace-compatible orthofinder run is expecting fasta files to have specific naming (e.g. AH3Ma.fa, BCMa.fa, GRMa.fa, etc)
tar -zcf myproteins.tar.gz *.fa
aws s3 cp myproteins.tar.gz s3://path/to/orthofinder/dir/
```

### run orthofinder 

<details>
<summary> version info that may be important </summary>
Using snakemake/snakemake:v7.25.0 for Tibanna jobs. <br>
2023-10-03 22:49:03 : Started OrthoFinder version 2.5.5 <br>
</details>

```
nohup snakemake -p -j 1 --use-conda --tibanna \
    --snakefile Snakefile \
    --default-remote-prefix path/to/orthofinder/dir/ \
    --default-resources mem_mb=30000 disk_mb=1000000 \
    --tibanna-config root_ebs_size=16 log_bucket=log_dir \
    --config \
        proteins=myproteins.tar.gz \
        outbase=orthofinder \
        threads=32 \
    > orthofinder.log 2>&1 &
```

```
aws s3 cp s3://path/to/orthofinder/dir/orthofinder.orthofinder.tar.gz .
tar -xvzf orthofinder.orthofinder.tar.gz
```

### prep directories for genespace

```
python scripts/prepareGenespaceDirs.v4.py data/assemblyIDs.txt
chmod +x genespace_dir_commands.sh
./genespace_dir_commands.sh
```

### generate code to orient scaffolds -- specifically helpful for unoriented cannabis genomes, but not strictly necessary

```
python scripts/invertChroms.v4.py data/assemblyIDs.txt data/csat_orientations.tsv AH3Ma
invertTheseChrs = data.frame(genome = c("BCMa","KOMPb","BCMa","GRMa","KOMPa","AH3Mb","BCMa","BCMb","GRMb","KOMPb","BCMa","GRMa","AH3Mb","BCMa","BCMb","KOMPa","KOMPb","AH3Mb","BCMb","GRMb","KOMPa","BCMa","BCMb","GRMb","KOMPa","KOMPb","AH3Mb","BCMa","KOMPb","KOMPa","BCMa","GRMb"), 
chr = c("BCMa.chr1","KOMPb.chr1","BCMa.chr2","GRMa.chr2","KOMPa.chr2","AH3Mb.chr3","BCMa.chr3","BCMb.chr3","GRMb.chr3","KOMPb.chr3","BCMa.chr4","GRMa.chr4","AH3Mb.chr5","BCMa.chr5","BCMb.chr5","KOMPa.chr5","KOMPb.chr5","AH3Mb.chr6","BCMb.chr6","GRMb.chr6","KOMPa.chr6","BCMa.chr7","BCMb.chr7","GRMb.chr7","KOMPa.chr7","KOMPb.chr7","AH3Mb.chr8","BCMa.chr8","KOMPb.chr8","KOMPa.chr9","BCMa.chrX","GRMb.chrX"))
```

### run genespace in RStudio (Version 2022.12.0+353 (2022.12.0+353))

```
# GENESPACE_0.9.3
scripts/genespace.R
```

<details>
<summary>R session info</summary>
> sessionInfo()
R version 4.2.2 (2022-10-31)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Ventura 13.2

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] GENESPACE_0.9.3   data.table_1.14.8


loaded via a namespace (and not attached):
 [1] Rcpp_1.0.10                 lattice_0.20-45             prettyunits_1.1.1           Rsamtools_2.14.0           
 [5] ps_1.7.2                    Biostrings_2.66.0           digest_0.6.31               mime_0.12                  
 [9] R6_2.5.1                    GenomeInfoDb_1.34.7         stats4_4.2.2                zlibbioc_1.44.0            
[13] rlang_1.1.0                 rstudioapi_0.14             miniUI_0.1.1.1              callr_3.7.3                
[17] urlchecker_1.0.1            S4Vectors_0.36.1            R.utils_2.12.2              R.oo_1.25.0                
[21] Matrix_1.5-4.1              devtools_2.4.5              BiocParallel_1.32.5         stringr_1.5.0              
[25] htmlwidgets_1.6.1           igraph_1.3.5                RCurl_1.98-1.9              DelayedArray_0.24.0        
[29] shiny_1.7.4                 compiler_4.2.2              httpuv_1.6.8                rtracklayer_1.58.0         
[33] pkgconfig_2.0.3             BiocGenerics_0.44.0         pkgbuild_1.4.0              htmltools_0.5.4            
[37] SummarizedExperiment_1.28.0 GenomeInfoDbData_1.2.9      matrixStats_0.63.0          IRanges_2.32.0             
[41] codetools_0.2-18            XML_3.99-0.14               crayon_1.5.2                later_1.3.0                
[45] GenomicAlignments_1.34.0    bitops_1.0-7                R.methodsS3_1.8.2           grid_4.2.2                 
[49] xtable_1.8-4                lifecycle_1.0.3             magrittr_2.0.3              cli_3.6.0                  
[53] stringi_1.7.12              cachem_1.0.7                XVector_0.38.0              dbscan_1.1-11              
[57] fs_1.6.1                    promises_1.2.0.1            remotes_2.4.2               ellipsis_0.3.2             
[61] vctrs_0.6.4                 rjson_0.2.21                restfulr_0.0.15             tools_4.2.2                
[65] Biobase_2.58.0              glue_1.6.2                  purrr_1.0.1                 MatrixGenerics_1.10.0      
[69] processx_3.8.0              pkgload_1.3.2               parallel_4.2.2              fastmap_1.1.1              
[73] yaml_2.3.7                  BiocManager_1.30.19         GenomicRanges_1.50.2        sessioninfo_1.2.2          
[77] memoise_2.0.1               profvis_0.3.7               BiocIO_1.8.0                usethis_2.1.6 
</details>
