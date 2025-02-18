# CannabisSequenceEntropy

<details>
<summary>Download data</summary>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/Orthogroups.tsv <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/pangenome_chemotype_population_IDs.csv <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 mRNA_Table (for summaryFileList.txt) <br>
</details>

### Calculate entropy values for orthofinder-based multiple sequence alignments
```
nohup python scripts/analyzeOGs.py pangenome_chemotype_population_IDs.csv Orthogroups.tsv msaFileList.txt summaryFileList.txt full_run > full_run.log 2>&1 &
```
<details>
<summary>Output file names</summary>
all_population_average_entropies.txt <br>
</details>

### Remove repeat-associated gene groups
```
python scripts/filterGroups.v3.py pangenome_chemotype_population_IDs.csv all_population_average_entropies.txt all_orthogroups
```
<details>
<summary>Output file names</summary>
all_population_average_entropies_all_orthogroups.repeat_filtered.txt <br>
</details>

### Create joint histogram to visualize population-level entropy comparisons
```
python scripts/jointHistogram.v4.py all_population_average_entropies.txt all_population_average_entropies_all_orthogroups.repeat_filtered.txt 5
```
<details>
<summary>Output information</summary>
pop1    pop2    statistic       pvalue <br>
asian_hemp feral 0.8375639347856714 0.0 <br>
asian_hemp hc_hemp 0.8395177728884496 0.0 <br>
asian_hemp hemp 0.8274982870253889 0.0 <br>
asian_hemp mj 0.8392839827683827 0.0 <br>
F1 asian_hemp 0.8390746784938936 0.0 <br>
F1 feral 0.8657673167634514 0.0 <br>
F1 hc_hemp 0.9024173330520182 0.0 <br>
F1 hemp 0.8618989236995579 0.0 <br>
F1 mj 0.9148663951998351 0.0 <br>
feral hc_hemp 0.8659468310917966 0.0 <br>
feral hemp 0.8739148670378099 0.0 <br>
feral mj 0.8609754305102137 0.0 <br>
hc_hemp hemp 0.8547936333921432 0.0 <br>
hc_hemp mj 0.9109642807961587 0.0 <br>
hemp mj 0.8553464635053353 0.0 <br>
</details>

### Create visualization of entropy values for each position in multiple sequence alignment
```
python scripts/evaluateMSA.v2.py pangenome_chemotype_population_IDs.csv Orthogroups.tsv OG0003894.fa summaryFileList.txt all_population_average_entropies.txt OG0003894 PRR
python scripts/evaluateMSA.v2.py pangenome_chemotype_population_IDs.csv Orthogroups.tsv OG0004389.fa summaryFileList.txt all_population_average_entropies.txt OG0004389 PRR
python scripts/evaluateMSA.v2.py pangenome_chemotype_population_IDs.csv Orthogroups.tsv OG0016523.fa summaryFileList.txt all_population_average_entropies.txt OG0016523 GI
python scripts/evaluateMSA.v2.py pangenome_chemotype_population_IDs.csv Orthogroups.tsv OG0006134.fa summaryFileList.txt all_population_average_entropies.txt OG0006134 FT
python scripts/evaluateMSA.v2.py pangenome_chemotype_population_IDs.csv Orthogroups.tsv OG0000392.fa summaryFileList.txt all_population_average_entropies.txt OG0000392 sesquiterpene
```
<details>
<summary>Output file names</summary>
OG0003894_PRR.png <br>
OG0004389_PRR.png <br> 
OG0016523_GI.png <br> 
OG0006134_FT.png <br> 
OG0000392_sesquiterpene.png <br> 
</details>
