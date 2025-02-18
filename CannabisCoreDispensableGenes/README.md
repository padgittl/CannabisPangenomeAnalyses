# CannabisCoreDispensableGenes

### Get putative orthogroup functions from eggnog annotations
<details>
<summary>Download data</summary>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/Orthogroups.tsv <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/genomeIDs.txt <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 mRNA_Table (for summaryFileList.txt) <br>
</details>

```
ls *transcript_summary.tsv > summaryFileList.txt
python scripts/getOGFunctionFromSummary.v2.py Orthogroups.tsv summaryFileList.txt genomeIDs.txt
```
<details>
<summary>Output file names</summary>
orthogroups_with_annotations.txt <br>
</details>


### Get core and dispensable genes and visualize with stacked bar chart 
<details>
<summary>Download data</summary>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/pangenome_chemotype_population_IDs.csv <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/allGenomeIDs.v2.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/orthogroups_with_annotations.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/Orthogroups_UnassignedGenes.tsv <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/scaffoldedFileIDs.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/publicScaffoldedIDs.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/filteredContigIDs.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/core_dispensable_genes/assembly_focused_core_shell_genes.txt <br>
</details>

```
python scripts/assessOGMembership.v5.py pangenome_chemotype_population_IDs.csv allGenomeIDs.v2.txt orthogroups_with_annotations.txt Orthogroups.tsv Orthogroups_UnassignedGenes.tsv scaffoldedFileIDs.txt publicScaffoldedIDs.txt filteredContigIDs.txt > assembly_focused_core_shell_genes.txt
python scripts/stackedBarChart.py assembly_focused_core_shell_genes.txt allGenomeIDs.v2.txt pangenome_chemotype_population_IDs.csv scaffoldedFileIDs.txt unscaffoldedFileIDs.txt publicFileIDs.txt
```

<details>
<summary>Output file names</summary>
stackedBarChart.png <br>
fullPangenome_pieChart.png <br>
scaffolded_pieChart.png <br>
unscaffolded_pieChart.png <br>
publics_pieChart.png <br>
allGenomes_mj_pieChart.png <br>
allGenomes_feral_pieChart.png <br>
allGenomes_hc_hemp_pieChart.png <br>
allGenomes_F1_pieChart.png <br>
allGenomes_hemp_pieChart.png <br>
allGenomes_asian_hemp_pieChart.png <br>
chromLevelHapResolved_mj_pieChart.png <br>
chromLevelHapResolved_feral_pieChart.png <br>
chromLevelHapResolved_hc_hemp_pieChart.png <br>
chromLevelHapResolved_F1_pieChart.png <br>
chromLevelHapResolved_hemp_pieChart.png <br>
chromLevelHapResolved_asian_hemp_pieChart.png <br>
contigScaffoldLevel_mj_pieChart.png <br>
contigScaffoldLevel_hc_hemp_pieChart.png <br>
contigScaffoldLevel_hemp_pieChart.png <br>
contigScaffoldLevel_asian_hemp_pieChart.png <br>
contigScaffoldLevel_F1_pieChart.png <br>
contigScaffoldLevel_feral_pieChart.png <br>
publics_hc_hemp_pieChart.png <br>
publics_hemp_pieChart.png <br>
publics_asian_hemp_pieChart.png <br>
publics_F1_pieChart.png <br>
publics_mj_pieChart.png <br>
allGenomes_barchart.png <br>
chromLevelHapResolved_barchart.png <br>
contigScaffoldLevel_barchart.png <br>
publics_barchart.png <br>
</details>
