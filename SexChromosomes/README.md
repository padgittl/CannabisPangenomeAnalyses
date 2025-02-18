# Cannabis Sex Chromosome Analysis 
## This analysis produces the results shown in Fig.2 and Extended Data Fig. 2


<details>
<summary>Download data</summary>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/AH3Ma_vs_AH3Mb.collinearity  <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/AH3M_tpm0.0.csv <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/AH3Ma_vs_AH3Mb.mcscanx.filtered.blastp <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/AH3Mb_vs_AH3Ma.mcscanx.filtered.blastp <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/AH3M_uniprotFunctionalAnnotations.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/ordered_syntenic_and_non_syntenic_genes.txt <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 AH3Ma/b GFF3 <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 AH3Ma/b mRNA_Table <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 BCMa/b mRNA_Table <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 EH23a/b mRNA_Table <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 GRMa/b mRNA_Table <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 KOMPa/b mRNA_Table <br>
</details>

### Prep commands
```
cat AH3Ma.primary_high_confidence.gff3 AH3Mb.primary_high_confidence.gff3 > AH3M.primary_high_confidence.gff3
cat AH3Ma.transcript_summary.tsv AH3Mb.transcript_summary.tsv > AH3M.transcript_summary.tsv
cat AH3Ma.transcript_summary.tsv AH3Mb.transcript_summary.tsv BCMa.transcript_summary.tsv BCMb.transcript_summary.tsv EH23a.transcript_summary.tsv EH23b.transcript_summary.tsv GRMa.transcript_summary.tsv GRMb.transcript_summary.tsv KOMPa.transcript_summary.tsv KOMPb.transcript_summary.tsv > compiled_transcript_summary.tsv
```

### Calculate average expression for male and female plants
```
python scripts/getAverageExpression.py AH3Ma_vs_AH3Mb.collinearity AH3M_tpm0.0.csv AH3Ma_vs_AH3Mb.mcscanx.filtered.blastp AH3Mb_vs_AH3Ma.mcscanx.filtered.blastp AH3M_uniprotFunctionalAnnotations.txt AH3M.transcript_summary.tsv AH3M.primary_high_confidence.gff3 AH3Ma AH3Mb 0.0
```
<details>
<summary>Output file names</summary>
female_vs_male_all_average_exp.txt
</details>


### Identify syntenic genes from genespace
<details>
<summary>Download data</summary>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/assembliesToTest.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/gffWithOgs.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/all_annotations.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/terpeneAnnotations.v6.tsv <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/sex_chromosomes/female_vs_male_all_average_exp.txt <br>
</details>

```
python scripts/getSyntenicGenes.py gffWithOgs.txt compiled_transcript_summary.tsv assembliesToTest.txt terpeneAnnotations.v6.tsv all_annotations.txt female_vs_male_all_average_exp.txt
```
<details>
<summary>Output file names</summary>
present_in_chrX_and_chrY.txt <br>
present_in_only_chrX.txt <br>
present_in_only_chrY.txt <br>
present_in_only_autosomes.txt <br>
terpene_synthase_related_syntenic_genes.txt <br>
ordered_syntenic_and_non_syntenic_genes.txt <br>
genome_total_synteny_counts.txt <br>
</details>


### Identify biased gene expression in male vs female plants
```
python scripts/collectSyntenicExpressedGenes.py AH3Ma_vs_AH3Mb.collinearity AH3M_tpm0.0.csv AH3Ma_vs_AH3Mb.mcscanx.filtered.blastp AH3Mb_vs_AH3Ma.mcscanx.filtered.blastp AH3M_uniprotFunctionalAnnotations.txt AH3M.transcript_summary.tsv AH3M.primary_high_confidence.gff3 ordered_syntenic_and_non_syntenic_genes.txt AH3Ma AH3Mb 0.0 1.0
```

<details>
<summary>Output file names</summary>
flower_male_only_expression.txt <br>
flower_female_only_expression.txt  <br>
flower_male_higher_exp_min5TPM_diffs.txt  <br>
flower_female_higher_exp_min5TPM_diffs.txt  <br>
flower_balanced_male_female_expression.txt  <br>
leaf_male_only_expression.txt  <br>
leaf_female_only_expression.txt  <br>
leaf_male_higher_exp_min5TPM_diffs.txt <br>  
leaf_female_higher_exp_min5TPM_diffs.txt  <br>
leaf_balanced_male_female_expression.txt  <br>
</details>


### Create bar chart to visualize gene counts
```
python scripts/plotCounts.py flower_balanced_male_female_expression.txt leaf_balanced_male_female_expression.txt flower_male_only_expression.txt leaf_male_only_expression.txt flower_female_only_expression.txt leaf_female_only_expression.txt flower_male_higher_exp_min5TPM_diffs.txt leaf_male_higher_exp_min5TPM_diffs.txt flower_female_higher_exp_min5TPM_diffs.txt leaf_female_higher_exp_min5TPM_diffs.txt
```
<details>
<summary>Output file names</summary>
male_female_barchart.png <br>
syntenic_percent.txt <br>
</details>


### Concatenate biased expression files for downstream analyses
```
cat flower_male_only_expression.txt flower_male_higher_exp_min5TPM_diffs.txt > flower_male_biased_expression.txt
cat flower_female_only_expression.txt flower_female_higher_exp_min5TPM_diffs.txt > flower_female_biased_expression.txt
cat leaf_male_only_expression.txt leaf_male_higher_exp_min5TPM_diffs.txt > leaf_male_biased_expression.txt
cat leaf_female_only_expression.txt leaf_female_higher_exp_min5TPM_diffs.txt > leaf_female_biased_expression.txt
```


### Create joint histogram and scatter plot to show differences in gene expression in the PAR and SDR
```
python scripts/geneExpressionJointHist.py flower_male_biased_expression.txt flower_male male
```
<details>
<summary>Output file names</summary>
AH3Mb_chrY_flower_male_pos_hist.png <br>
AH3Ma_chrX_flower_male_pos_hist.png <br>
</details>


