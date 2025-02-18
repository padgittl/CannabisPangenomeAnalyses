### Prep input gene expression file
<details>
<summary>Download data</summary>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/haplotype_biased_expression/quant.tsv <br>
</details>

```
python scripts/createTSVFile.v2.py quant.tsv EH23 0
```

### Find haplotype-biased gene expression in EH23
<details>
<summary>Download data</summary>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/haplotype_biased_expression/EH23_tpm0.0.csv <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/haplotype_biased_expression/EH23_uniprotFunctionalAnnotations.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/haplotype_biased_expression/EH23a_EH23b_ortholog_blast_fractionation.txt <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/haplotype_biased_expression/EH23a_vs_EH23b.collinearity <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/haplotype_biased_expression/EH23a_vs_EH23b.mcscanx.filtered.blastp <br>
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/haplotype_biased_expression/EH23b_vs_EH23a.mcscanx.filtered.blastp <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 EH23a/b GFF3 <br>
https://resources.michael.salk.edu/root/home.html --> Genes V1 EH23a/b mRNA_Table <br>
</details>

```
python scripts/assessHaplotypeBiasedExpression.v3.py EH23_tpm0.0.csv EH23_uniprotFunctionalAnnotations.txt EH23.transcript_summary.tsv EH23.primary_high_confidence.gff3 EH23a_EH23b_ortholog_blast_fractionation.txt EH23a_vs_EH23b.collinearity EH23a_vs_EH23b.mcscanx.filtered.blastp EH23b_vs_EH23a.mcscanx.filtered.blastp EH23a EH23b 0.0
```

<details>
<summary>Output file names</summary>
average_haplotype_biased_expression.txt <br>
average_haplotype_biased_expression_geneIDs_EH23a.txt	<br>
average_haplotype_biased_expression_geneIDs_EH23b.txt	<br>
earlyflower.sr_haplotype_biased_expression.txt	<br>
earlyflower.sr_haplotype_biased_expression_EH23a.txt	<br>
earlyflower.sr_haplotype_biased_expression_EH23b.txt	<br>
foliage.sr_haplotype_biased_expression.txt	<br>
foliage.sr_haplotype_biased_expression_EH23a.txt	<br>
foliage.sr_haplotype_biased_expression_EH23b.txt	<br>
foliage12.sr_haplotype_biased_expression.txt	<br>
foliage12.sr_haplotype_biased_expression_EH23a.txt	<br>
foliage12.sr_haplotype_biased_expression_EH23b.txt	<br>
lateflower.sr_haplotype_biased_expression.txt	<br>
lateflower.sr_haplotype_biased_expression_EH23a.txt	<br>
lateflower.sr_haplotype_biased_expression_EH23b.txt	<br>
roots.sr_haplotype_biased_expression.txt	<br>
roots.sr_haplotype_biased_expression_EH23a.txt	<br>
roots.sr_haplotype_biased_expression_EH23b.txt	<br>
shoottips.sr_haplotype_biased_expression.txt	<br>
shoottips.sr_haplotype_biased_expression_EH23a.txt	<br>
shoottips.sr_haplotype_biased_expression_EH23b.txt	<br>
</details>