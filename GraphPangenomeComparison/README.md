### Collect lengths for Minigraph-Cactus variants

<details>
<summary>Download data</summary>
https://resources.michael.salk.edu/root/home.html --> 78csatHaps_minigraphcactus.gaf <br>
</details>

```
python scripts/assessMGCVariants.v3.py 78csatHaps_minigraphcactus.gaf
```

<details>
<summary>Output file names</summary>
mgc_variant_length_info.txt <br>
</details>


### Collect lengths for Syri structural variants
<details>
<summary>Download data</summary>
https://doi.org/10.25452/figshare.plus.25909024.v1 --> *_query_coord.bed <br>
</details>

```
python scripts/collectSyriSVLengths.py DUP_query_coord.bed dups
python scripts/collectSyriSVLengths.py INVs_query_coord.bed invs
python scripts/collectSyriSVLengths.py INVTR_query_coord.bed invtrs
python scripts/collectSyriSVLengths.py TRANS_query_coord.bed trans
```

<details>
<summary>Output file names</summary>
syri_dups_lengths.txt <br>
syri_invs_lengths.txt <br>
syri_invtrs_lengths.txt <br>
syri_trans_lengths.txt <br>
</details>

<details>
<summary>Download data</summary>
https://resources.michael.salk.edu/root/home.html --> *vcf <br>
</details>

```
cat syri_dups_lengths.txt syri_invs_lengths.txt syri_invtrs_lengths.txt syri_trans_lengths.txt > allSyriLengths.txt
for f in ../13csatSexChroms_pggb-oOrient_*vcf; do echo python scripts/filterPGGBVCF.py $f `basename $f .vcf`_filtered.vcf `basename $f .vcf`_filtered.bed; done > filterVCF.sh
cat *_filtered.vcf > allFilteredCsatSexChroms.vcf
python scripts/compareLengths.py allFilteredCsatSexChroms.vcf all16FilteredCsatChroms.vcf allSyriLengths.txt basic_variant_length_info.txt
```

<details>
<summary>Output file names</summary>
graph_vs_syri_length_comparison.png <br>
</details>
