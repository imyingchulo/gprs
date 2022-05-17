# About gprs package
This package aims to generate the PRS model from GWAS summary statistics. 
It is designed to deal with the data format based on the GWAS catalog and GeneATLAS database GWAS summary statistics data.

- Understanding the workflow of this package:

1. Filter GWAS summary statistics files (remove duplicate SNPID and select significant SNPs by P-value)
2. Generate bfiles by Plink1.9
3. Do clumping by Plink1.9
4. Generate PRS model by Plink2.0
5. Calculate statistic value of PRS model

## Environment setup

1. Setup virtualenv

```shell
$ python3 -m venv venv
```

2. Activate virtualenv

```shell
$ source ./venv/bin/activate
```

3. Install this package

```shell
$ pip install -e .
```

- Twelve commands in gprs:

1. `geneatlas-filter-data`

2. `gwas-filter-data`

3. `generate-plink-bfiles`

4. `clump`

5. `select-clump-snps`  

6. `build-prs`

7. `combine-prs`

8. `prs-statistics`

9. `combine-prs-stat`

10. `transfer_atcg` (optional)

11. `sub-setpop` (optional)

12. `generate_plink_bfiles_w_individual_info` (optional)

13. `random_draw_samples_from_fam` (optional)

14. `subset_vcf_w_random_sample` (optional)