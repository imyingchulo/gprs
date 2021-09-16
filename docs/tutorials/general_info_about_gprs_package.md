# About gprs package:
This package aims to generate a PRS model from GWAS summary statistics data-set. 
It is designed to use GWAS summary statistics from the GWAS catalog and GeneATLAS database.

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

- Ten commands in gprs:

1. `geneatlas-filter-data`

2. `gwas-filter-data`

3. `generate-plink-bfiles`

4. `clump`

5. `select-clump-snps`  

6. `build-prs`

7. `transfer_atcg` (optional)

8. `combine-prs`

9. `prs-statistics`

10. `combine-prs-stat`