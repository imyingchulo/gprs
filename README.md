# Genome-wide Polygenic Risk Score (GPRS) 

---

This package is aiming to generate a PRS model from GWAS summary statistics and is designed to deal with GWAS summary statistics from the GWAS catalog and GeneATLAS database.

Download GWAS summary statistics from
GWAS catalog: https://www.ebi.ac.uk/gwas/
GeneATLAS: http://geneatlas.roslin.ed.ac.uk/
(The GWAS summary statistics should have contains the information: SNPID, ALLELE, BETA, P-value, and StdErr)

:octocat: Understanding the workflow of this package:
1. Filter GWAS summary statistics files (remove duplicate SNPID and select significant SNPs by P-value)
2. Generate bfiles by Plink1.9
3. Do clumping by Plink1.9
4. Generate PRS model by Plink2.0


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

## Additional requirements
Plink version 1.9 and 2.0
install plink: https://zzz.bwh.harvard.edu/plink/download.shtml


## Execute script

```shell
$gprs geneatlas-filter-data --data_dir [str] --result_dir [str] --snp_id_header [str] --allele_header [str] --beta_header --pvalue_header [str] --pvalue [float/scientific notation] --output_name [str]  
$gprs gwas-filter-data --data_dir [str] --result_dir [str] --snp_id_header [str] --allele_header  [str] --beta_header --pvalue_header [str] --pvalue [float/scientific notation] --output_name [str]  
$gprs generate-plink-bfiles --ref [str] --plink_bfiles_dir [str] --snplists_dir [str] --output_name [str]
$gprs clump --data_dir [str] --plink_bfiles_dir [str] --clump_output_dir [str] --qc_dir [str] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_p2 [float/scientific notation] --clump_r2 [float] --clump_field [str] --clump_snp_field [str] --output_name [output name]
$gprs select-clump-snps --qc_dir [str] --clump_output_dir [str] --qc_clump_snplists_dir [str] --output_name [output name]
$gprs build-prs --vcf_input [str] --prs_output_dir [str] --qc_clump_snplists_dir [str] --columns [int] --plink_modifier [str] --output_name [output name]
```


## Commands in gprs package:

:octocat: Five commands in gprs:

1. `geneatlas-filter-data`

2. `gwas-filter-data`

3. `generate-plink-bfiles`

4. `clump`

5. `select-clump-snps`  

6. `build-prs`


### `geneatlas-filter-data`

Filter GeneAtlas csv file by P-value and unify the data format as following order:
SNPID, ALLELE,  BETA,  StdErr, Pvalue

This option generate two folders and output files:

./result/qc/

./result/snplists/

QC files (.QC.csv)

snplist (.csv)
```
Options:
  --data_dir                The directory of GeneAtlas csv files (all 1-24 chr) [required]
  --result_dir              path to output folder; default:[./result]
  --snp_id_header           SNP ID column name in GeneAtlas original file  [required]
  --allele_header           ALLELE column name in GeneAtlas original file  [required]
  --beta_header             BETA column name in GeneAtlas original file [required]
  --se_header               StdErr column name in GeneAtlas original file  [required]
  --pvalue_header           P-value column name in GeneAtlas original file  [required]
  --output_name             output name; default: geneatlas; the output file will be [chrnb]_[output_name].csv
  --pvalue                  P-value threshold
  --help                    Show this message and exit.
````
 
### `gwas-filter-data`
Filter GeneAtlas csv file by P-value and unify the data format as following order:
SNPID, ALLELE,  BETA,  StdErr, Pvalue

This option generate two folders and output files:

./result/qc/

./result/snplists/

QC files (.QC.csv)

snplist (.csv)
````
Options:
  --data_dir                 path to GWAS catalog summary statistic csv file (all 1-24 chr)  [required]
  --result_dir               path to output folder; default:[./result]
  --snp_id_header            SNP ID column name in GWAS catalog original file  [required]
  --allele_header            ALLELE column name in GWAS catalog original file  [required]
  --beta_header              BETA column name in GWAS catalog original file  [required]
  --se_header                StdErr column name in GWAS catalog original file  [required]
  --pvalue_header            P-value column name in GWAS catalog original file  [required]
  --output_name              output name; default: gwas ; the output file will be named as [chrnb]_[output_name].csv
  --pvalue                   P-value threshold for filtering SNPs
  --help                     Show this message and exit.
````

### `generate-plink-bfiles`
This option will generate ./result/plink/bfiles folder and .bim .bed and .fam files

````
Options:
  --ref                path to population reference panel  [required]
  --plink_bfiles_dir   plink bfiles output folder, default:"./result/plink/bfiles"
  --snplists_dir       snplists folder, default: "./result/snplists", snplists the name of the file should be chrnb_[output_name].csv i.e. chr1_[output_name].csv;
                       please check "grps gwas_filter_data --help or grps geneatlas_filter_data --help  [required]
  --output_name        output name
  --help               Show this message and exit.
````

### `clump`
This options will generate ./result/plink/clump folder and .clump file
 ````
Options:
  --data_dir                 path to GWAS catalog/GeneAtlas .csv file  [required]
  --plink_bfiles_dir         plink bfiles output folder, default:"./result/plink/bfiles"
  --clump_output_dir         path to clump output folder, default:"./result/plink/clump"
  --qc_dir                   path to qc folder, default: "./result/qc", the qc files were generated from gwas_filter_data or geneatlas_filter_data options
  --clump_kb                 distance(kb) parameter for clumping [required]
  --clump_p1                 first set of P-value for clumping [required]
  --clump_p2                 should equals to p1 reduce the snps [required]
  --output_name              output name should remain consistent as output_name to plink and filtered data; output format: [chrnb]_[output_name]_clumped_snplist.csv [required]
  --clump_r2                 r2 value for clumping, default = 0.1
  --clump_field              P-value column name, default = Pvalue
  --clump_snp_field          SNP ID column name, default = SNPID
  --help                     Show this message and exit.
````

### `select-clump-snps`
This options will generate ./result/plink/clump folder and .clump file
```` 
Options:
  --qc_dir             path to qc folder, default: "./result/qc", the qc files were generated from gwas_filter_data or geneatlas_filter_data options
  --clump_output_dir   path to clump output folder, default: "./result/plink/clump"
  --output_name        output name should remain consistent as output_name to plink and filtered data; output format is: [chrnb]_[output_name]_clumped_snplist.csv [required]
  --help               Show this message and exit.
````

### `build-prs`
This options will generate ./result/plink/prs folder and .score file
````
Options:
  --ref <str>                    path to population reference  [required]
  --prs_output_dir <str>         path to prs model output, default: "./result/plink/prs"
  --qc_clump_snplists_dir <str>  path to snpslist (after qc and clumping), default:"./result/plink/qc_and_clump_snpslist"
  --columns <int>                a column index indicate the [SNPID] [ALLELE] [BETA] position; column nb starts from 1
  --plink_modifier <str>         no-mean-imputation as default in here, get more info by searching plink2.0 modifier
  --output_name <str>            output name should remain consistent as output_name to plink and filtered data [required]
  --help                         Show this message and exit.
````


