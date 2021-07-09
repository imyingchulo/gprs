# Genome-wide Polygenic Risk Score (GPRS) 

---

This package is aiming to generate a PRS model from GWAS summary statistics and is designed to deal with GWAS summary statistics from the GWAS catalog and GeneATLAS database.

Download GWAS summary statistics from
[GWAS catalog](https://www.ebi.ac.uk/gwas/) and [GeneATLAS](http://geneatlas.roslin.ed.ac.uk/). 
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


Install both version 1.9 and 2.0 plink.  https://zzz.bwh.harvard.edu/plink/download.shtml


## Usage

### 1. Use python

For example, if you want to run GeneAtlas, create a python file, e.g: `gene_atlas.py`. Paste the following code in the file and change the configuration to your settings. 

```python
from gprs import GeneAtlas


gene_atlas = GeneAtlas( ref='/home1/ylo40816/1000genomes/hg19', data_dir='/home1/ylo40816/Projects/GPRS/data/Gene_ATLAS/selfReported_n_1526' )
gene_atlas.filter_data( snp_id_header='SNP',
                        allele_header='ALLELE',
                        beta_header='NBETA-selfReported_n_1526',
                        se_header='NSE-selfReported_n_1526',
                        pvalue_header='PV-selfReported_n_1526' )

gene_atlas.generate_plink_bfiles(output_name='geneatlas')
gene_atlas.clump(output_name='geneatlas',
                clump_kb='10000',
                clump_p1='1e-3', clump_p2='1e-2')
gene_atlas.select_clump_snps(output_name='geneatlas')
gene_atlas.build_prs( vcf_input= '/home1/ylo40816/1000genomes/hg19',
                      output_name ='geneatlas')
```

### 2. Use Commandline Interface

```shell
$ gprs geneatlas-filter-data --data_dir [str] --result_dir [str] --snp_id_header [str] --allele_header [str] --beta_header --pvalue_header [str] --pvalue [float/scientific notation] --output_name [str]  
$ gprs gwas-filter-data --data_dir [str] --result_dir [str] --snp_id_header [str] --allele_header  [str] --beta_header --pvalue_header [str] --pvalue [float/scientific notation] --output_name [str]  
$ gprs generate-plink-bfiles --ref [str] --output_name [str]
$ gprs clump --data_dir [str] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_p2 [float/scientific notation] --clump_r2 [float] --clump_field [str] --clump_snp_field [str] --output_name [output name]
$ gprs select-clump-snps --output_name [output name]
$ gprs build-prs --vcf_input [str] --columns [int] --plink_modifier [str] --output_name [output name]
```


## Commands in gprs package:

:octocat: Five commands in gprs:

1. `geneatlas-filter-data`

2. `gwas-filter-data`

3. `generate-plink-bfiles`

4. `clump`

5. `select-clump-snps`  

6. `build-prs`

### Result folder
In the first step, you need to indicate the path to creating the result folder.
Five folders will automatically generate under the result folder by script. 

- qc folder: `./result/qc/`
- snplists folder : `./result/snplists/`
- bfile folder: `./result/plink/bfiles`
- clump folder: `./result/plink/clump`
- qc_and_clump_snpslist folder: `./result/plink/qc_and_clump_snpslist`
- prs folder: `./result/plink/prs`

:heavy_exclamation_mark: Users have to indicate reference and result directories every time when using the command interface.

### `gprs geneatlas-filter-data`

Filter GeneAtlas csv file by P-value and unify the data format as following order:
SNPID, ALLELE,  BETA,  StdErr, Pvalue

#### Options:
```
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

#### Result:

This option generate output files in qc and snplists folders:

- `*.QC.csv` (QC files )
- `*.csv` (snplist)

 
### `gprs gwas-filter-data`

Filter GeneAtlas csv file by P-value and unify the data format as following order:
SNPID, ALLELE,  BETA,  StdErr, Pvalue

#### Options:

````
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

#### Result:

This option generate two folders and output files:

- `*.QC.csv` (QC files )
- `*.csv` (snplist)


### `gprs generate-plink-bfiles`

#### Options:
````
  --ref                path to population reference panel  [required]
  --plink_bfiles_dir   plink bfiles output folder, default:"./result/plink/bfiles"
  --snplists_dir       snplists folder, default: "./result/snplists", snplists the name of the file should be chrnb_[output_name].csv i.e. chr1_[output_name].csv;
                       please check "grps gwas_filter_data --help or grps geneatlas_filter_data --help  [required]
  --output_name        output name
  --help               Show this message and exit.
````

#### Result:

This option will generate three files in bfiles folder:

- `*.bim`
- `*.bed`
- `*.fam`

### `gprs clump`

#### Options:
 ````
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

#### Result:

This option will generate one files in clump folder:

- `*.clump`

### `gprs select-clump-snps`

#### Options:
```` 
  --qc_dir                 path to qc folder, default: "./result/qc", the qc files were generated from gwas_filter_data or geneatlas_filter_data options
  --clump_output_dir       path to clump output folder, default: "./result/plink/clump"
  --qc_clump_snplists_dir  path to snpslist (after qc and clumping), default:"./result/plink/qc_and_clump_snpslist"
  --output_name            output name should remain consistent as output_name to plink and filtered data; output format is: [chrnb]_[output_name]_clumped_snplist.csv [required]
  --help                   Show this message and exit.
````

#### Result:
This options will generate one file in qc_and_clump_snpslist folder:
- `*.qc_clump_snpslist.csv`


### `gprs build-prs`
#### Options:
````
  --ref <str>                    path to population reference  [required]
  --prs_output_dir <str>         path to prs model output, default: "./result/plink/prs"
  --qc_clump_snplists_dir <str>  path to snpslist (after qc and clumping), default:"./result/plink/qc_and_clump_snpslist"
  --columns <int>                a column index indicate the [SNPID] [ALLELE] [BETA] position; column nb starts from 1
  --plink_modifier <str>         no-mean-imputation as default in here, get more info by searching plink2.0 modifier
  --output_name <str>            output name should remain consistent as output_name to plink and filtered data [required]
  --help                         Show this message and exit.
````
#### Result:
This options will generate one file in prs folder:
- `*.psam`


