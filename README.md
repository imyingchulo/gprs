# Genome-wide Polygenic Risk Score (GPRS) 

[![Documentation Status](https://readthedocs.org/projects/gprs/badge/?version=latest)](https://gprs.readthedocs.io/en/latest/?badge=latest)

---

This package aims to generate the PRS model from GWAS summary statistics. 
It is designed to deal with GWAS summary statistics from the GWAS catalog and GeneATLAS database or similar format datasets.

:octocat: Understanding the workflow of this package:

1. Filter GWAS summary statistics files (unify the data format. optional: remove duplicate SNPID and select significant SNPs by P-value)
2. Generate bfiles by Plink1.9
3. Do clumping and thresholding by Plink1.9
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
$ pip install -r requirements.txt
$ pip install -e .
```

If you are USC users
```shell
$ pip install . --user 
```

## Additional requirements
Install both version 1.9 and 2.0 plink.  https://zzz.bwh.harvard.edu/plink/download.shtml

If you are USC users:
Please load the plink modules by using:
```shell
$ module load plink2
```
```shell
$ module load gcc/11.2.0
$ module load plink/1.9-beta6.24
```

## Prepare dataset
Download GWAS summary statistics from
[GWAS catalog](https://www.ebi.ac.uk/gwas/) and [GeneATLAS](http://geneatlas.roslin.ed.ac.uk/). 
(The GWAS summary statistics should have contains the information: SNPID, ALLELE, BETA, P-value, and StdErr)


## Usage guidance
### Before starting:
If your GWAS summary statistics data are zipped, please unzip your .gz file first.

## Two ways to run the package: Python or Commandline Interface
### 1. Use Python

For example, if you want to run GeneAtlas, create a python file, e.g: `gene_atlas_model.py`. 
Paste the following code in the file and change the configuration to your settings. 

```python
from gprs.gene_atlas_model import GeneAtlasModel

if __name__ == '__main__':
    geneatlas = GeneAtlasModel( ref='/home/1000genomes/GRCh38',
                    data_dir='/home/Projects/GPRS/data/2014_GWAS_Height' )

    geneatlas.filter_data( snp_id_header='MarkerName',
                            allele_header='Allele1',
                            beta_header='b',
                            se_header ='SE',
                            pvalue_header='p',
                            output_name='2014height')

    geneatlas.generate_plink_bfiles(snplist_name='2014height', output_name='2014height',extra_commands="--vcf-half-call r" ,symbol='.genotypes')

    geneatlas.clump(output_name='2014height',plink_bfile_name='2014height',
                    qc_file_name='2014height',clump_kb='250',
                    clump_p1='0.02',clump_p2='0.02', clump_r2='0.02')

    geneatlas.select_clump_snps(output_name='2014height',
                                clump_file_name='2014height',
                                qc_file_name='2014height',
                                clump_kb='250',
                                clump_p1='0.02',
                                clump_r2='0.5',clumpfolder_name='2014height')

    geneatlas.build_prs( vcf_input= 'home/1000genomes/hg19',
                         output_name ='2014height', memory='1000',clump_kb='250',
                         clump_p1='0.02', clump_r2='0.02', qc_clump_snplist_foldername='2014height')

    geneatlas.combine_prs(filename="2014height",
                          clump_r2="0.5",clump_kb="250",clump_p1="0.02")

    geneatlas.prs_statistics(output_name='2014height', score_file = "home/GPRS/tmp/2014height_250_0.02_0.1.sscore",
                             pheno_file = "prs/2014height_pheno.csv",
                             r_command='/bin/Rscript',
                             prs_stats_R="/GPRS/gprs/prs_stats_quantitative_phenotype.R", 
                             data_set_name="2014height",
                             clump_kb='250',
                             clump_p1='0.02',
                             clump_r2='0.1')
    
    geneatlas.subset_pop(input_data="home/1000genomes/hg19/integrated_call_samples_v3.20130502.ALL.panel",
                         column_name="super_pop",pop_info = "EUR", output_name = "2014height")

    geneatlas.generate_plink_bfiles_w_individual_info(popfile_name="2014height", output_name="2014height",bfile_name="2014height")
    
    geneatlas.subset_vcf_w_random_sample(fam_dir="/home/result/plink/bfiles", fam_filename="2014height", samplesize=50, vcf_input="1000genomes/hg19", symbol=".")
    
    geneatlas.random_draw_samples_from_fam(fam_dir="/home/result/plink/bfiles",fam_filename="2014height", samplesize= "500", tag="LD_reference") 

```
 
### 2. Use Commandline Interface

```shell
$ gprs geneatlas-filter-data --ref [str] --data_dir [str] --result_dir [str] --snp_id_header [str] --allele_header [str] --beta_header [str] --se_header [str] --pvalue_header [str] --pvalue [float/scientific notation] --output_name [str]  
$ gprs gwas-filter-data --ref [str] --data_dir [str] --result_dir [str] --snp_id_header [str] --allele_header  [str] --beta_header [str] --se_header [str] --pvalue_header [str] --pvalue [float/scientific notation] --output_name [str]  
$ gprs generate-plink-bfiles --ref [str] --snplist_name [str] --output_name [str] --symbol [str] --extra_commands [str] 
$ gprs clump --plink_bfile_name [str] --output_name [str] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_p2 [float/scientific notation] --clump_r2 [float] --clump_field [str] --qc_file_name [str] --clump_snp_field [str]   
$ gprs select-clump-snps --qc_file_name [str] --clump_file_name [str] --output_name [output name] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_r2 [float] --clumpfolder_name [str]
$ gprs build-prs --vcf_input [str] --output_name [str] --qc_clump_snplist_foldername [str] --memory [int] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_r2 [float] --symbol [str/int] --columns [int] --plink_modifier [str]  
$ gprs combine-prs --filename [str] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_r2 [float]
$ gprs prs-statistics --score_file [str] --pheno_file [str] --output_name [str] --data_set_name [str] --prs_stats_R [str] --r_command [str] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_r2 [float]
$ gprs combine-prs-stat --data_set_name [str] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_r2 [float]

```

### optional function
If alleles are a, t, c, g instead of capital A, T, C, G it might affect the further analysis. 
```shell
$ gprs transfer_atcg --qc_file_name [str]
$ gprs subset_pop --input_data [str] --column_name [str] --pop_info [str] --output_name [str]
$ gprs generate-plink-bfiles-w-individual-info --popfile_name [str] --bfile_name [str] --output_name [str]
$ gprs subset-vcf-w-random-sample --fam_dir [str] --fam_filename [str] --samplesize [int] --vcf_input [str] --symbol [str/int]
$ gprs random_draw_samples_from_fam --fam_dir [str] --fam_filename [str] --samplesize [int] --tag [str]
```


## Commands in gprs package:

:octocat: Thirteen commands in gprs:

1. `geneatlas-filter-data`

2. `gwas-filter-data`

3. `generate-plink-bfiles`

4. `clump`

5. `select-clump-snps`  

6. `build-prs`

7. `combine-prs`

8. `prs-statistics`

9. `combine-prs-stat`

10. `subset_pop` (optional)

11. `generate_plink_bfiles_w_individual_info` (optional)

12. `transfer_atcg` (optional)

13. `subset_vcf_w_random_sample` (optional)

### Result folder
In the first step, you need to indicate the path to creating the result folder.
Five folders will automatically generate under the result folder by script. 

- qc folder: `./result/qc/`
- snplists folder : `./result/snplists/`
- bfile folder: `./result/plink/bfiles`
- clump folder: `./result/plink/clump`
- qc_and_clump_snpslist folder: `./result/plink/qc_and_clump_snpslist`
- prs folder: `./result/plink/prs`
- pop folder: `./result/pop/`
- random_draw_sample folder: `./result/random_draw_sample/`

:heavy_exclamation_mark: Users have to indicate reference and result directories every time when using the command interface.

:heavy_exclamation_mark: Users have to provide output_name every time when they execute the commands. The output_name should be the same in every execution.

### Output file format
This package will generate output files below:
- `*.QC.csv` 
- `*.csv` 
- `*.bim`
- `*.bed`
- `*.fam`
- `*.clump`
- `*.qc_clump_snpslist.csv`
- `*.sscore`
- `*_stat.txt`
- `*_combined_stat.txt`

All output files will be named as: `[chrnb]_[name].[extension]`. 
The chrnb will given automatically, users only have to give `[name]` while using the package.
Thus, it is better use the same output name to generate all files.

- `--output_name` 
- `--snplist_name`
- `--qc_file_name`
- `--clump_file_name`
- `--plink_bfile_name`


### `gprs geneatlas-filter-data`

Filter GeneAtlas csv file by P-value and unify the data format as following order:
SNPID, ALLELE,  BETA,  StdErr, Pvalue

#### Options:
```
  --ref                     path to population reference panel  [required]
  --data_dir                The directory of GeneAtlas csv files (all 1-24 chr) [required]
  --result_dir              path to output folder; default:[./result]
  --snp_id_header           SNP ID column name in GeneAtlas original file  [required]
  --allele_header           ALLELE column name in GeneAtlas original file  [required]
  --beta_header             BETA column name in GeneAtlas original file [required]
  --se_header               StdErr column name in GeneAtlas original file  [required]
  --pvalue_header           P-value column name in GeneAtlas original file  [required]
  --output_name             output name; default: "geneatlas"; the output file name is [chrnb]_[output_name].csv and [chrnb]_[output_name].QC.csv
  --pvalue                  P-value threshold
  --help                    Show this message and exit.
````

#### Result:

This option generates two types of output in `qc` and `snplists` folders:

- `*.QC.csv` (QC files )
- `*.csv` (snplist)

 
### `gprs gwas-filter-data`

Filter GeneAtlas csv file by P-value and unify the data format as following order:
SNPID, ALLELE,  BETA,  StdErr, Pvalue

#### Options:

````
  --ref                      path to population reference panel  [required]
  --data_dir                 path to GWAS catalog summary statistic csv file (all 1-24 chr)  [required]
  --result_dir               path to output folder; default:[./result]
  --snp_id_header            SNP ID column name in GWAS catalog original file  [required]
  --allele_header            ALLELE column name in GWAS catalog original file  [required]
  --beta_header              BETA column name in GWAS catalog original file  [required]
  --se_header                StdErr column name in GWAS catalog original file  [required]
  --pvalue_header            P-value column name in GWAS catalog original file  [required]
  --output_name              output name; default: "gwas"; the output file named is [chrnb]_[output_name].csv and [chrnb]_[output_name].QC.csv
  --pvalue                   P-value threshold for filtering SNPs
  --help                     Show this message and exit.
````

#### Result:

This option generates two types of output in `qc` and `snplists` folders:

- `*.QC.csv` (QC files )
- `*.csv` (snplist)


### `gprs generate-plink-bfiles`
This option encodes plink1.9 make-bed function
```
plink --vcf [ref] --extract [snplists after qc] --make-bed --out [bfile folder/output_name]
```
snplists and bfiles folders will automatically be filled in the script.
Users have to indicate ref and output_name only.

#### Options:
````
  --ref                path to population reference panel  [required]
  --output_name        output name
  --symbol             indicate the symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz" 
  --snplist_name       snplist_name is [output_name] from [chrnb]_[output_name].csv [required]
  --extra_commands     put your extract function if needed, otherwise no need to specify this argument,
  --help               Show this message and exit.
````

#### Result:

This option will generate three files in `bfiles` folder:

- `*.bim`
- `*.bed`
- `*.fam`

### `gprs clump`
This option encodes plink1.9 clump function
```
plink --bfile [bfiles] --clump [qc snpslists] --clump-p1  --clump-p2  --clump-r2  --clump-kb  --clump-field  --clump-snp-field  --out 
```
The plink_bfiles_dir, qc snpslists and clump_output_dir will automatically be filled in the script.
Users have to indicate the options below.

#### Options:
 ````
  --clump_kb                 distance(kb) parameter for clumping [required]
  --clump_p1                 first set of P-value for clumping [required]
  --clump_p2                 should equals to p1 reduce the snps [required]
  --qc_file_name             qc_file_name is [output_name] from [chrnb]_[output_name].QC.csv [required]
  --plink_bfile_name         plink_bfile_name is [output_name] from [chrnb]_[output_name].bim/bed/fam [required]
  --output_name              it is better if the output_name remain the same. The clump output: [chrnb]_[output_name]_clumped_snplist.csv [required]
  --clump_r2                 r2 value for clumping, default = 0.1
  --clump_field              P-value column name, default = Pvalue
  --clump_snp_field          SNP ID column name, default = SNPID
  --help                     Show this message and exit.
````

#### Result:

This option will generate one file in `clump` folder:

- `*.clump`


### `gprs select-clump-snps`

#### Options:
```` 
  --qc_dir                       path to qc folder, default: "./result/qc", the qc files were generated from gwas_filter_data or geneatlas_filter_data options
  --clump_output_dir             path to clump output folder, default: "./result/plink/clump"
  --qc_clump_snplists_dir        path to snpslist (after qc and clumping), default:"./result/plink/qc_and_clump_snpslist"
  --clump_file_name              clump_file_name is [output_name] from [chrnb]_[output_name].clump [required]
  --output_name                  it is better if the output_name remain the same. output: [chrnb]_[output_name]_clumped_snplist.csv [required]
  --clump_kb                     distance(kb) parameter for clumping [required]
  --clump_p1                     first set of P-value for clumping [required]
  --clump_r2                     r2 value for clumping, default = 0.1
  --help                         Show this message and exit.
````

#### Result:
This option will generate one file in `qc_and_clump_snpslist` folder:
- `*.qc_clump_snpslist.csv`


### `gprs build-prs`
This option encodes plink2.0 function
```
plink2 --vcf [vcf input] dosage=DS --score [snplists afte clumped and qc]  --out 
```
The clumped qc snpslists and prs_output_dir will automatically be filled in the script.
Users have to indicate the options below.

#### Options:
````
  --symbol                       indicate the symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz"
  --vcf_input                    path to vcf files  [required]
  --qc_clump_snplist_foldername  the folder name of .qc_clump_snpslist.csv file. The folder name should be 
                                 the same as the output name in select_clump_snps step
  --columns                      a column index indicate the [SNPID] [ALLELE] [BETA] position; column nb starts from 1
  --plink_modifier               no-mean-imputation as default in here, get more info by searching plink2.0 modifier
  --output_name                  it is better if the output_name remain the same. output: [chrnb]_[output_name].sscore
  --clump_kb                     distance(kb) parameter for clumping [required]
  --clump_p1                     first set of P-value for clumping [required]
  --clump_r2                     r2 value for clumping, default = 0.1
  --help                         Show this message and exit.
````
#### Result:
This option will generate `.sscore` files in `prs` folder:
- `*.sscore`


### `gprs combine-prs`

Combine-prs option will combine all .sscore files as one .sscore file.

#### Options:
````
  --filename          name of .sscore, i.e.  chr10_geneatlas_500_1e-7_0.05.sscore, The file name here is "geneatlas"
  --clump_kb          distance(kb) parameter for clumping [required]
  --clump_p1          first set of P-value for clumping [required]
  --clump_r2          r2 value for clumping, default = 0.1
  --help              Show this message and exit.
````
#### Result:
This option will generate one file in `prs` folder:
- `*.sscore`


### `gprs prs-statistics`

After obtained combined sscore file, `prs-statistics` calculate BETA, AIC, AUC, PseudoR2 and OR ratio 

#### Options:
````
  --score_file             the absolute path to combined .sscore file [required]
  --pheno_file             the absolute path to pheno file  [required]
  --output_name            the output name  [required]
  --data_set_name          the name of the data-set i.e. gout_2019_GCST008970
                           [required]
  --clump_kb               distance(kb) parameter for clumping [required]
  --clump_p1               first set of P-value for clumping [required]
  --clump_r2               r2 value for clumping, default = 0.1
  --prs_stats_R            the absolute path to "prs_stats.R"  [required]
  --r_command              use "which R" in linux, and copy the path after
                           --r_command  [required]
  --help                   Show this message and exit.
````
#### Result:
This option will generate .txt file in `stat` folder:
- `*_stat.txt`


### `gprs combine-prs-stat`

If you have more than one trained PRS model, `combine-prs-stat` function is designed to combine statistics results.
For instance: the first PRS model was filtered with P < 0.05, the second PRS model was filtered with P < 0.0005. You will have DATA_0.05_stat.txt/DATA_0.0005_stat.txt
Combining two statistic tables allows users easy to compare between PRS models

#### Options:
````
  --data_set_name        the name of the data-set i.e. gout_2019_GCST008970 [required]
  --clump_kb             distance(kb) parameter for clumping [required]
  --clump_p1             first set of P-value for clumping [required]
  --clump_r2             r2 value for clumping, default = 0.1
  --help                 Show this message and exit.
````
#### Result:
This option will generate .txt file in `stat` folder:
- `*_combined_stat.txt`