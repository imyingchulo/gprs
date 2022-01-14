# Get started: understanding the data format


## GWAS template

|Chr|Pos|RSID|Allele1|Allele2|Freq1|Effect|StdErr|P-value|n_total_sum|
|---|---|---|---|---|---|---|---|---|---|
|1 |10539 |rs537182016 |a |c |0.0013 |-5.1213 |20.0173 |0.7981 |7043|
|1 |11008 |rs575272151 |c |g |0.9215 |0.1766 |0.2610 |0.4985 |7042.99|
|1 |11012 |rs544419019 |c |g |0.9215 |0.1766 |0.2610 |0.4985 |7042.99| 
|1 |14674 |rs561913721 |a |g |0.0032 |0.8040 |0.8364 |0.3364 |7043|


## GeneAtlas template

|SNP |ALLELE |NBETA-selfReported_n_1526 |NSE-selfReported_n_1526 |PV-selfReported_n_1526
|---|---|---|---|---|
|rs1110052 |T |-0.00032386 |0.000269 |0.2286|
|rs112164716 |T |6.5121e-05 |0.001126 |0.95388|
|rs11240779|A |-0.00022416 |0.00028937 |0.43855|
|rs11260596 |C |0.00025168 |0.00024248 |0.29931|


## How to choose model template?

|-|chr info in the file name| chr info not in the file name|
|---|---|---|
|chr info in the header|gene_atlas_model|gwas_model|
|chr info absent in the header|gene_atlas_model|gene_atlas_model|

# Before step1:

Please unzip your .gz file first.

# Step1: Unify the data format
After knowing the data format, users can choose the model (gwas or geneatlas) to unify the data format and filter out SNPs(optional).
:heavy_exclamation_mark: SNPs are extract out by RSID not chromosome position

## Function: `gprs geneatlas-filter-data`

Filter GeneAtlas csv file by P-value and unify the data format as following order:
SNPID, ALLELE,  BETA,  StdErr, Pvalue

## How to use it?

Shell:

```shell
$ gprs geneatlas-filter-data --ref [str] --data_dir [str] --result_dir [str] --snp_id_header [str] --allele_header [str] --beta_header [str] --se_header [str] --pvalue_header [str] --pvalue [float/scientific notation] --output_name [str]  
$ gprs gwas-filter-data --ref [str] --data_dir [str] --result_dir [str] --snp_id_header [str] --allele_header  [str] --beta_header [str] --se_header [str] --pvalue_header [str] --pvalue [float/scientific notation] --output_name [str]  
```

Python:

```python
from gprs.gene_atlas_model import GeneAtlasModel
if __name__ == '__main__':
    geneatlas = GeneAtlasModel( ref='1000genomes/hg19',
                    data_dir='data/2014_GWAS_Height' )

    geneatlas.filter_data( snp_id_header='MarkerName',
                            allele_header='Allele1',
                            beta_header='b',
                            se_header ='SE',
                            pvalue_header='p',
                            output_name='2014height')
   
from gprs.gwas_model import GwasModel
if __name__ == '__main__':
    gwas = GwasModel( ref='/home1/ylo40816/1000genomes/hg19',
                 data_dir='/home1/ylo40816/Projects/GPRS/data/2019_GCST008970')

    gwas.filter_data( snp_id_header='RSID',
                   allele_header='Allele1',
                   beta_header='Effect',
                   se_header='StdErr',
                   pvalue_header='P-value',
                   output_name='GCST008970',
                   file_name='gout_chr1_22_LQ_IQ06_mac10_all_201_rsid.csv')
```

## output files
- `*.QC.csv` (QC files )
- `*.csv` (snplist)