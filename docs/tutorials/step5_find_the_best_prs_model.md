# Setp 5-1: Calculate statistics value of PRS model
(Rscript create by Soyoung Jeon; update by Ying-Chu Lo)

## Function: `prs-statistics`

After obtained combined sscore file, `prs-statistics` calculates BETA, AIC, AUC, PseudoR2 and OR ratio 

## How to use it?

Shell:

```shell
$ gprs prs-statistics --score_file [str] --pheno_file [str] --data_set_name [str] --prs_stats_R [str] --r_command [str] --output_name [str] 
```

Python:

```python
from gprs.gene_atlas_model import GeneAtlasModel

if __name__ == '__main__':
    geneatlas = GeneAtlasModel( ref='1000genomes/hg19',
                    data_dir='data/2014_GWAS_Height' )
    
    geneatlas.prs_statistics(output_name='2014height', score_file = "/home1/ylo40816/Projects/GPRS/tmp/2014height_250_0.02_0.1.sscore",
        pheno_file = "Projects/GPRS/tmp/result/plink/prs/2014height_pheno.csv",
        r_command='/spack/apps/linux-centos7-x86_64/gcc-8.3.0/r-4.0.0-jfy3icn4kexk7kyabcoxuio2iyyww3o7/bin/Rscript',
        prs_stats_R="Projects/GPRS/gprs/prs_stats_quantitative_phenotype.R", data_set_name="2014height",
                             clump_kb='250',
                             clump_p1='0.02',
                             clump_r2='0.1'
                             )
```

## output files

- `*_stat.txt`

| data       | filter_condition | snps_nb | P        | BETA     | Degree_of_freedom | AIC      | AUC      | PseudoR2 | OR_top1_to_middle20 | OR_top2_to_middle20 | OR_top5_to_middle20 | OR_top10_to_middle20 | 
|------------|------------------|---------|----------|----------|-------------------|----------|----------|----------|---------------------|---------------------|---------------------|----------------------|
| 2014height | 500_5e-6_0.5     | 143     | 3.97E-11 | 2.01E-01 | 1007              | 8.62E+02 | 7.47E-01 | 2.57E-01 | 4.50E+00            | 1.14E-01            | 9.59E-02            | 1.45E-01             |

# Setp 5-2: Combined PRS statistics results (Optional)

## Function:`combine-prs-stat`

If you have more than one trained PRS model, `combine-prs-stat` function is designed to combine statistics results.
For instance: the first PRS model was filtered with P < 0.05, the second PRS model was filtered with P < 0.0005. You will have DATA_0.05_stat.txt/DATA_0.0005_stat.txt
Combining two statistic tables allows users easy to compare between PRS models

## How to use it?

Shell:

```shell
$ gprs combine-prs-stat --data_set_name [str]
```

Python:

```python
from gprs.gene_atlas_model import GeneAtlasModel

if __name__ == '__main__':
    geneatlas = GeneAtlasModel( ref='1000genomes/hg19',
                    data_dir='data/2014_GWAS_Height' )

    geneatlas.combine_prs_stat(data_set_name='2014height')
```

## output files

- `*_combined_stat.txt`

|data|filter_conditions |P| BETA|PseudoR2|OR_top1_to_middle20|OR_top2_to_middle20|OR_top5_to_middle20|OR_top10_to_middle20|
|---|---|---|---|---|---|---|---|---|
|2014height |250_0.02_0.1    | 0.5632    |0.0019|0.0002|0.5396|1.3899|0.9408| 1.2972|
|2014height |250_0.02_0.1    | 0.5632    |0.0019|0.0002|0.5396|1.3899|0.9408| 1.2972|



