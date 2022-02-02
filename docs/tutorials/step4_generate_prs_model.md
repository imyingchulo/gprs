# Setp 4-1: Generate PRS model
Generate PRS model by using Dosage by plink2.0

## Function: `gprs build-prs`

```
plink2 --vcf [vcf input] dosage=DS --score [snplists afte clumped and qc]  --out 
```

The clumped qc snpslists and prs_output_dir will automatically be filled in the script.
Users have to indicate the options below.

## How to use it?

Shell:

```shell
$ gprs build-prs --vcf_input [str] --qc_clump_snplist_foldername [str] --symbol [str/int] --columns [int] --plink_modifier [str] --memory [int] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_r2 [float] --output_name [output name]
```

Python:

```python
from gprs.gene_atlas_model import GeneAtlasModel

if __name__ == '__main__':
    geneatlas = GeneAtlasModel( ref='1000genomes/hg19',
                    data_dir='data/2014_GWAS_Height' )

    geneatlas.build_prs( vcf_input= '1000genomes/hg19',
                          output_name ='2014height', memory='1000',clump_kb='250',
                    clump_p1='0.02', clump_r2='0.02', qc_clump_snplist_foldername='2014height')
```

## output files

- `*.sscore`

|IID|NMISS_ALLELE_CT |NAMED_ALLELE_DOSAGE_SUM |SCORE1_AVG|
|---|---|---|---|
|HG00096 |130     |116    |-0.00131555|
|HG00097 |130     |114     |-0.00133961|
|HG00099 |130     |119    |-0.00137718|
|HG00100 |130     |110     |-0.00125486|


# Setp 4-2: Combined PRS model
Combined PRS model (python script create by Soyoung Jeon; update by Ying-Chu Lo))

## Function:`gprs combine-prs`

Combine-prs will combine all .sscore files as one .sscore file.
And calculate score average and sum per individual.

## How to use it?

Shell:

```shell
$ gprs combine-prs --ref [str] --result_dur [str] 
```

Python:

```python
from gprs.gene_atlas_model import GeneAtlasModel

if __name__ == '__main__':
    geneatlas = GeneAtlasModel( ref='1000genomes/hg19',
                    data_dir='data/2014_GWAS_Height' )

    geneatlas.combine_prs(filename="2014height",clump_r2="0.5",clump_kb="250",clump_p1="0.02")
```

## output files

- `*.sscore`

|id      |ALLELE_CT       |SCORE_AVG       |SCORE_SUM|
|---|---|---|---|
|HG03270 |1872.0  |-0.00109666512273       |-2.2826939078|
|HG03271 |1872.0  |-0.00111419935  |-2.2831272058|
|NA19670 |1872.0  |-0.00117191923182       |-2.4014961794|
|HG03279 |1872.0  |-0.00115016386364       |-2.3057819|