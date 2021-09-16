# Setp 3-1: Clumping (remove linked SNPs)

## Function: `gprs clump`

This option encodes plink1.9 clump function

```
plink --bfile [bfiles] --clump [qc snpslists] --clump-p1  --clump-p2  --clump-r2  --clump-kb  --clump-field  --clump-snp-field  --out 
```

The plink_bfiles_dir, qc snpslists and clump_output_dir will automatically be filled in the script.
Users have to indicate the options below.

## How to use it?

Shell:

```shell
$ gprs clump --ref [str] --data_dir [str] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_p2 [float/scientific notation] --clump_r2 [float] --clump_field [str] --clump_snp_field [str] --plink_bfile_name [str] --qc_file_name [str] --output_name [output name]
```

Python:

```python
from gprs.gene_atlas_model import GeneAtlasModel

if __name__ == '__main__':
    geneatlas = GeneAtlasModel( ref='1000genomes/hg19',
                    data_dir='data/2014_GWAS_Height' )

    geneatlas.clump(output_name='2014height',
                    clump_kb='250',
                    clump_p1='0.02', clump_p2='0.02',
                    qc_file_name='2014height',
                    plink_bfile_name='2014height')


```

## output files
- `*.clump`

|CHR|F|SNP|BP|P|TOTAL|NSIG|S05|S01|S001|S0001|SP2|
|---|---|---|---|---|---|---|---|---|---|---|---|
|1   | 1   |rs1967017 | 145723645  | 3.72e-16   |    19  |    0  |    2   |   6  |    3  |    8 |rs11590105(1),rs17352281(1),rs9728345(1),rs11587821(1)|
|1   | 1  |  rs760077  |155178782 |  7.45e-10    |   12   |   0   |   2  |    2   |   1 |     7 |rs11589479(1),rs3766918(1),rs4625273(1),rs4745(1),rs12904(1)|


# Setp 3-2: Filter SNPs depends on `.clump` 
After clumping, we have to filter SNPs again, to remove linked SNPs.
In this step, we will have new SNPs list, and use it for generate PRS model.

## Function: `gprs select-clump-snps`

## How to use it?

Shell:

```shell
$ gprs select-clump-snps --result_dir [str] --qc_file_name [str] --clump_file_name [str] --clump_kb [int] --clump_p1 [float/scientific notation] --clump_r2 [float] --output_name [output name]
```

Python:

```python
from gprs.gene_atlas_model import GeneAtlasModel

if __name__ == '__main__':
    geneatlas = GeneAtlasModel( ref='1000genomes/hg19',
                    data_dir='data/2014_GWAS_Height' )

    geneatlas.select_clump_snps(output_name='2014height',clump_file_name='2014height',
                           qc_file_name='2014height')

```

## output files

- `*.qc_clump_snpslist.csv`

|CHR |SNP|
|---|---|
|1 |rs1967017|
|1 |rs760077|