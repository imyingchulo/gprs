

For example:

```
REF = "/home1/ylo40816/1000genomes/hg19"
GWAS = "/home1/ylo40816/Projects/GPRS/data/2019_GCST008970"
GENE_ATLAS = "/home1/ylo40816/Projects/GPRS/data/Gene_ATLAS/selfReported_n_1526"
```

```shell
$ python abc.py --ref='1000genomes/hg19' --gwas='PRS_project/gout/Gene_ATLAS/selfReported_n_1526' --gene_atlas='PRS_project/gout/Gene_ATLAS/selfReported_n_1526'
``` 

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

## Execute script

```shell
$ gprs
```

#### Help message

```shell
$ gprs --help
```

