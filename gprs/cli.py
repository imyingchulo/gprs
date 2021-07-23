import click
from gprs import GPRS
from gprs.gwas import Gwas
from gprs.gene_atlas import GeneAtlas


@click.group()
def main():
    pass

@click.command()
@click.option( '--ref', metavar='<str>', help='path to population reference panel' )
@click.option( '--data_dir', metavar='<str>', required=True, help='The directory of GeneAtlas .csv files (all 1-24 chr)' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder; default:[./result]' )
@click.option( '--snp_id_header', metavar='<str>', required=True, help='SNP ID column name in GeneAtlas original file' )
@click.option( '--allele_header', metavar='<str>', required=True, help='ALLELE column name in GeneAtlas original file' )
@click.option( '--beta_header', metavar='<str>', required=True, help='BETA column name in GeneAtlas original file' )
@click.option( '--se_header', metavar='<str>', required=True, help='StdErr column name in GeneAtlas original file' )
@click.option( '--pvalue_header', metavar='<str>', required=True, help='P-value column name in GeneAtlas original file' )
@click.option( '--output_name', metavar='<str>', default='geneatlas', help='output name; default: "geneatlas"; the output file name is [chrnb]_[output_name].csv and [chrnb]_[output_name].QC.csv' )
@click.option( '--pvalue', metavar='<float/scientific notation>', default=0.05, help='P-value threshold' )
def geneatlas_filter_data(ref,data_dir, result_dir, snp_id_header,
                          allele_header,
                          beta_header,
                          se_header,
                          pvalue_header,
                          output_name,
                          pvalue):
    gene_atlas = GeneAtlas( ref=ref, data_dir=data_dir, result_dir=result_dir )
    gene_atlas.filter_data( snp_id_header=snp_id_header,
                            allele_header=allele_header,
                            beta_header=beta_header,
                            se_header=se_header,
                            pvalue_header=pvalue_header,
                            output_name=output_name,
                            pvalue=pvalue )

@click.command()
@click.option( '--ref', metavar='<str>', help='path to population reference panel' )
@click.option( '--data_dir', metavar='<str>', required=True, help='path to GWAS catalog summary statistic .csv file (all 1-24 chr)' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder, default: "./result"' )
@click.option( '--snp_id_header', metavar='<str>', required=True, help='SNP ID column name in GWAS catalog original file' )
@click.option( '--allele_header', metavar='<str>', required=True, help='ALLELE column name in GWAS catalog original file' )
@click.option( '--beta_header', metavar='<str>', required=True, help='BETA column name in GWAS catalog original file' )
@click.option( '--se_header', metavar='<str>', required=True, help='StdErr column name in GWAS catalog original file' )
@click.option( '--pvalue_header', metavar='<str>', required=True, help='P-value column name in GWAS catalog original file' )
@click.option( '--output_name', metavar='<str>', default='gwas', help='output name; default: "gwas"; the output file name is [chrnb]_[output_name].csv and [chrnb]_[output_name].QC.csv' )
@click.option( '--pvalue', metavar='<float/scientific notation>', default=0.05, help='P-value threshold for filtering SNPs' )
def gwas_filter_data(ref, data_dir, result_dir, snp_id_header,
                     allele_header,
                     beta_header,
                     se_header,
                     pvalue_header,
                     output_name,
                     pvalue):
    gwas = Gwas( ref=ref, data_dir=data_dir, result_dir=result_dir )
    gwas.filter_data( snp_id_header=snp_id_header,
                      allele_header=allele_header,
                      beta_header=beta_header,
                      se_header=se_header,
                      pvalue_header=pvalue_header,
                      output_name=output_name,
                      pvalue=pvalue )

@click.command()
@click.option( '--ref', metavar='<str>', required=True, help='path to population reference panel' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder, default: "./result"' )
@click.option( '--qc_file_name', metavar='<str>', help='qc_file_name is [output_name] from [chrnb]_[output_name].QC.csv' )
def transfer_atcg(ref, result_dir, qc_file_name):
    gprs = GPRS( ref=ref, result_dir=result_dir )
    gprs.transfer_atcg( qc_file_name=qc_file_name )


@click.command()
@click.option( '--ref', metavar='<str>', required=True, help='path to population reference panel' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder, default: "./result"' )
@click.option( '--symbol', metavar='<str/int>', default='.', help='indicate the symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz"' )
@click.option( '--snplist_name', metavar='<str>', required=True, help='snplist_name is [output_name] from [chrnb]_[output_name].csv' )
@click.option( '--output_name', metavar='<str>', help='it is better if the output name should be the same as snplist file name' )
def generate_plink_bfiles(ref, result_dir, snplist_name, output_name, symbol):
    gprs = GPRS( ref=ref, result_dir=result_dir )
    gprs.generate_plink_bfiles( output_name=output_name, symbol = symbol, snplist_name=snplist_name )


@click.command()
@click.option( '--ref', metavar='<str>', help='path to population reference panel' )
@click.option( '--data_dir', metavar='<str>', required=True, help='path to GWAS catalog/ GeneAtlas .csv file' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder, default: "./result"' )
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_p2', metavar='<float/scientific notation>', required=True, help='should equals to p1 reduce the snps' )
@click.option( '--qc_file_name', metavar='<str>', required=True, help='qc_file_name is [output_name] from [chrnb]_[output_name].QC.csv' )
@click.option( '--plink_bfile_name', metavar='<str>', required=True, help='plink_bfile_name is [output_name] from [chrnb]_[output_name].bim/bed/fam' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output_name remain the same. The clump output: [chrnb]_[output_name]_clumped_snplist.csv' )
@click.option( '--clump_r2', metavar='<float>', default=0.1, help='r2 value for clumping, default = 0.1' )
@click.option( '--clump_field', metavar='<str>', default='Pvalue', help='P-value column name, default = Pvalue' )
@click.option( '--clump_snp_field', metavar='<str>', default='SNPID', help='SNP ID column name, default = SNPID' )
def clump(ref, data_dir, qc_file_name, result_dir, plink_bfile_name, clump_kb, clump_p1, clump_p2, output_name, clump_r2, clump_field, clump_snp_field):
    gprs = GPRS( ref=ref, data_dir=data_dir, result_dir=result_dir)
    gprs.clump( qc_file_name=qc_file_name,
                plink_bfile_name=plink_bfile_name,
                clump_kb=clump_kb,
                clump_p1=clump_p1,
                clump_p2=clump_p2,
                output_name=output_name,
                clump_r2=clump_r2,
                clump_field=clump_field,
                clump_snp_field=clump_snp_field )


@click.command()
@click.option( '--ref', metavar='<str>', help='path to population reference panel' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder, default: "./result"' )
@click.option( '--qc_file_name', metavar='<str>', required=True, help='qc_file_name is [output_name] from [chrnb]_[output_name].QC.csv' )
@click.option( '--clump_file_name', metavar='<str>', required=True, help='clump_file_name is [output_name] from [chrnb]_[output_name].clump' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output_name remain the same. output: [chrnb]_[output_name]_clumped_snplist.csv' )
def select_clump_snps(ref, result_dir, clump_file_name, qc_file_name,output_name):
    gprs = GPRS( ref=ref, result_dir=result_dir )
    gprs.select_clump_snps( qc_file_name=qc_file_name, clump_file_name=clump_file_name, output_name=output_name )


@click.command()
@click.option( '--ref', metavar='<str>', help='path to population reference panel' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder, default: "./result"' )
@click.option( '--vcf_input', metavar='<str>', required=True, help='path to vcf files' )
@click.option( '--symbol', metavar='<str/int>', default='.', help='indicate the symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz"' )
@click.option( '--columns', metavar='<int>', default='1 2 3', help='a column index indicate the [SNPID] [ALLELE] [BETA] position; column nb starts from 1 ' )
@click.option( '--plink_modifier', metavar='<str>', default='no-mean-imputation', help='no-mean-imputation as default in here, get more info by searching plink2.0 modifier ' )
@click.option( '--qc_file_name', metavar='<str>', required=True, help='qc_file_name is [output_name] from [chrnb]_[output_name].QC.csv' )
@click.option( '--memory', metavar='<int>', help='number of memory use' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output_name remain the same. output: [chrnb]_[output_name].sscore' )
def build_prs(ref, result_dir, vcf_input, qc_file_name, columns, plink_modifier, output_name, symbol,memory ):
    gprs = GPRS( ref=ref, result_dir=result_dir )
    gprs.build_prs( memory=memory,
                    vcf_input=vcf_input,
                    qc_file_name=qc_file_name,
                    symbol = symbol,
                    columns=columns,
                    plink_modifier=plink_modifier,
                    output_name=output_name )


main.add_command( geneatlas_filter_data )
main.add_command( gwas_filter_data )
main.add_command( transfer_atcg )
main.add_command( generate_plink_bfiles )
main.add_command( clump )
main.add_command( select_clump_snps )
main.add_command( build_prs )


