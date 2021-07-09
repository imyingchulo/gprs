import click
from gprs import GPRS
from gprs.gwas import Gwas
from gprs.gene_atlas import GeneAtlas


@click.group()
def main():
    pass


@click.command()
@click.option( '--data_dir', metavar='<str>', required=True, help='a directory of GeneAtlas csv files' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder; default:[./result]' )
@click.option( '--snp_id_header', metavar='<str>', required=True, help='SNP ID column name' )
@click.option( '--allele_header', metavar='<str>', required=True, help='ALLELE column name' )
@click.option( '--beta_header', metavar='<str>', required=True, help='BETA column name' )
@click.option( '--se_header', metavar='<str>', required=True, help='StdErr column name' )
@click.option( '--pvalue_header', metavar='<str>', required=True, help='P-value column name' )
@click.option( '--output_name', metavar='<str>', default='geneatlas',
               help='output name; default: geneatlas; the output file will be [chrnb]_[output_name].csv' )
@click.option( '--pvalue', metavar='<float/scientific notation>', default=0.05, help='P-value threshold' )
def geneatlas_filter_data(data_dir, result_dir, snp_id_header,
                          allele_header,
                          beta_header,
                          se_header,
                          pvalue_header,
                          output_name,
                          pvalue):
    gene_atlas = GeneAtlas( data_dir=data_dir, result_dir=result_dir )
    gene_atlas.filter_data( snp_id_header=snp_id_header,
                            allele_header=allele_header,
                            beta_header=beta_header,
                            se_header=se_header,
                            pvalue_header=pvalue_header,
                            output_name=output_name,
                            pvalue=pvalue )


@click.command()
@click.option( '--data_dir', metavar='<str>', required=True,
               help='path to GWAS catalog summary statistic csv file [1-24 chr]' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder; default:[./result]' )
@click.option( '--snp_id_header', metavar='<str>', required=True, help='SNP ID column name' )
@click.option( '--allele_header', metavar='<str>', required=True, help='ALLELE column name' )
@click.option( '--beta_header', metavar='<str>', required=True, help='BETA column name' )
@click.option( '--se_header', metavar='<str>', required=True, help='StdErr column name' )
@click.option( '--pvalue_header', metavar='<str>', required=True, help='P-value column name' )
@click.option( '--output_name', metavar='<str>', default='gwas',
               help='output name; default: gwas ; the output file will be named as [chrnb]_[output_name].csv' )
@click.option( '--pvalue', metavar='<float/scientific notation>', default=0.05,
               help='P-value threshold for filtering SNPs' )
def gwas_filter_data(data_dir, result_dir, snp_id_header,
                     allele_header,
                     beta_header,
                     se_header,
                     pvalue_header,
                     output_name,
                     pvalue):
    gwas = Gwas( data_dir=data_dir, result_dir=result_dir )
    gwas.filter_data( snp_id_header=snp_id_header,
                      allele_header=allele_header,
                      beta_header=beta_header,
                      se_header=se_header,
                      pvalue_header=pvalue_header,
                      output_name=output_name,
                      pvalue=pvalue )


@click.command()
@click.option( '--ref', metavar='<str>', required=True, help='path to population reference panel' )
@click.option( '--plink_bfiles_dir', metavar='<str>', required=True, help='output folder' )
@click.option( '--snplists_dir', metavar='<str>', required=True,
               help='snplists the name of the file should be chrnb_[output_name].csv '
                    'i.e. chr1_[output_name].csv; please check "grps gwas_filter_data --help or grps geneatlas_filter_data --help' )
@click.option( '--output_name', metavar='<str>', help='output name' )
def generate_plink_bfiles(ref, plink_bfiles_dir, snplists_dir, output_name):
    gprs = GPRS( ref=ref )
    gprs.generate_plink_bfiles( plink_bfiles_dir=plink_bfiles_dir,
                                snplists_dir=snplists_dir,
                                output_name=output_name )


@click.command()
@click.option( '--data_dir', metavar='<str>', required=True,
               help='path to GWAS catalog summary statistics/ GeneAtlas csv file' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder' )
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<int>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_p2', metavar='<int>', required=True, help='second set of P-value for clumping' )
@click.option( '--output_name', metavar='<str>', required=True,
               help='output name should remain consistent as output_name to plink and filtered data; output format is: [chrnb]_[output_name]_clumped_snplist.csv' )
@click.option( '--clump_r2', metavar='<int>', default=0.1, required=True, help='r2 value for clumping, default = 0.1' )
@click.option( '--clump_field', metavar='<str>', default='Pvalue', required=True,
               help='p value column name, default = Pvalue' )
@click.option( '--clump_snp_field', metavar='<str>', default='SNPID', required=True,
               help='SNP ID column name, default = SNPID' )
def clump(data_dir, result_dir, clump_kb, clump_p1, clump_p2, output_name, clump_r2, clump_field, clump_snp_field):
    gprs = GPRS( data_dir=data_dir, result_dir=result_dir )
    gprs.clump( clump_kb=clump_kb,
                clump_p1=clump_p1,
                clump_p2=clump_p2,
                output_name=output_name,
                clump_r2=clump_r2,
                clump_field=clump_field,
                clump_snp_field=clump_snp_field )
    gprs.select_clump_snps( output_name=output_name )


@click.command()
@click.option( '--ref', metavar='<str>', required=True, help='path to genome reference' )
@click.option( '--prs_output_dir', metavar='<str>', required=True, help='path to prs model output' )
@click.option( '--qc_clump_snplists_dir', metavar='<str>', required=True, help='path to snpslist' )
@click.option( '--columns', metavar='<int>', default='1 2 3', required=True,
               help='a column index indicate the [SNPID] [ALLELE] [BETA] position ' )
@click.option( '--plink_modifier', metavar='<str>', default='no-mean-imputation', required=True,
               help='no-mean-imputation as default in here, get more info by searching plink2.0 modifier ' )
@click.option( '--output_name', metavar='<str>', required=True,
               help='output name should remain consistent as output_name to plink and filtered data' )
def build_prs(ref, prs_output_dir, qc_clump_snplists_dir, columns, plink_modifier, output_name, ):
    gprs = GPRS( ref=ref )
    gprs.build_prs( prs_output_dir=prs_output_dir,
                    qc_clump_snplists_dir=qc_clump_snplists_dir,
                    columns=columns,
                    plink_modifier=plink_modifier,
                    output_name=output_name )


main.add_command( geneatlas_filter_data )
main.add_command( gwas_filter_data )
main.add_command( generate_plink_bfiles )
main.add_command( clump )
main.add_command( build_prs )
