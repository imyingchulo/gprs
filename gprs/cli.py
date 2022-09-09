import click
from gprs.gprs_main import GPRS
from gprs.gwas_model import GwasModel
from gprs.gene_atlas_model import GeneAtlasModel

@click.group()
def main():
    pass

# @click.command()
# def test():
#     print("hello world")

@click.command()
@click.option( '--ref', metavar='<str>', help='path to population reference panel' )
@click.option( '--data_dir', metavar='<str>', required=True, help='The directory of GeneAtlas .csv files (all chromosomes)' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder; default:[./result]' )
@click.option( '--snp_id_header', metavar='<str>', required=True, help='SNP ID column name in GeneAtlas original file' )
@click.option( '--allele_header', metavar='<str>', required=True, help='ALLELE column name in GeneAtlas original file' )
@click.option( '--beta_header', metavar='<str>', required=True, help='BETA column name in GeneAtlas original file' )
@click.option( '--se_header', metavar='<str>', required=True, help='StdErr column name in GeneAtlas original file' )
@click.option( '--pvalue_header', metavar='<str>', required=True, help='P-value column name in GeneAtlas original file' )
@click.option( '--output_name', metavar='<str>', default='geneatlas_model', help='output name; default: "geneatlas"; the output file name is [chrnb]_[output_name].csv and [chrnb]_[output_name].QC.csv' )
@click.option( '--pvalue', metavar='<float/scientific notation>', default=1, help='P-value threshold' )
def geneatlas_filter_data(ref,data_dir, result_dir, snp_id_header,
                          allele_header,
                          beta_header,
                          se_header,
                          pvalue_header,
                          output_name,
                          pvalue):
    gene_atlas = GeneAtlasModel( ref=ref, data_dir=data_dir, result_dir=result_dir )
    gene_atlas.filter_data( snp_id_header=snp_id_header,
                            allele_header=allele_header,
                            beta_header=beta_header,
                            se_header=se_header,
                            pvalue_header=pvalue_header,
                            output_name=output_name,
                            pvalue=pvalue )

@click.command()
@click.option( '--ref', metavar='<str>', help='path to population reference panel' )
@click.option( '--data_dir', metavar='<str>', required=True, help='path to GWAS catalog summary statistic .csv file (all chromosomes)' )
@click.option( '--result_dir', metavar='<str>', default='./result', help='path to output folder, default: "./result"' )
@click.option( '--snp_id_header', metavar='<str>', required=True, help='SNP ID column name in GWAS catalog original file' )
@click.option( '--allele_header', metavar='<str>', required=True, help='ALLELE column name in GWAS catalog original file' )
@click.option( '--beta_header', metavar='<str>', required=True, help='BETA column name in GWAS catalog original file' )
@click.option( '--se_header', metavar='<str>', required=True, help='StdErr column name in GWAS catalog original file' )
@click.option( '--pvalue_header', metavar='<str>', required=True, help='P-value column name in GWAS catalog original file' )
@click.option( '--file_name', metavar='<str>', default='gwas_model', help='raw data file name' )
@click.option( '--output_name', metavar='<str>', default='gwas_model', help='output name; default: "gwas"; the output file name is [chrnb]_[output_name].csv and [chrnb]_[output_name].QC.csv' )
@click.option( '--pvalue', metavar='<float/scientific notation>', default=1, help='P-value threshold for filtering SNPs' )
def gwas_filter_data(ref, data_dir, result_dir, snp_id_header,
                     allele_header,
                     beta_header,
                     se_header,
                     pvalue_header,
                     output_name,
                     pvalue,
                     file_name):
    gwas = GwasModel( ref=ref, data_dir=data_dir, result_dir=result_dir )
    gwas.filter_data( snp_id_header=snp_id_header,
                      allele_header=allele_header,
                      beta_header=beta_header,
                      se_header=se_header,
                      pvalue_header=pvalue_header,
                      output_name=output_name,
                      pvalue=pvalue,
                      file_name=file_name)

@click.command()
@click.option( '--ref', metavar='<str>', required=True, help='path to population reference panel' )
@click.option( '--snplist_name', metavar='<str>', required=True, help='snplist_name is [output_name] from [chrnb]_[output_name].csv' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output name should be the same as snplist file name' )
@click.option( '--symbol', metavar='<str/int>', required=True, default='.', help='indicate the symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz"' )
@click.option( '--extra_commands', metavar='<str>', default=' ', help='a space to add new functions for generate the plink bfiles' )
def generate_plink_bfiles(ref, snplist_name, output_name, symbol,extra_commands):
    gprs = GPRS( ref=ref)
    gprs.generate_plink_bfiles( output_name=output_name, symbol = symbol, snplist_name=snplist_name, extra_commands=extra_commands )

@click.command()
@click.option( '--plink_bfile_name', metavar='<str>', required=True, help='plink_bfile_name is [output_name] from [chrnb]_[output_name].bim/bed/fam' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output_name remain the same. The clump output: [chrnb]_[output_name]_clumped_snplist.csv' )
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_p2', metavar='<float/scientific notation>', required=True, help='should equals to p1 reduce the snps' )
@click.option( '--clump_r2', metavar='<float>', default=0.1, help='r2 value for clumping, default = 0.1' )
@click.option( '--clump_field', metavar='<str>', default='Pvalue', help='P-value column name, default = Pvalue' )
@click.option( '--qc_file_name', metavar='<str>', required=True, help='qc_file_name is [output_name] from [chrnb]_[output_name].QC.csv' )
@click.option( '--clump_snp_field', metavar='<str>', default='SNPID', help='SNP ID column name, default = SNPID' )
def clump(qc_file_name, plink_bfile_name, clump_kb, clump_p1, clump_p2, output_name, clump_r2, clump_field, clump_snp_field):
    gprs = GPRS()
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
@click.option( '--qc_file_name', metavar='<str>', required=True, help='qc_file_name is [output_name] from [chrnb]_[output_name].QC.csv' )
@click.option( '--clump_file_name', metavar='<str>', required=True, help='clump_file_name is [output_name] from [chrnb]_[output_name].clump' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output_name remain the same. output: [chrnb]_[output_name]_clumped_snplist.csv' )
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_r2', metavar='<float>', required=True, help='r2 value for clumping' )
@click.option('--clumpfolder_name',metavar='<str>', required=True, help='folder name for .clump files')
def select_clump_snps(clump_file_name, qc_file_name,output_name,clump_kb,clump_p1,clump_r2, clumpfolder_name):
    gprs = GPRS()
    gprs.select_clump_snps( qc_file_name=qc_file_name, clump_file_name=clump_file_name, output_name=output_name,
                            clump_kb=clump_kb,clump_p1=clump_p1,clump_r2=clump_r2,clumpfolder_name=clumpfolder_name)

@click.command()
@click.option( '--vcf_input', metavar='<str>', required=True, help='path to vcf files' )
@click.option( '--vcf_extension', metavar='<str>', required=True, help='vcf file extension' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output_name remain the same. output: [chrnb]_[output_name].sscore' )
@click.option('--qc_clump_snplist_foldername',metavar='<str>', required=True, help='folder name for .qc_clump_snpslist.csv files, i.e. LAT is the name of LAT_250_1e-5_0.5 folder')
@click.option( '--memory', metavar='<int>', help='number of memory use' )
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_r2', metavar='<float>', required=True, help='r2 value for clumping' )
@click.option( '--symbol', metavar='<str/int>', required=True,default='.', help='indicate the symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz"' )
@click.option( '--columns', metavar='<int>', default='1 2 3', help='a column index indicate the [SNPID] [ALLELE] [BETA] position; column nb starts from 1 ' )
@click.option( '--plink_modifier', metavar='<str>', default='no-mean-imputation', help='no-mean-imputation as default in here, get more info by searching plink2.0 modifier ' )
def build_prs(vcf_input, vcf_extension, columns, plink_modifier, output_name, symbol,memory,clump_kb,clump_p1,clump_r2,qc_clump_snplist_foldername ):
    gprs = GPRS()
    gprs.build_prs( memory=memory,
                    vcf_input=vcf_input,
                    vcf_extension=vcf_extension,
                    symbol = symbol,
                    columns=columns,
                    plink_modifier=plink_modifier,
                    output_name=output_name,
                    clump_kb=clump_kb,
                    clump_p1=clump_p1,
                    clump_r2=clump_r2,
                    qc_clump_snplist_foldername=qc_clump_snplist_foldername)

@click.command()
@click.option( '--filename', metavar='<str>', required=True, help='name of .sscore, i.e.  chr10_geneatlas_500_1e-7_0.05.sscore, The file name here is "geneatlas"')
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_r2', metavar='<float>', required=True, help='r2 value for clumping' )
def combine_prs(filename,clump_kb,clump_p1,clump_r2):
    gprs = GPRS()
    gprs.combine_prs( filename=filename,clump_kb=clump_kb,clump_p1=clump_p1,clump_r2=clump_r2)

@click.command()
@click.option( '--score_file', metavar='<str>', required=True, help='the absolute path to combined .sscore file')
@click.option( '--pheno_file', metavar='<str>', required=True, help='the absolute path to pheno file')
@click.option( '--output_name', metavar='<str>', required=True, help='the output name')
@click.option( '--data_set_name', metavar='<str>', required=True, help='the name of the data-set i.e. gout_2019_GCST008970 ')
@click.option( '--prs_stats_R', metavar='<str>', required=True, help='the absolute path to "prs_stats_quantitative_phenotype.R')
@click.option( '--r_command', metavar='<str>', required=True, help='use "which R" in linux, and copy the path after --r_command')
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_r2', metavar='<float>', default=0.1, help='r2 value for clumping, default = 0.1' )
def prs_statistics(score_file, pheno_file, output_name, data_set_name, prs_stats_R, r_command,clump_kb,clump_p1,clump_r2):
    gprs = GPRS()
    gprs.prs_statistics( score_file=score_file,
                         pheno_file=pheno_file,
                         output_name=output_name,
                         data_set_name=data_set_name,
                         prs_stats_R=prs_stats_R,
                         r_command=r_command,
                         clump_kb=clump_kb,clump_p1=clump_p1,clump_r2=clump_r2)

@click.command()
@click.option( '--data_set_name', metavar='<str>', required=True, help='the name of the data-set i.e. gout_2019_GCST008970' )
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_r2', metavar='<float>', default=0.1, help='r2 value for clumping, default = 0.1' )
def combine_prs_stat(data_set_name,clump_kb,clump_p1,clump_r2):
    gprs = GPRS()
    gprs.combine_prs_stat( data_set_name=data_set_name,
                           clump_kb=clump_kb,clump_p1=clump_p1,clump_r2=clump_r2)

## Optional function here
@click.command()
@click.option('--qc_file_name', metavar='<str>', help='qc_file_name is [output_name] from [chrnb]_[output_name].QC.csv')
def transfer_atcg(qc_file_name):
    gprs = GPRS()
    gprs.transfer_atcg(qc_file_name=qc_file_name)

@click.command()
@click.option('--input_data', metavar='<str>', required=True, help='The full path to phenotype file (contains the population information)')
@click.option('--column_name', metavar='<str>', required=True, help='The header of population column')
@click.option('--pop_info', metavar='<str>', required=True, help='The target population to extract')
@click.option('--output_name', metavar='<str>', required=True, help='file name for subset population')
def subset_pop(input_data, pop_info, output_name, column_name):
    gprs = GPRS()
    gprs.subset_pop(column_name=column_name,
                    input_data=input_data,
                    output_name=output_name,
                    pop_info=pop_info)

@click.command()
@click.option('--popfile_name', metavar='<str>', required=True, help='the pop file name')
@click.option('--bfile_name', metavar='<str>', required=True, help='the bfile name')
@click.option('--plink_command', metavar='<str>', required=True, help='--keep or --exclude')
@click.option('--output_name', metavar='<str>', help='output name for sub-set population bfile')
def generate_plink_bfiles_w_individual_info(popfile_name, bfile_name, output_name, plink_command):
    gprs = GPRS()
    gprs.generate_plink_bfiles_w_individual_info(output_name=output_name, bfile_name=bfile_name,
                                                 popfile_name=popfile_name, plink_command=plink_command)

@click.command()
@click.option('--fam_dir', metavar='<str>', required=True, help='the path to .fam')
@click.option('--fam_filename', metavar='<str>', required=True, help='name of fam file(without chr number)')
@click.option('--samplesize', metavar='<int>', required=True, help='number of subset samples')
@click.option('--vcf_input', metavar='<str>', required=True, help='path to vcf files')
@click.option('--symbol', metavar='<str/int>', required=True, default='.', help='indicate the symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz"' )
def subset_vcf_w_random_sample(fam_dir, fam_filename, samplesize, vcf_input, symbol):
    gprs = GPRS()
    gprs.subset_vcf_w_random_sample(fam_dir=fam_dir,
                                    fam_filename=fam_filename,
                                    samplesize=samplesize,
                                    vcf_input=vcf_input,
                                    symbol=symbol)

@click.command()
@click.option('--fam_dir', metavar='<str>', required=True, help='the path to .fam')
@click.option('--fam_filename', metavar='<str>', required=True, help='name of fam file(without chr number)')
@click.option('--samplesize', metavar='<int>', required=True, help='number of subset samples')
@click.option('--tag', metavar='<str>', required=True, help='extract tag for output file name, i.e. LD_preference, training_individuals' )
def random_draw_samples_from_fam(fam_dir, fam_filename, samplesize, tag):
    gprs = GPRS()
    gprs.random_draw_samples_from_fam(fam_dir=fam_dir,
                                    fam_filename=fam_filename,
                                    samplesize=samplesize,tag=tag)

@click.command()
@click.option('--index_of_chr_id', metavar='<int>', required=True, help='index of chr column in summary file')
@click.option('--index_of_pos', metavar='<int>', required=True, help='index of pos column in summary file')
@click.option('--index_of_allele1', metavar='<int>', required=True, help='index of REF column in summary file')
@click.option('--index_of_allele2', metavar='<int>', required=True, help='index of ALT column in summary file' )
@click.option('--index_of_beta', metavar='<int>', required=True, help='index of beta column in summary file' )
@click.option('--index_of_se', metavar='<int>', required=True, help='index of se column in summary file' )
@click.option('--index_of_pvalue', metavar='<int>', required=True, help='index of pvalue column in summary file' )
@click.option('--summary_file', metavar='<str>', required=True, help='full GWAS summary file name with extension' )
@click.option('--output_name', metavar='<str>', required=True, help='output name for uniq snplist and qc files' )
@click.option('--data_dir', metavar='<str>', required=True, help='path to summary data' )
def create_new_marker(index_of_chr_id,index_of_pos,index_of_allele1, index_of_allele2, summary_file, index_of_beta, index_of_se, index_of_pvalue, output_name, data_dir):
    gprs = GPRS(data_dir=data_dir)
    gprs.create_new_marker(index_of_chr_id=index_of_chr_id,
                           index_of_pos=index_of_pos,index_of_allele1=index_of_allele1,index_of_allele2=index_of_allele2,
                           summary_file=summary_file,output_name=output_name, index_of_beta=index_of_beta, index_of_se=index_of_se, index_of_pvalue=index_of_pvalue)

@click.command()
@click.option( '--data_set_name', metavar='<str>', required=True, help='the name of the data-set i.e. gout_2019_GCST008970' )
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='P-value for clumping' )
@click.option( '--clump_r2', metavar='<float>',required=True, help='r2 value for clumping' )
@click.option('--indv', metavar='<int>', required=True, help='full path to individuals list (individuals keep in sscore)' )
@click.option('--output_name', metavar='<str>', required=True, help='output name' )
def filtered_sscore_w_indv(data_set_name, clump_kb, clump_p1, clump_r2, indv, output_name):
    gprs = GPRS()
    gprs.filtered_sscore_w_indv(data_set_name=data_set_name, clump_kb=clump_kb,
                           clump_p1=clump_p1, clump_r2=clump_r2, indv=indv, output_name=output_name)

# main.add_command( test )
main.add_command( build_prs )
main.add_command( clump )
main.add_command( combine_prs )
main.add_command( combine_prs_stat)
main.add_command( geneatlas_filter_data )
main.add_command( generate_plink_bfiles )
main.add_command( generate_plink_bfiles_w_individual_info )
main.add_command( gwas_filter_data )
main.add_command( prs_statistics )
main.add_command( select_clump_snps )
main.add_command( transfer_atcg )
main.add_command( subset_vcf_w_random_sample )
main.add_command( subset_pop )
main.add_command( random_draw_samples_from_fam )
main.add_command( create_new_marker )
main.add_command( filtered_sscore_w_indv )