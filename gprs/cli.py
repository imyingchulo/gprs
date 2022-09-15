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
@click.option('--file/--dir', default=True, help='Whether summary statistics is given as one file, or as a directory with 22 chromosome files with --sumstat')
@click.option('--sumstat', metavar='<str>', required=True, help='Path to one summary statistic file(default), or a directory with 22 chromosome files (use with --dir flag in this case)')
@click.option('--comment', metavar='<str>', default='', help='In summary statistic file(s), indicate the text for lines that should be skipped (for example, "#" for snptest results)')
@click.option('--symbol', metavar='<str>', default='.', help='When giving summary statistics DIRECTORY, indicate the symbol or text after chromosome number in each file, default = "." ')
@click.option('--out', metavar='<str>', required=True, help='Output prefix for 22 processed summary statistics, deposited in sumstat folder')
@click.option('--snpid', metavar='<str>', default=None, help='Column header name for SNP ID in sumstat')
@click.option('--chr', metavar='<str>', default=None, help='Column header name for CHROMOSOME in sumstat')
@click.option('--pos', metavar='<str>', default=None, help='Column header name for POSITION in sumstat')
@click.option('--ea', metavar='<str>', default=None, help='Column header name for EFFECT ALLELE in sumstat')
@click.option('--nea', metavar='<str>', default=None, help='Column header name for NON-EFFECT ALLELE in sumstat')
@click.option('--beta', metavar='<str>', default=None, help='Column header name for BETA(EFFECT SIZE) for EFFECT ALLELE in sumstat')
@click.option('--se', metavar='<str>', default=None, help='Column header name for STANDARD ERROR in sumstat')
@click.option('--pval', metavar='<str>', default=None, help='Column header name for P-VALUE in sumstat')
@click.option('--neff', metavar='<str>', default=None, help='Column header name for EFFECTIVE SAMPLE SIZE in sumstat')
@click.option('--total', metavar='<int>', default=0, help='Total sample size for continuous trait; DO NOT use with --Neff or --case_control')
@click.option('--case_control', metavar='<int>', nargs=2, default=(0,0), help='Case and control sample size for binary trait, separated by a space (order does not matter); DO NOT use with --Neff or --total')
def prepare_sumstat( file, sumstat, comment, symbol, out, snpid, chr, pos, ea, nea,  beta, se, pval, neff, total, case_control ):
   gprs = GPRS()
   gprs.prepare_sumstat(file=file,
                        sumstat=sumstat,
                        comment=comment,
                        out=out,
                        symbol=symbol,
                        snpid=snpid, chr=chr, pos=pos, ea=ea, nea=nea, beta=beta, se=se, pval=pval, neff=neff,
                        total=total, case_control=case_control)

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
@click.option('--merge/--no-merge', default=True, help='Whether to keep or skip merging step; use with --no-merge flag if not using LDPred2 model')
@click.option( '--ref', metavar='<str>', required=True, help='path to population reference panel' )
@click.option( '--sumstat', metavar='<str>', required=True, help='prefix to summary statistics file from perepare_sumstat function.' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output name should be the same as snplist file name' )
@click.option( '--symbol', metavar='<str/int>', required=True, default='.', help='indicate the symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz"' )
@click.option( '--extra_commands', metavar='<str>', default=' ', help='a space to add new functions for generate the plink bfiles' )
def generate_plink_bfiles( merge, ref, sumstat, output_name, symbol,extra_commands):
    gprs = GPRS( ref=ref)
    gprs.generate_plink_bfiles( merge = merge, output_name = output_name, symbol = symbol, sumstat=sumstat, extra_commands=extra_commands )

@click.command()
@click.option( '--plink_bfile_name', metavar='<str>', required=True, help='plink_bfile_name is [output_name] from [chrnb]_[output_name].bim/bed/fam' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output_name remain the same. The clump output: [chrnb]_[output_name]_clumped_snplist.csv' )
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_p2', metavar='<float/scientific notation>', required=True, help='should equals to p1 reduce the snps' )
@click.option( '--clump_r2', metavar='<float>', default=0.1, help='r2 value for clumping, default = 0.1' )
@click.option( '--clump_field', metavar='<str>', default='Pvalue', help='P-value column name, default = Pvalue' )
@click.option( '--sumstat', metavar='<str>', required=True, help='[output_name] from [output_name]_[chrnb].csv in sumstat directory' )
@click.option( '--clump_snp_field', metavar='<str>', default='SNPID', help='SNP ID column name, default = SNPID' )
def clump(sumstat, plink_bfile_name, clump_kb, clump_p1, clump_p2, output_name, clump_r2, clump_field, clump_snp_field):
    gprs = GPRS()
    gprs.clump( sumstat=sumstat,
                plink_bfile_name=plink_bfile_name,
                clump_kb=clump_kb,
                clump_p1=clump_p1,
                clump_p2=clump_p2,
                output_name=output_name,
                clump_r2=clump_r2,
                clump_field=clump_field,
                clump_snp_field=clump_snp_field )
@click.command()
@click.option( '--sumstat', metavar='<str>', required=True, help='[output_name] from [output_name]_[chrnb].csv in sumstat directory' )
@click.option( '--clump_file_name', metavar='<str>', required=True, help='clump_file_name is [output_name] from [chrnb]_[output_name].clump' )
@click.option( '--output_name', metavar='<str>', required=True, help='it is better if the output_name remain the same. output: [chrnb]_[output_name]_clumped_snplist.csv' )
@click.option( '--clump_kb', metavar='<int>', required=True, help='distance(kb) parameter for clumping' )
@click.option( '--clump_p1', metavar='<float/scientific notation>', required=True, help='first set of P-value for clumping' )
@click.option( '--clump_r2', metavar='<float>', required=True, help='r2 value for clumping' )
@click.option('--clumpfolder_name',metavar='<str>', required=True, help='folder name for .clump files')
def select_clump_snps(clump_file_name, sumstat,output_name,clump_kb,clump_p1,clump_r2, clumpfolder_name):
    gprs = GPRS()
    gprs.select_clump_snps( sumstat=sumstat, clump_file_name=clump_file_name, output_name=output_name,
                            clump_kb=clump_kb,clump_p1=clump_p1,clump_r2=clump_r2,clumpfolder_name=clumpfolder_name)

@click.command()
@click.option( '--bfile', metavar='<str>', required=True, help='prefix to training PLINK files including full path')
@click.option( '--LDref', metavar='<str>', default='', help='If using external LD reference, provide directory for a PLINK file with all chromsome merged')
@click.option( '--LDmatrix', metavar='<str>', default='./tmp-data/LD_matrix', help='Path to LD matrix directory')
@click.option( '--sumstat', metavar='<str>', required=True, help='prefix to GWAS by-chromosome summary statistics including full path')
@click.option( '--output_dir', metavar='<str>', required=True, help='directory name to output beta files')
@click.option( '--h2', metavar='<float>', default='', help='heritability estimate')
def ldpred2_train(bfile, ldref, ldmatrix, sumstat, output_dir, h2):
    gprs=GPRS()
    gprs.ldpred2_train(bfile=bfile, 
                        ldref=ldref,
                        ldmatrix=ldmatrix,
                        sumstat=sumstat,
                        output_dir=output_dir,
                        h2=h2)


@click.command()
@click.option( '--beta_dirs', metavar='<str>', required=True, help='list of beta directories to compute PRS with, separated by space')
@click.option( '--out',metavar='<str>', required=True, help='prefix for output .list file')
def beta_list(beta_dirs, out):
    gprs = GPRS()
    gprs.beta_list( beta_dirs=beta_dirs, out=out)

@click.command()
@click.option( '--vcf_dir', metavar='<str>', required=True, help='path to directories containing vcf files')
@click.option( '--beta_dir_list', metavar='<str>', required=True, help='list of beta directories created from beta-list function. If not moved, it is in ./result/prs ')
@click.option( '--slurm_name', metavar='<str>', required=True, help='slurm job name')
@click.option( '--slurm_account', metavar='<str>', required=True, default='chia657_28', help='slurm job account; default="chia657_28" ')
@click.option( '--slurm_time', metavar='<str>', required=True, default='12:00:00', help='slurm job time; default="12:00:00" ')
@click.option( '--memory', metavar='<int>',required=True, default=10, help='slurm job memory in GB; default="10" ')
@click.option( '--symbol', metavar='<str>', required=True, default='.', help='symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz"')
@click.option( '--columns', metavar='<int>', default='1 4 6', help='a column index indicate the [SNPID] [ALLELE] [BETA] position; column nb starts from 1; default="1 4 6"' )
@click.option( '--plink_modifier', metavar='<str>', required=True, default="'no-mean-imputation' 'cols='nmissallele,dosagesum,scoresums", help='plink2 modifier for score function')
@click.option( '--combine', metavar='<str>', required=True, default='T', help='whether to combine scores per chromosomes to generate a final genome-wide PRS (T/F); default="T" ')
@click.option( '--out', metavar='<str>', required=True, default='', help='directory name to output PRS')
def multiple_prs(vcf_dir, beta_dir_list, slurm_name, slurm_account, slurm_time, memory, symbol, columns, plink_modifier, combine, out):
    gprs=GPRS()
    gprs.multiple_prs( vcf_dir=vcf_dir,
                        beta_dir_list=beta_dir_list,
                        slurm_name=slurm_name,
                        slurm_account=slurm_account,
                        slurm_time=slurm_time,
                        memory=memory,
                        symbol=symbol,
                        columns=columns,
                        plink_modifier=plink_modifier,
                        combine=combine,
                        out=out)

@click.command()
@click.option( '--vcf_dir', metavar='<str>', required=True, help='path to vcf files' )
@click.option( '--model', metavar='<str>', required=True, help='model to use to generate PRS')
@click.option( '--beta_dir_list', metavar='<str>', required=True, help='list of beta directories created from beta-list function, used to look up the path for specified PRS model')
@click.option( '--memory', metavar='<int>', help='number of memory use' )
@click.option( '--out', metavar='<str>', required=True, default='', help='directory name to output PRS')
@click.option( '--symbol', metavar='<str/int>', required=True,default='.', help='indicate the symbol or text after chrnb in vcf file, default = "." ; i.e. ALL.chr8.vcf.gz, you can put "." or ".vcf.gz"' )
@click.option( '--columns', metavar='<int>', default='1 2 3', help='a column index indicate the [SNPID] [ALLELE] [BETA] position; column nb starts from 1 ' )
@click.option( '--plink_modifier', metavar='<str>', default='no-mean-imputation', help='no-mean-imputation as default in here, get more info by searching plink2.0 modifier ' )
@click.option( '--combine', metavar='<str>', required=True, default='T', help='whether to combine scores per chromosomes to generate a final genome-wide PRS (T/F); default="T" ')
def build_prs(vcf_dir, model, beta_dir_list, memory, out, symbol, columns, plink_modifier, combine):
    gprs = GPRS()
    gprs.build_prs( vcf_dir=vcf_dir,
                    model=model,
                    beta_dir_list=beta_dir_list,
                    memory=memory,
                    out=out,
                    symbol=symbol,
                    columns=columns,
                    plink_modifier=plink_modifier,
                    combine=combine)

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
@click.option( '--model_name', metavar='<str>', required=True, help='the model name to be used for output')
@click.option( '--r_command', metavar='<str>', required=True, help='use "which R" in linux, and copy the path after --r_command')
@click.option('--binary/--quantitative', default=False, help='whether phenotype is binary or quantitative; default: --quantitative')
@click.option( '--pop_prev', metavar='<str>', default='NA', help='population prevalence for binary trait. Required for binary trait but leave it blank or enter NA for quantitative trait')
@click.option( '--plotroc/--no_plot', metavar='<str>', default=False, help='whether to plot ROC curve for binary trait. Leave it blank for --no_plot for quantitative trait')
def prs_statistics(score_file, pheno_file, model_name, r_command, binary, pop_prev, plotroc):
    gprs = GPRS()
    gprs.prs_statistics( score_file=score_file,
                         pheno_file=pheno_file,
                         model_name=model_name,
                         r_command=r_command,
                         binary=binary,pop_prev=pop_prev,plotroc=plotroc)

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
@click.option('--output_name', metavar='<str>', help='output name for sub-set population bfile')
def generate_plink_bfiles_w_individual_info(popfile_name, bfile_name, output_name):
    gprs = GPRS()
    gprs.generate_plink_bfiles_w_individual_info(output_name=output_name, bfile_name=bfile_name,
                                                 popfile_name=popfile_name)

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

# main.add_command( test )
main.add_command( beta_list )
main.add_command( build_prs )
main.add_command( clump )
main.add_command( combine_prs )
main.add_command( combine_prs_stat)
main.add_command( geneatlas_filter_data )
main.add_command( generate_plink_bfiles )
main.add_command( generate_plink_bfiles_w_individual_info )
main.add_command( gwas_filter_data )
main.add_command( ldpred2_train )
main.add_command( multiple_prs )
main.add_command( prepare_sumstat )
main.add_command( prs_statistics )
main.add_command( select_clump_snps )
main.add_command( transfer_atcg )
main.add_command( subset_vcf_w_random_sample )
main.add_command( subset_pop )

