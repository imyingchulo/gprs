import sys
import os
import glob
import random
import pandas as pd
from pathlib import Path
from subprocess import call
from timeit import default_timer as timer
from collections import defaultdict

from sklearn.utils import column_or_1d


class GPRS(object):
    def __init__(self,
                 ref='',
                 data_dir='',
                 result_dir='./result'):
        """
        At the beginning of gprs package, result and sub-folders will generate automatically when the user first run the package.
        """
        self.ref = ref
        self.data_dir = data_dir
        self.result_dir = Path('{}/{}'.format(os.getcwd(), result_dir)).resolve()
        self.sumstat_dir = '{}/{}'.format(self.result_dir, 'sumstat')
        self.plink_dir = '{}/{}'.format(self.result_dir, 'plink')
        self.pop_dir = '{}/{}'.format(self.result_dir, 'pop')
        self.random_draw_sample_dir = '{}/{}'.format(self.result_dir, 'random_draw_sample')
        self.stat_dir = '{}/{}'.format(self.result_dir, 'stat')
        self.plink_bfiles_dir = '{}/{}'.format(self.plink_dir, 'bfiles')
        self.plink_clump_dir = '{}/{}'.format(self.plink_dir, 'clump')
        self.prs_dir = '{}/{}'.format(self.result_dir, 'prs')
        self.ct_dir = '{}/{}'.format(self.plink_dir, 'ct')
        self.ldpred2_dir = '{}/{}'.format(self.result_dir, 'ldpred2')
        self.setup_dir()

    def setup_dir(self):  # The setup_dir function is automatically create 10 folders
        self.create_result_dir()
        self.create_sumstat_dir()
        self.create_plink_dir()
        self.create_plink_bfiles_dir()
        self.create_plink_clump_dir()
        self.create_prs_dir()
        self.create_ct_dir()
        self.create_stat_dir()
        self.create_random_draw_sample_dir()
        self.create_ldpred2_dir()

    def create_result_dir(self):  # A function to create result folder
        if not os.path.exists(self.result_dir):
            os.mkdir(self.result_dir)

    def create_sumstat_dir(self): # A function to create sumstat folder
        if not os.path.exists(self.sumstat_dir):
            os.mkdir(self.sumstat_dir)

    def create_plink_dir(self):  # A function to create plink folder
        if not os.path.exists(self.plink_dir):
            os.mkdir(self.plink_dir)

    def create_stat_dir(self):  # A function to create stat folder
        if not os.path.exists(self.stat_dir):
            os.mkdir(self.stat_dir)

    def create_plink_bfiles_dir(self):  # A function to create plink bfile folder
        if not os.path.exists(self.plink_bfiles_dir):
            os.mkdir(self.plink_bfiles_dir)

    def create_plink_clump_dir(self):  # A function to create plink clump folder
        if not os.path.exists(self.plink_clump_dir):
            os.mkdir(self.plink_clump_dir)

    def create_prs_dir(self):  # A function to create prs folder
        if not os.path.exists(self.prs_dir):
            os.mkdir(self.prs_dir)

    def create_ct_dir(self):  # A function to create qc'd clump snplists folder
        if not os.path.exists(self.ct_dir):
            os.mkdir(self.ct_dir)

    def create_random_draw_sample_dir(self):
        if not os.path.exists(self.random_draw_sample_dir):
            os.mkdir(self.random_draw_sample_dir)

    def create_ldpred2_dir(self):
        if not os.path.exists(self.ldpred2_dir):
            os.mkdir(self.ldpred2_dir)

    # Unify sumstat format
    def prepare_sumstat(self, file, sumstat, out, symbol='.', comment='',
                        snpid=None, chr=None, pos=None, ea=None, nea=None, beta=None, se=None, pval=None, neff=None,
                        total=0, case_control=(0,0)):
        # dict for column name mapping        
        columns = {
            'SNPID': snpid,
            'CHR': chr ,
            'POS': pos,
            'Effect_Allele': ea,
            'NonEffect_Allele': nea,
            'Beta': beta,
            'SE': se,
            'Pvalue': pval,
            'N_eff': neff
        }       
        col_exist = {value: key for (key, value) in columns.items() if value != None }
        
        # function to unify format
        def unify(df, columns, neff, total, case_control):
            #fill in N_eff if given as number
            if neff == None:
                if total != 0:
                    print('Inserting {} as effective sample size. Make sure it is TOTAL Sample size if you are using continuous trait, or EFFECTIVE Sample size if binary trait\n'.format(total))
                    df['N_eff'] = total
                elif case_control != (0,0):
                    cs, ct = case_control 
                    print('\nUsing {} and {} to calculate effective sample size. Make sure you are using BINARY trait'.format(cs,ct))
                    print('Formula: 4 / (1 / case + 1 / control)\n')
                    df['N_eff'] = int(4/(1/cs + 1/ct))
            # fill in missing column as NA
            for col in columns.keys():
                if col not in list(df):
                    print('WARNING: Header for {} not provided (will be recorded as NA)'.format(col))
                    df[col]='NA'                

            # filter, re-order, unify dtypes columns
            df = df[['SNPID','CHR','POS','Effect_Allele','NonEffect_Allele','Beta','SE','Pvalue','N_eff' ]]
            pd.set_option('mode.chained_assignment', None)
            df[['Effect_Allele','NonEffect_Allele']] = df[['Effect_Allele','NonEffect_Allele']].astype(str).apply(lambda x: x.str.upper())
            df = df.astype( dtype= {'CHR': int, 'POS': int, 'N_eff': int}, errors='raise')  
            return df
        
        # if sumstat given as one file
        if file:
            if not os.path.isfile(sumstat):
                sys.exit('ERROR: {} is not a file. Check help page.\n'.format(sumstat))
            print('Processing one summary statistics file...')
            if len(comment) == 0:         
                df=pd.read_csv(sumstat, delim_whitespace=True)
            else:
                df=pd.read_csv(sumstat, delim_whitespace=True, comment=comment)
            df.rename(columns=col_exist, inplace=True)
            #chr column missing --> look for it in SNP ID
            if chr == None:
                if snpid == None:
                    sys.exit('ERROR: both chromosome and SNP ID header not provided')
                else:
                    print('Looking for chromosome information from SNP ID...Make sure it is in chr:pos:allele1:allele2 format.\n')                    
                    df[['CHR','POS','A1','A2']] =  df['SNPID'].str.split(pat=':',expand=True)
                    if sum( df['CHR'].str.contains("chr") ) == len(df):
                       df['CHR']= df['CHR'].apply(lambda x: x[3:])
                    elif sum( df['CHR'].str.contains("chr") ) != 0:
                        print('ERROR: SNP ID format not consistent.\n')
            df = unify(df, columns, neff, total, case_control)
            #write into 22 sumstat files
            for chrnb in range(1,23):
                outfile="{}/{}_chr{}.csv".format(self.sumstat_dir, out, chrnb)
                df[df.CHR == chrnb].to_csv(outfile, index=False, sep='\t', header=True)

        #sumstat given as directory  
        else:
            if not os.path.isdir(sumstat):
                sys.exit('Error: {} is not a directory. Check help page.\n'.format(sumstat))
            print('Directory for summary statistic files provided...')
            for chrnb in range(1,23):
                chr_symbol='chr{}{}'.format(chrnb,symbol)
                for i in os.listdir(sumstat):
                    if chr_symbol in i:
                        sumstat_chr=i
                        if len(comment) == 0:
                            df = pd.read_csv( "{}/{}".format(sumstat, sumstat_chr), delim_whitespace=True)
                        else:
                            df = pd.read_csv( "{}/{}".format(sumstat, sumstat_chr), delim_whitespace=True, comment=comment)                
                        df.rename(columns=col_exist, inplace=True)

                        #fill in chromosome from filename
                        if chr == None:
                            print('WARNING: chromosome header not provided. Looking for it in filename..')
                            df['CHR']=chrnb
                        df = unify(df, columns, neff, total, case_control)
                        print('Processing chromosome {} done\n'.format(chrnb))
                        outfile="{}/{}_chr{}.csv".format(self.sumstat_dir, out, chrnb)
                        df[df.CHR == chrnb].to_csv(outfile, index=False, sep='\t', header=True)

        print('\nProcessing Done. 22 summary statistics saved in result/sumstat folder!\n')
    
    # Using plink to generate bfiles fam/bim/bed.
    def generate_plink_bfiles(self, merge, sumstat, output_name, symbol='.', extra_commands=" "):
        def run_plink_bfiles():
            visited = set()
            if "{}_{}".format(chrnb, output_name)not in visited:
                print("start to generate {}_{} bfiles".format(chrnb, output_name))
                os.system("plink --vcf {} --extract {} {} --make-bed --out {}".format(vcfinput, snp, extra_commands, output))
            print("{}_{} bfile is ready!".format(chrnb, output_name))
            visited.add("{}_{}".format(chrnb, output_name))

        if any("chr" in file and "{}".format(sumstat) in file for file in os.listdir(self.sumstat_dir)):
            # Generate chr number (chr1-chr22)
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                snp = "{}/{}_{}.csv".format(self.sumstat_dir, sumstat, chrnb)
                output = "{}/{}_{}".format(self.plink_bfiles_dir, chrnb, output_name)
                for i in os.listdir(self.ref):
                        # The if condition here is to exclude chr X, Y, and MT, and make sure all inputs are consistent
                        # ex: The input should be chr1 snps-list and chr1 vcf reference
                        if i.endswith('.vcf.gz') and chrnb != "chrX" and chrnb != "chrY" and chrnb != "chrMT" and "{}{}".format(chrnb, symbol) in i:
                            vcfinput = "{}/{}".format(self.ref, i)
                            print("summary statistics: {}, output: {}, vcfinput:{}\n start to generate bfiles".format(snp, output, vcfinput))
                            run_plink_bfiles()
        else:
            print("ERROR: chromosome information are NOT found in snplists")
        if merge:
            print("Merging 22 plink files...")
            with open("{}/merge.list".format(self.plink_bfiles_dir),'w') as o:
                for chrnb in range(2,23):
                   o.write("{}/chr{}_{}\n".format(self.plink_bfiles_dir, chrnb, output_name))
            os.system("plink --bfile {}/chr1_{} --merge-list {}/merge.list --make-bed --out {}/merged_{}".format(
                                    self.plink_bfiles_dir, output_name, self.plink_bfiles_dir, self.plink_bfiles_dir, output_name ))
            print("Merged file saved!")

    def clump(self, sumstat, plink_bfile_name, output_name, clump_kb, clump_p1, clump_p2, clump_r2='0.1',
              clump_field='Pvalue', clump_snp_field='SNPID'):
        # Create a C+T tag
        output_name_with_conditions = "{}_{}_{}_{}".format(output_name, clump_kb, clump_p1, clump_r2)

        # check the dir exists for output
        print("check the directory")
        if os.path.exists("{}/{}".format(self.plink_clump_dir, output_name_with_conditions)):
            print("{}/{} exists".format(self.plink_clump_dir, output_name_with_conditions))
        else:
            # Create a folder with C+T tag under plink/clump
            print("{}/{} not exists, going to create a folder".format(self.plink_clump_dir, output_name_with_conditions))
            os.mkdir("{}/{}".format(self.plink_clump_dir, output_name_with_conditions))

        def run_plink():
            visited = set()
            if sumstat_files not in visited:
                print("start {}/{} clumping".format(self.sumstat_dir, sumstat_files))
                os.system("plink --bfile {} --clump {}/{} --clump-p1 {} --clump-p2 {} --clump-r2 {} --clump-kb {} --clump-field {} --clump-snp-field {} --out {} ".format(
                        plinkinput,
                        self.sumstat_dir, sumstat_files,
                        clump_p1, clump_p2, clump_r2, clump_kb, 
                        clump_field, clump_snp_field, output))
            visited.add("{}".format(sumstat_files))
            print("finished {}/{} clumping".format(self.sumstat_dir, sumstat_files))

        if any("chr" in file and "{}".format(sumstat) in file for file in os.listdir(self.sumstat_dir)):
            # Generate chr number (chr1-chr22)
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                sumstat_files = "{}_{}.csv".format(sumstat, chrnb)
                output = "{}/{}/{}_{}".format(self.plink_clump_dir, output_name_with_conditions, chrnb, output_name_with_conditions)
                plinkinput = "{}/{}_{}.bim".format(self.plink_bfiles_dir, chrnb, plink_bfile_name).split(".")[0]
                if os.path.exists("{}.bim".format(plinkinput)):
                    print("sumstat_files:{} \noutput:{} \nplinkinput:{}".format(sumstat_files, output, plinkinput))
                    run_plink()
                else:
                    print("Warning: {}.bim not found. Moving on to next file".format(plinkinput))
        else:
             print("ERROR: cannot file summary statistic files")
        #     for nb in range(1, 23):
        #         chrnb = "chr{}".format(nb)
        #         qc_files = "{}.QC.csv".format(qc_file_name)
        #         output = "{}/{}/{}_{}".format(self.plink_clump_dir, output_name_with_conditions, chrnb, output_name_with_conditions)
        #         plinkinput = "{}/{}_{}.bim".format(self.plink_bfiles_dir, chrnb, plink_bfile_name).split(".")[0]
        #         if os.path.exists("{}.bim".format(plinkinput)):
        #             print("qc_files:{} \noutput:{} \nplinkinput:{}".format(qc_files, output, plinkinput))
        #             run_plink()
        #         else:
        #             print("{}.bim not found. Move to next file".format(plinkinput))
        print("All chromosome clumping finished!")

    def select_clump_snps(self, sumstat, clump_file_name, clumpfolder_name, output_name, clump_kb, clump_p1, clump_r2):
        # Create a C+T tag
        clump_conditions = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2)

        # check dir
        if os.path.exists("{}/{}_{}".format(self.plink_clump_dir, clumpfolder_name, clump_conditions)):
            pass
        else:
            os.mkdir("{}/{}_{}".format(self.plink_clump_dir, clumpfolder_name, clump_conditions))

        if os.path.exists("{}/{}_{}".format(self.ct_dir, clumpfolder_name, clump_conditions)):
            pass
        else:
            os.mkdir("{}/{}_{}".format(self.ct_dir, clumpfolder_name, clump_conditions))

        # STEP1: generate [filename]_clumped_snplist.csv
        def generate_clumped_snplist():
            visited = set()
            if "{}_{}_{}".format(chrnb, clump_file_name, clump_conditions) not in visited:
                print("starting to generate {}_{}_{}_clumped_snplist.csv".format(chrnb, clump_file_name, clump_conditions))
                new_df = df.loc[(df['CHR'] == nb)]
                new_df_2 = new_df[['CHR', 'SNP']]
                new_df_2.to_csv("{}/{}_{}/{}_{}_{}_clumped_snplist.csv".format(self.plink_clump_dir, clumpfolder_name,
                                                                               clump_conditions, chrnb, clump_file_name,
                                                                               clump_conditions), sep=' ', index=False, header=True)
                visited.add("{}_{}_{}".format(chrnb, clump_file_name, clump_conditions))
                print("{}_{}_{} clumped snplist created".format(chrnb, clump_file_name, clump_conditions))

        for nb in range(1, 23):
            chrnb = "chr{}".format(nb)
            if os.path.exists("{}/{}_{}/{}_{}_{}.clumped".format(self.plink_clump_dir, clumpfolder_name, clump_conditions, chrnb, clump_file_name, clump_conditions)):
                df = pd.read_csv("{}/{}_{}/{}_{}_{}.clumped".format(self.plink_clump_dir, clumpfolder_name, clump_conditions, chrnb, clump_file_name, clump_conditions), delim_whitespace=True)
                generate_clumped_snplist()
            else:
                print("{}/{}_{}/{}_{}_{}.clumped not found. Move to next .clumped file".format(self.plink_clump_dir, clumpfolder_name, clump_conditions, chrnb, clump_file_name, clump_conditions))
        print("STEP 1 finished. Move to STEP2")

        # STEP2: generate snplist with qc information (Allele Beta SE Pvalue)
        def generate_qc_snplist():
            visited = set()
            if clump_snp_file not in visited:
                    print("start extracting {} and {}".format(sumstat_files, clump_snp_file))
                    clump_snp = pd.read_csv(clump_snp_file, delim_whitespace=True)
                    clump_snp.rename(columns={'SNP': 'SNPID'}, inplace=True)
                    qc_snp = pd.read_csv("{}".format(sumstat_files), delim_whitespace=True)
                    newsnplist = qc_snp[qc_snp["SNPID"].isin(clump_snp["SNPID"])]
                    newsnplist.to_csv("{}/{}_{}/{}.weight".format(self.ct_dir, clumpfolder_name, clump_conditions, output), sep=' ', index=False, header=True)
                    print("{}.weight created".format(output))
            else:
                pass
            visited.add("{}".format(clump_snp_file))

        if any("chr" in file and "{}".format(sumstat) in file for file in os.listdir(self.sumstat_dir)):
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                clump_snp_file = ("{}/{}_{}/{}_{}_{}_clumped_snplist.csv".format(self.plink_clump_dir,
                                                                                 clumpfolder_name,
                                                                                 clump_conditions,
                                                                                 chrnb, clump_file_name,
                                                                                 clump_conditions))
                sumstat_files = "{}/{}_{}.csv".format(self.sumstat_dir, sumstat, chrnb)
                output = "{}_{}_{}".format(chrnb, output_name, clump_conditions)
                if os.path.exists("{}".format(clump_snp_file)):
                    print("clump_snp_file:{} \noutput:{} \nsumstat_files:{}".format(clump_snp_file, output, sumstat_files))
                    generate_qc_snplist()
                else:
                    print("{} not found skip".format(clump_snp_file))
        else:
            print("ERROR: chr information are not found")
            # for nb in range(1, 23):
            #     chrnb = "chr{}".format(nb)
            #     clump_snp_file = ("{}/{}_{}/{}_{}_{}_clumped_snplist.csv".format(self.plink_clump_dir,
            #                                                                      clumpfolder_name,
            #                                                                      clump_conditions,
            #                                                                      chrnb, clump_file_name,
            #                                                                      clump_conditions))
            #     output = "{}_{}_{}".format(chrnb, output_name, clump_conditions)
            #     qc_files = "{}/{}.QC.csv".format(self.qc_dir, qc_file_name)
            #     if os.path.exists("{}".format(clump_snp_file)):
            #         print("clump_snp_file:{} \noutput:{} \nqc_files:{}".format(clump_snp_file, output, qc_files))
            #         generate_qc_snplist()
            #     else:
            #         print("{} not found skip".format(clump_snp_file))
        print("All jobs are completed")

    def ldpred2_train(self, bfile, sumstat, output_dir, h2='', ldref='', ldmatrix='./tmp-data/LD_matrix'):
        command="Rscript --vanilla ./gprs/ldpred2.R --train {} --sumstat {} --output_dir {}/{}".format(bfile, sumstat,self.ldpred2_dir, output_dir)                                                                                        
        if len(ldref) > 0 :
            command += " --LDref {}".format(ldref)
        if ldmatrix != './tmp-data/LD_matrix':
            command += " --LDmatrix {}".format(ldmatrix)
        if len(h2) > 0 :
            command += " --h2 {}".format(h2)

        call(command, shell=True)

    #make beta list for multiple_prs function.
    def beta_list(self, beta_dirs, out):
        beta_dirs = beta_dirs.split()
        print('Iterating {} beta directoreis..'.format( len(beta_dirs)))
        allmodels={}
        for beta_dir in beta_dirs:
            #get list of directories in given path
            models=[ x for x in os.listdir(beta_dir) if os.path.isdir( os.path.join(beta_dir, x)) ]
            print('- {} models found in {}\n\n'.format( len(models), beta_dir))
            for model in sorted(models):
                allmodels[model] = '{}/{}'.format(beta_dir, model)
        #write .list
        with open('{}/{}.list'.format(self.prs_dir, out), 'w') as o:
            for x in allmodels.keys():
                o.write('{}\t{}\n'.format(x, allmodels[x]))
        print('File {}.list saved in ./result/prs!\n'.format(out))

    def multiple_prs(self, vcf_dir, beta_dir_list,
                slurm_name, slurm_account='chia657_28', slurm_time='12:00:00', memory=10,
                symbol='.',
                columns='1 4 6', plink_modifier='no-mean-imputation cols=nmissallele,dosagesum,scoresums',
                combine='T', out=''):
    # read in list of beta directories and path dictionary
        print('Reading list of beta directories..{}'.format(beta_dir_list))
        beta_list_file = pd.read_csv(beta_dir_list, header=None, sep='\t', index_col=0, squeeze=True)
        beta_list = beta_list_file.T.to_dict()
        print('{} models found'.format( len(beta_list)))

        # create output folder
        if out == '':
            out = self.prs_dir
            print('Warning: output directory not defined. Will be saved in {}'.format(self.prs_dir))
        else:
            out = "{}/{}".format( self.prs_dir, out )
            os.makedirs(out, exist_ok=True)
            print('Output will be saved in {}'.format(out))

        # compile args for build_prs
        args = "--beta_dir_list {} --memory {} --symbol {} --columns '{}' --plink_modifier '{}' --combine {} --out {}".format(
                beta_dir_list, memory*1000, symbol, columns, plink_modifier, combine, out)

        # write slurm script
        command = """#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --account={}
#SBATCH --time={}
#SBATCH --mem-per-cpu={}gb
#SBATCH --job-name={}
#SBATCH --output=slurm.%x.%j.out
#SBATCH --error=slurm.%x.%j.err
#SBATCH --array=0-{}\n
module load usc
module load plink2
arr=( {} )\n

source ./venv/bin/activate\n

gprs build-prs --vcf_dir {} --model""".format(
                    slurm_account, slurm_time, memory, slurm_name, len(beta_list)-1, 
                    ' '.join(beta_list.keys()), 
                    vcf_dir)
        command+= ' ${arr[$SLURM_ARRAY_TASK_ID]} '
        command+= args+'\n'

        with open('build-prs.sh','w') as o:
            o.writelines(command)
        os.system('sbatch build-prs.sh')

    def build_prs(self, vcf_dir, model, beta_dir_list, memory, out,
                    symbol='.',
                    columns='1 4 6', plink_modifier="no-mean-imputation cols=nmissallele,dosagesum,scoresums",
                    combine='T'):            
        # create output folder per model with same name
        os.mkdir("{}/{}".format(out, model))
        # read in list of beta directories and path dictionary
        beta_list_file = pd.read_csv(beta_dir_list, header=None, sep='\t', index_col=0, squeeze=True)
        beta_list = beta_list_file.T.to_dict()
        print('Building PRS for "{}" model in {}...'.format(model, beta_list[model]))
        # master score dataframe
        all=pd.DataFrame(columns=['#IID'])

        for nb in range(1, 23):
            chrnb = "chr{}".format(nb)
            # Define and check number of beta file
            beta_file = list( filter( lambda x: x.startswith("{}{}".format(chrnb, "_")) and x.endswith(".weight"), os.listdir(beta_list[model])))
            if len(beta_file) < 1 :
                print("{} beta file in {} not found. skip".format(chrnb, model))
            elif len(beta_file) > 1 :
                raise Exception("Multiple {} beta files in {}. Skip".format(chrnb, model))
            else:
                # only if there's one beta file
                beta_file=beta_file[0]
                # Define vcf file
                for vcf_file in os.listdir(vcf_dir):
                    if vcf_file.endswith('.vcf.gz') and "{}{}".format(chrnb, symbol) in vcf_file:
                        os.system("plink2 --vcf {}/{} dosage=DS --score {}/{} {} {} --memory {} --out {}/{}/{}_{}".format(
                                                                            vcf_dir, vcf_file,
                                                                            beta_list[model], beta_file, columns, plink_modifier,
                                                                            memory,
                                                                            out,model, chrnb, model))
                        print("{}_{}.sscore saved in {}/{}".format(chrnb, model, out, model))
                        file = pd.read_csv('{}/{}/{}_{}.sscore'.format(out, model, chrnb, model), sep='\t')
                        file.rename(columns={'SCORE1_SUM':'SCORE_{}'.format(chrnb), 'NMISS_ALLELE_CT':'ALLELE_CT_{}'.format(chrnb)}, inplace=True)
                        if len(all) == 0:
                            all=all.merge(file[['#IID','SCORE_{}'.format(chrnb), 'ALLELE_CT_{}'.format(chrnb)]], how='right')
                        else:
                            all=all.merge(file[['#IID','SCORE_{}'.format(chrnb), 'ALLELE_CT_{}'.format(chrnb)]], on='#IID')
                
        if combine == 'T':
            print('\nOption to combine scores per chromsome is ON')
            all['SCORE_SUM']=all.filter(like='SCORE').sum(axis=1)
            all['TOTAL_ALLELE_CT']=all.filter(like='ALLELE_CT').sum(axis=1)
            print('Summary for scores...')
            print(all[['SCORE_SUM']].describe())
            all[['#IID', 'SCORE_SUM', 'TOTAL_ALLELE_CT']].to_csv('{}/{}.sscore'.format(out, model), index=False, sep='\t')
            print('Done! Combined score for "{}" model saved'.format(model))

    # Calculate the PRS statistical results and output the statistics summary
    def prs_stat(self, score, pheno, data, model, binary, pop_prev, plotroc, r):
        if os.path.exists(score):
            if binary:
                family = 'binary'
            else:
                family = 'quantitative'
            
            if plotroc:
                plotroc = 'plotroc'
            else:
                plotroc = 'no_plot'

            if not os.path.exists("{}/{}".format(self.stat_dir, data)):
                os.mkdir("{}/{}".format(self.stat_dir, data))
            else:
                print("Warning: statistics directory for {} already exists. Result will be overwritten with same model name. \n".format(data))
            # The R script is written by Soyoung Jeon
            call("{0}script --vanilla ./gprs/prs_stat.R {1} {2} {3} {4} {5} {6} {7}/{8}/{1}".format(r, model, score, pheno,
                                                                            family, pop_prev, plotroc, 
                                                                            self.stat_dir, data), shell=True)
        else:
            print("{} not found. Please check the sscore again".format(score))


    # In combine_prs_stat function is to combine PRS statistical results as one file
    def combine_stat(self, data):
        target_files = [i for i in glob.glob("{}/{}/*.stat".format(self.stat_dir, data))]
        print("Combining : \n\t{}".format("\n\t".join(target_files)))
        # Combine all the csv files
        combined_csv = pd.concat([pd.read_csv(f) for f in target_files])
        # Export to csv
        combined_csv.to_csv("{}/{}_combined.stat".format(self.stat_dir, data), index=False, header=True, sep='\t', quoting=None)
        print("Combining Done. Result saved in {}/{}_combined.stat.".format(self.stat_dir,data))
