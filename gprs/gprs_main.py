import glob
import os
import random
from collections import defaultdict
from pathlib import Path
from subprocess import call

import pandas as pd


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
        self.plink_dir = '{}/{}'.format(self.result_dir, 'plink')
        self.pop_dir = '{}/{}'.format(self.result_dir, 'pop')
        self.random_draw_sample_dir = '{}/{}'.format(self.result_dir, 'random_draw_sample')
        self.stat_dir = '{}/{}'.format(self.result_dir, 'stat')
        self.plink_bfiles_dir = '{}/{}'.format(self.plink_dir, 'bfiles')
        self.plink_clump_dir = '{}/{}'.format(self.plink_dir, 'clump')
        self.prs_dir = '{}/{}'.format(self.plink_dir, 'prs')
        self.qc_dir = '{}/{}'.format(self.result_dir, 'qc')
        self.snplists_dir = '{}/{}'.format(self.result_dir, 'snplists')
        self.qc_clump_snpslist_dir = '{}/{}'.format(self.plink_dir, 'qc_and_clump_snpslist')
        self.setup_dir()

    def setup_dir(self):  # The setup_dir function is automatically create 10 folders
        self.create_result_dir()
        self.create_plink_dir()
        self.create_plink_bfiles_dir()
        self.create_plink_clump_dir()
        self.create_prs_dir()
        self.create_qc_dir()
        self.create_snplists_dir()
        self.create_qc_clump_snpslist_dir()
        self.create_stat_dir()
        self.create_pop_dir()
        self.create_random_draw_sample_dir()

    def create_result_dir(self):  # A function to create result folder
        if not os.path.exists(self.result_dir):
            os.mkdir(self.result_dir)

    def create_plink_dir(self):  # A function to create plink folder
        if not os.path.exists(self.plink_dir):
            os.mkdir(self.plink_dir)

    def create_pop_dir(self):  # A function to create pop folder
        if not os.path.exists(self.pop_dir):
            os.mkdir(self.pop_dir)

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

    def create_qc_dir(self):  # A function to create qc folder
        if not os.path.exists(self.qc_dir):
            os.mkdir(self.qc_dir)

    def create_snplists_dir(self):  # A function to create snplists folder
        if not os.path.exists(self.snplists_dir):
            os.mkdir(self.snplists_dir)

    def create_qc_clump_snpslist_dir(self):  # A function to create qc'd clump snplists folder
        if not os.path.exists(self.qc_clump_snpslist_dir):
            os.mkdir(self.qc_clump_snpslist_dir)

    def create_random_draw_sample_dir(self):
        if not os.path.exists(self.random_draw_sample_dir):
            os.mkdir(self.random_draw_sample_dir)

    # Using plink to generate bfiles fam/bim/bed.
    def generate_plink_bfiles(self, snplist_name, output_name, symbol='.', extra_commands=" "):
        def run_plink_bfiles():
            visited = set()
            if "{}_{}".format(chrnb, output_name)not in visited:
                print("strat to generate {}_{} bfiles".format(chrnb, output_name))
                os.system("plink --vcf {} --extract {} {} --make-bed --out {}".format(vcfinput, snp, extra_commands, output))
            print("{}_{} bfile is ready!".format(chrnb, output_name))
            visited.add("{}_{}".format(chrnb, output_name))

        if any("chr" in file and "{}".format(snplist_name) in file for file in os.listdir(self.snplists_dir)):
            print("chromosome information are found in snplists")
            # Generate chr number (chr1-chr22)
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                snp = "{}/{}_{}.csv".format(self.snplists_dir, chrnb, snplist_name)
                output = "{}/{}_{}".format(self.plink_bfiles_dir, chrnb, output_name)
                for i in os.listdir(self.ref):
                        # The if condition here is to exclude chr X, Y, and MT, and make sure all inputs are consistent
                        # ex: The input should be chr1 snps-list and chr1 vcf reference
                        if i.endswith('.vcf.gz') and chrnb != "chrX" and chrnb != "chrY" and chrnb != "chrMT" and "{}{}".format(chrnb, symbol) in i:
                            vcfinput = "{}/{}".format(self.ref, i)
                            print("snplist: {}, output: {}, vcfinput:{}\n start to generate bfiles".format(snp, output, vcfinput))
                            run_plink_bfiles()
        else:
            print("chromosome information are NOT found in snplists")
            snp = "{}/{}.csv".format(self.snplists_dir, snplist_name)
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                output = "{}/{}_{}".format(self.plink_bfiles_dir, chrnb, output_name)
                for i in os.listdir(self.ref):
                    if i.endswith('.vcf.gz') and chrnb != "chrX" and chrnb != "chrY" and chrnb != "chrMT" and "{}{}".format(chrnb, symbol) in i:
                        vcfinput = "{}/{}".format(self.ref, i)
                        print("snplist: {}, output: {}, vcfinput:{}\n start to generate bfiles".format(snp, output, vcfinput))
                        run_plink_bfiles()
        print("all jobs completed!")

    def clump(self, qc_file_name, plink_bfile_name, output_name, clump_kb, clump_p1, clump_p2, clump_r2='0.1',
              clump_field='Pvalue', clump_snp_field='SNPID'):
        # Create a C+T tag
        output_name_with_conditions = "{}_{}_{}_{}".format(output_name, clump_kb, clump_p1, clump_r2)

        # check the dir exists
        print("check the directory")
        if os.path.exists("{}/{}".format(self.plink_clump_dir, output_name_with_conditions)):
            print("{}/{} exists".format(self.plink_clump_dir, output_name_with_conditions))
        else:
            # Create a folder with C+T tag under plink/clump
            print("{}/{} not exists, going to create a folder".format(self.plink_clump_dir, output_name_with_conditions))
            os.mkdir("{}/{}".format(self.plink_clump_dir, output_name_with_conditions))

        def run_plink():
            visited = set()
            if qc_files not in visited:
                print("start {}/{} clumping".format(self.qc_dir, qc_files))
                os.system("plink --bfile {} --clump {}/{} --clump-p1 {} --clump-p2 {} --clump-r2 {} --clump-kb {} --clump-field {} --clump-snp-field {} --out {} ".format(
                        plinkinput,
                        self.qc_dir, qc_files,
                        clump_p1, clump_p2, clump_r2,
                        clump_kb, clump_field, clump_snp_field, output))
            visited.add("{}".format(qc_files))
            print("finished {}/{} clumping".format(self.qc_dir, qc_files))

        if any("chr" in file and "{}".format(qc_file_name) in file for file in os.listdir(self.qc_dir)):
            print("chromosome information are found in .QC.csv")
            # Generate chr number (chr1-chr22)
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                qc_files = "{}_{}.QC.csv".format(chrnb, qc_file_name)
                output = "{}/{}/{}_{}".format(self.plink_clump_dir, output_name_with_conditions, chrnb, output_name_with_conditions)
                plinkinput = "{}/{}_{}.bim".format(self.plink_bfiles_dir, chrnb, plink_bfile_name).split(".")[0]
                if os.path.exists("{}.bim".format(plinkinput)):
                    print("qc_files:{} \noutput:{} \nplinkinput:{}".format(qc_files, output, plinkinput))
                    run_plink()
                else:
                    print("{}.bim not found. Move to next file".format(plinkinput))
        else:
            print("chromosome information are not found in .QC.csv")
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                qc_files = "{}.QC.csv".format(qc_file_name)
                output = "{}/{}/{}_{}".format(self.plink_clump_dir, output_name_with_conditions, chrnb, output_name_with_conditions)
                plinkinput = "{}/{}_{}.bim".format(self.plink_bfiles_dir, chrnb, plink_bfile_name).split(".")[0]
                if os.path.exists("{}.bim".format(plinkinput)):
                    print("qc_files:{} \noutput:{} \nplinkinput:{}".format(qc_files, output, plinkinput))
                    run_plink()
                else:
                    print("{}.bim not found. Move to next file".format(plinkinput))
        print("All chromosome clumping finished!")

    def select_clump_snps(self, qc_file_name, clump_file_name, clumpfolder_name, output_name, clump_kb, clump_p1, clump_r2):
        # Create a C+T tag
        clump_conditions = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2)

        # check dir
        if os.path.exists("{}/{}_{}".format(self.plink_clump_dir, clumpfolder_name, clump_conditions)):
            pass
        else:
            os.mkdir("{}/{}_{}".format(self.plink_clump_dir, clumpfolder_name, clump_conditions))

        if os.path.exists("{}/{}_{}".format(self.qc_clump_snpslist_dir, clumpfolder_name, clump_conditions)):
            pass
        else:
            os.mkdir("{}/{}_{}".format(self.qc_clump_snpslist_dir, clumpfolder_name, clump_conditions))

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
                    print("start extracting {} and {}".format(qc_files, clump_snp_file))
                    clump_snp = pd.read_csv(clump_snp_file, delim_whitespace=True)
                    clump_snp.rename(columns={'SNP': 'SNPID'}, inplace=True)
                    qc_snp = pd.read_csv("{}".format(qc_files), delim_whitespace=True)
                    newsnplist = qc_snp[qc_snp["SNPID"].isin(clump_snp["SNPID"])]
                    newsnplist.to_csv("{}/{}_{}/{}.qc_clump_snpslist.csv".format(self.qc_clump_snpslist_dir, clumpfolder_name, clump_conditions, output), sep=' ', index=False, header=True)
                    print("{}.qc_clump_snpslist.csv created".format(output))
            else:
                pass
            visited.add("{}".format(clump_snp_file))

        if any("chr" in file and "{}".format(qc_file_name) in file for file in os.listdir(self.qc_dir)):
            print("chr information are found")
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                clump_snp_file = ("{}/{}_{}/{}_{}_{}_clumped_snplist.csv".format(self.plink_clump_dir,
                                                                                 clumpfolder_name,
                                                                                 clump_conditions,
                                                                                 chrnb, clump_file_name,
                                                                                 clump_conditions))
                qc_files = "{}/{}_{}.QC.csv".format(self.qc_dir, chrnb, qc_file_name)
                output = "{}_{}_{}".format(chrnb, output_name, clump_conditions)
                if os.path.exists("{}".format(clump_snp_file)):
                    print("clump_snp_file:{} \noutput:{} \nqc_files:{}".format(clump_snp_file, output, qc_files))
                    generate_qc_snplist()
                else:
                    print("{} not found skip".format(clump_snp_file))
        else:
            print("chr information are not found")
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                clump_snp_file = ("{}/{}_{}/{}_{}_{}_clumped_snplist.csv".format(self.plink_clump_dir,
                                                                                 clumpfolder_name,
                                                                                 clump_conditions,
                                                                                 chrnb, clump_file_name,
                                                                                 clump_conditions))
                output = "{}_{}_{}".format(chrnb, output_name, clump_conditions)
                qc_files = "{}/{}.QC.csv".format(self.qc_dir, qc_file_name)
                if os.path.exists("{}".format(clump_snp_file)):
                    print("clump_snp_file:{} \noutput:{} \nqc_files:{}".format(clump_snp_file, output, qc_files))
                    generate_qc_snplist()
                else:
                    print("{} not found skip".format(clump_snp_file))
        print("All jobs are completed")

    def build_prs(self, vcf_input, output_name, qc_clump_snplist_foldername, memory, clump_kb, clump_p1, clump_r2, symbol='.', columns='1 2 3', plink_modifier='no-mean-imputation'):
        # Create a C+T tag
        clump_conditions = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2)

        # Check the folder exists or not, if not create the folder
        if os.path.exists("{}/{}_{}".format(self.prs_dir, output_name, clump_conditions)):
            print("{}/{}_{} exists".format(self.prs_dir, output_name, clump_conditions))
            pass
        else:
            print("{}/{}_{} not exists, going to create one".format(self.prs_dir, output_name, clump_conditions))
            os.mkdir("{}/{}_{}".format(self.prs_dir, output_name, clump_conditions))

        visited = set()
        for nb in range(1, 23):
            chrnb = "chr{}".format(nb)
            for vcf_file in os.listdir(vcf_input):
                # Define the input files (vcf and qc files)
                if vcf_file.endswith('.vcf.gz') and chrnb != "chrY" and chrnb != "chrX" and chrnb != "wgs" and "{}{}".format(chrnb, symbol) in vcf_file:
                    qc_file = "{}/{}_{}/{}_{}_{}.qc_clump_snpslist.csv".format(self.qc_clump_snpslist_dir,
                                                                               qc_clump_snplist_foldername, clump_conditions,
                                                                               chrnb, qc_clump_snplist_foldername, clump_conditions)
                    if os.path.exists(qc_file) and "{}".format(qc_file) not in visited:
                        os.system("plink2 --vcf {}/{} dosage=DS --score {} {} '{}' --memory {} --out {}/{}_{}/{}_{}_{}".format(
                                                                                vcf_input, vcf_file,
                                                                                qc_file, columns, plink_modifier,
                                                                                memory,
                                                                                self.prs_dir,
                                                                                output_name, clump_conditions,
                                                                                chrnb, output_name, clump_conditions))
                        print("{}_{}_{}.sscore completed!".format(chrnb, output_name, clump_conditions))
                    else:
                        print("{} not found. skip".format(qc_file))
                    visited.add("{}".format(qc_file))
        print("ALL work are complete!")

    def combine_prs(self, filename, clump_kb, clump_p1, clump_r2):
        # Create a C+T tag
        clump_conditions = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2)
        data_dir = "{}/{}_{}".format(self.prs_dir,filename, clump_conditions)
        table = defaultdict(list)

        if os.path.exists(data_dir):
            for i in os.listdir(data_dir):
                if i.endswith(".sscore"):
                    if os.path.exists("{}/{}".format(data_dir, i)):
                        print("start to extract the id, allele nb and score from {}/{}".format(data_dir, i))
                        with open("{}/{}".format(data_dir, i), 'r') as sscorein:
                            for index, rows in enumerate(sscorein.readlines()):
                                # use split to separate each row, and access the data by using row[0], row[1]...
                                row = rows.split()
                                # skip the header
                                if index > 0:
                                    key = "{}".format(row[0])
                                    # append data with select columns
                                    nb = float(row[1])
                                    score = float(row[3]) * nb
                                    table[key].append((nb, score))
                        print("finished to extract the id, allele nb and score from {}/{}".format(data_dir, i))
                    else:
                        print("{}/{} not found. Skip".format(data_dir, i))
        else:
            print("{} did not found".format(data_dir))

        result = ['id\tALLELE_CT\tSCORE_SUM']
        # take the value from dictionary
        for key, value in table.items():
            unzip_value = list(zip(*value))
            allele_ct = sum(list(unzip_value[0]))
            score_sum = sum(list(unzip_value[1]))
            result.append("{}\t{}\t{}".format(key, allele_ct, score_sum))
        combined_prs = '\n'.join(result)

        with open("{}/{}_{}_combined.sscore".format(self.prs_dir, filename, clump_conditions), 'w') as fout:
            fout.write(combined_prs)
        # call("rm -rf {}/{}_{}".format(self.prs_dir, filename, clump_conditions), shell=True)
        print('Combined all the sscore file. Original sscore files deleted')

    # Calculate the PRS statistical results and output the statistics summary
    def prs_statistics(self, score_file, pheno_file, output_name, data_set_name, prs_stats_R, r_command, clump_kb, clump_p1, clump_r2):
        filter_condition = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2)

        # Sum the SNPs number from .qc_clump_snpslist.csv file
        lines = []
        for i in range(1, 23):
            file = "{}/{}_{}/chr{}_{}_{}.qc_clump_snpslist.csv".format(self.qc_clump_snpslist_dir,
                                                                           data_set_name, filter_condition, i,
                                                                           data_set_name, filter_condition)
            if os.path.exists("{}".format(file)):
                print("sum th SNPs from {}".format(file))
                df = pd.read_csv("{}".format(file))
                index = df.index
                number_of_rows = len(index)
                lines.append(number_of_rows)
                print("{} SNPs nb: {}".format(file, number_of_rows))
            else:
                print("{} not found. Skip.".format(file))

        if os.path.exists(score_file):
            # The R script is written by Soyoung Jeon
            # USAGE Rscript --vanilla prs_stats_quantitative_phenotype.R [score file] [pheno file] [target pop for OR] [ref pop for OR] [graph pdf name]
            call("{0} --vanilla {1} {2} {3} {4} {5} {6} {8}/{7}".format(r_command, prs_stats_R, score_file, pheno_file,
                                                                        data_set_name, filter_condition, sum(lines),
                                                                        output_name, self.stat_dir), shell=True)

            # Read the statistics result and reformat it
            # Reformat: separator = tab, change float into scientific notation
            stat_data = pd.read_csv("{}/{}_{}_stat.txt".format(self.stat_dir, output_name, filter_condition))
            if os.path.exists("{}/{}_{}_stat.txt".format(self.stat_dir, output_name, filter_condition)):
                stat_data.to_csv("{}/{}_{}_stat.txt".format(self.stat_dir, output_name, filter_condition),
                    index=False,
                    header=True,
                    sep='\t',
                    float_format='%.2E')
            else:
                print("{}/{}_{}_stat.txt Not Found".format(self.stat_dir, output_name, filter_condition))
        else:
            print("{} not found. Please check the sscore again".format(score_file))

    # In combine_prs_stat function is to combine PRS statistical results as one file
    def combine_prs_stat(self, data_set_name, clump_kb, clump_p1, clump_r2):
        filter_condition = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2)
        extension = "_stat.txt"
        target_files = [i for i in glob.glob("{}/{}*{}".format(self.stat_dir, data_set_name, extension))]
        # Combine all the csv files
        combined_csv = pd.concat([pd.read_csv(f) for f in target_files])
        # Export to csv
        combined_csv.to_csv("{}/combined_{}.txt".format(self.stat_dir, data_set_name), index=False, header=True, sep=' ')

    # Optional functions
    def subset_pop(self, input_data, column_name, pop_info, output_name):
        file = pd.read_csv("{}".format(input_data), delim_whitespace=True)
        subset_file = file.loc[(file["{}".format(column_name)] == "{}".format(pop_info))]
        subset_file.to_csv("{}".format(output_name), sep='\t', index=False, header=True)

    def generate_plink_bfiles_w_individual_info(self, popfile_name, bfile_name, plink_command, output_name):
        for nb in range(1, 23):
            chrnb = "chr{}".format(nb)
            for j in os.listdir(self.pop_dir):
                if "{}.txt".format(popfile_name) in j:
                    for i in os.listdir(self.plink_bfiles_dir):
                        if "{}_{}".format(chrnb, bfile_name) in i:
                            bfile = i.split(".")[0]
                            os.system("plink --bfile {}/{} {} {}/{} --make-bed --out {}/{}_{}".format(
                                self.plink_bfiles_dir, bfile,
                                plink_command,
                                self.pop_dir, j,
                                self.plink_bfiles_dir,
                                chrnb, output_name))
                        print("with chr: {}_{} is finished!".format(chrnb, bfile_name))
        print("all jobs completed!")

    # Transfer a t c g into capital A T C G
    def transfer_atcg(self, qc_file_name):
        for nb in range(1, 23):
            chrnb = "chr{}".format(nb)
            df = pd.read_csv("{}/{}_{}.QC.csv".format(self.qc_dir, chrnb, qc_file_name), sep=' ')
            df.loc[:, 'Allele'] = df['Allele'].apply({'a': 'A', 't': 'T', 'c': 'C', 'g': 'G'}.get)
            df.to_csv("{}/{}_{}.QC.csv".format(self.qc_dir, chrnb, qc_file_name), sep=' ', index=False, header=True)
        print("transfer completed!")

    def subset_vcf_w_random_sample(self, fam_dir, fam_filename, samplesize, vcf_input, symbol):
        def random_draw_samples():
            list = []
            print("start to subset the samples")
            with open(fam_input, 'r') as fin:
                lines = fin.readlines()
                for index, cols in enumerate(lines):
                    col = cols.split(" ")
                    list.append("{}".format(col[1]))

            with open("{}/{}.txt".format(self.pop_dir, fam_name), 'w') as fin:
                random_sample = random.sample(list, int(samplesize))
                final_list = '\n'.join(random_sample)
                fin.write(final_list)
                print("random sample list created {}.txt".format(fam_name))
            fin.close()

        # get fam input
        if any("chr" in file for file in os.listdir(fam_dir)):
            print("chr information found, process to read the file")
            for nb in range(1, 23):
                chrnb = "chr{}_".format(nb)
                for f in os.listdir(fam_dir):
                    if f.endswith(".fam") and chrnb in f and fam_filename in f:
                        fam_input = "{}/{}".format(fam_dir, f)
                        fam_name = f.split(".fam")[0]
                        random_draw_samples()
        else:
            print("chr information not found, read file without chromosomes")
            for f in os.listdir(fam_dir):
                if f.endswith(".fam") and fam_filename in f:
                    fam_input = "{}/{}".format(fam_dir, f)
                    fam_name = f.split(".fam")[0]
                    random_draw_samples()

        # use subset file to create a new vcf file
        def build_new_vcf():
            print("start to subset vcf file")
            call("bcftools view -Oz -S {} {} > {}/subset_{}_{}.vcf.gz".format(sample_input,
                                                                    vcfinput,
                                                                    self.random_draw_sample_dir,
                                                                      sample_name, samplesize),shell=True)
            print("finished subset subset_{}_{}.vcf.gz file".format(sample_name, samplesize))

        vcfinput = ""
        for nb in range(1, 23):
            chrnb = "chr{}".format(nb)
            for vcf_file in os.listdir(vcf_input):
                if vcf_file.endswith('.vcf.gz') and chrnb != "chrY" and chrnb != "chrX" and chrnb != "wgs" and "{}{}".format(chrnb, symbol) in vcf_file:
                    vcfinput = "{}/{}".format(vcf_input, vcf_file)
            for j in os.listdir(self.pop_dir):
                if j.endswith(".txt") and "{}_".format(chrnb) in j and "{}".format(fam_filename) in j:
                    sample_input = "{}/{}".format(self.pop_dir, j)
                    sample_name = j.split(".txt")[0]
                    print("vcf input:{} sample input:{}".format(vcfinput, sample_input))
                    build_new_vcf()

    def random_draw_samples_from_fam(self, fam_dir, fam_filename, samplesize, tag):
        list =[]
        def process_data():
            print("start to subset the samples")
            with open(fam_input, 'r') as fin:
                lines = fin.readlines()
                for index, cols in enumerate(lines):
                    col = cols.split(" ")
                    list.append("{} {}".format(col[0], col[1]))

            with open("{}/{}.txt".format(self.pop_dir, fam_name), 'w') as fin:
                random_sample = random.sample(list, int(samplesize))
                final_list = '\n'.join(random_sample)
                fin.write(final_list)
                print("random sample list created {}.txt".format(fam_name))
            fin.close()

        if any("chr1" in file for file in os.listdir(fam_dir)):
            print("chr1 information found, process to read the file")
            chrnb = "chr1_"
            for f in os.listdir(fam_dir):
                if f.endswith(".fam") and chrnb in f and fam_filename in f:
                    fam_input = "{}/{}".format(fam_dir, f)
                    fam_name = "{}_{}".format(f.split(".fam")[0].split("chr1_")[1], tag)
                    print("input:{}".format(fam_input))
                    process_data()
        elif any("chr2" in file for file in os.listdir(fam_dir)):
            print("chr2 information found, process to read the file")
            chrnb = "chr2_"
            for f in os.listdir(fam_dir):
                if f.endswith(".fam") and chrnb in f and fam_filename in f:
                    fam_input = "{}/{}".format(fam_dir, f)
                    fam_name = "{}_{}".format(f.split(".fam")[0].split("chr2_")[1], tag)
                    print("input:{}".format(fam_input))
                    process_data()
        else:
            print("chr information not found, read file without chromosomes")
            for f in os.listdir(fam_dir):
                if f.endswith(".fam") and fam_filename in f:
                    fam_input = "{}/{}".format(fam_dir, f)
                    fam_name = "{}_{}".format(f.split(".fam")[0], tag)
                    print("input:{}".format(fam_input))
                    process_data()