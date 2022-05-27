import os
from gprs.gprs_main import GPRS
import pandas as pd


class GeneAtlasModel(GPRS):
    def filter_data(self, snp_id_header, allele_header, beta_header, se_header, pvalue_header, output_name, pvalue=1):
        # column name mapping
        columns = {
            snp_id_header: 'SNPID',
            allele_header: 'Allele',
            beta_header: 'Beta',
            se_header: 'SE',
            pvalue_header: 'Pvalue'
        }

        def unify_filter():
            print("Start process {}".format(input_file))
            df = pd.read_csv(input_file, delim_whitespace=True)
            if pvalue_header in df:
                df.rename(columns=columns, inplace=True)
                # IMPORTANT!!! Change column name before run the script
                # extract the SNP information, notice: column name have to check
                data = df.loc[(df['Pvalue'] <= pvalue)]
                # reorder column, and extract three columns only
                filtered_df = data[['SNPID', 'Allele', 'Beta', 'SE', 'Pvalue']]
                # drop duplicate SNPs
                filtered_df.drop_duplicates(subset=['SNPID'])
                filtered_df.to_csv(qc_output_file, sep=' ', index=False, header=True)
                print("QC: {}".format(qc_output_file))
                snps = data['SNPID']
                snps.to_csv(snp_output_file, sep=' ', index=False, header=True)
                print("SNP List: {}".format(snp_output_file))
            else:
                print("file: {} fall into condition but should be skipped".format(input_file))

        if any("chr" in file and (file.endswith(".csv") or file.endswith(".txt")) and ".gz" not in file for file in os.listdir(self.data_dir)):
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                print("chrnb: {}".format(chrnb))
                for i in os.listdir(self.data_dir):
                    if "chr" in i and i.endswith(".csv") or i.endswith(".txt") and ".gz" not in i and chrnb in i:
                        print("chromosome information are found")
                        input_file = "{}/{}".format(self.data_dir, i)
                        qc_output_file = "{}/{}_{}.QC.csv".format(self.qc_dir, chrnb, output_name)
                        snp_output_file = "{}/{}_{}.csv".format(self.snplists_dir, chrnb, output_name)
                        unify_filter()
        else:
            for i in os.listdir(self.data_dir):
                if i.endswith(".csv") or i.endswith(".txt") and ".gz" not in i:
                    print("chromosome information are not found")
                    input_file = "{}/{}".format(self.data_dir, i)
                    qc_output_file = "{}/{}.QC.csv".format(self.qc_dir, output_name)
                    snp_output_file = "{}/{}.csv".format(self.snplists_dir, output_name)
                    unify_filter()
        print("all Snps lists are ready!")

        # generate the summary file
        snp_summary = []
        def calculate_the_snps_number():
            visited = set()
            df = pd.read_csv(raw_input, delim_whitespace=True)
            if pvalue_header in df and "{}-{}".format(raw_input, cured_input) not in visited:
                raw_filter = open("{}".format(raw_input), "r")
                print("start to process {}".format(raw_input))
                Counter_raw_filter = 0
                Content_raw = raw_filter.read()
                CoList_raw = Content_raw.split("\n")

                for k in CoList_raw:
                    if k:
                        Counter_raw_filter += 1
                snp_summary.append("Raw file name:{} SNP BEFORE FILTERING:{}".format(raw_input, Counter_raw_filter))

                print("finished {} continuous to process {}".format(raw_filter, cured_input))
                cured_filter = open("{}".format(cured_input), "r")
                Counter_filter = 0
                Content = cured_filter.read()
                CoList = Content.split("\n")

                for j in CoList:
                    if j:
                        Counter_filter += 1
                snp_summary.append("Unified & filtered file name:{} SNP AFTER FILTERING:{}".format(cured_input, Counter_filter))
                visited.add("{}-{}".format(raw_input, cured_input))
            else:
                print("This .txt or .csv file is not a GWAS summary statistics data, Skip")


        if any("chr" in file and (file.endswith(".csv") or file.endswith(".txt")) and ".gz" not in file for file in os.listdir(self.data_dir)):
            print("chromosome information are found")
            for nb in range(1, 23):
                chrnb = "chr{}".format(nb)
                print("chrnb: {}".format(chrnb))
                for i in os.listdir(self.data_dir):
                    if "chr" in i and i.endswith(".csv") or i.endswith(".txt") and ".gz" not in i and chrnb in i:
                        raw_input = "{}/{}".format(self.data_dir, i)
                        cured_input = "{}/chr{}_{}.QC.csv".format(self.qc_dir, nb, output_name)
                        print("raw_input:{} \ncured_input:{}".format(raw_input, cured_input))
                        calculate_the_snps_number()
        else:
            for i in os.listdir(self.data_dir):
                if i.endswith(".csv") or i.endswith(".txt") and ".gz" not in i:
                    print("chromosome information are not found")
                    raw_input = "{}/{}".format(self.data_dir, i)
                    cured_input = "{}/{}.QC.csv".format(self.qc_dir, output_name)
                    print("raw_input:{} \ncured_input:{}".format(raw_input, cured_input))
                    calculate_the_snps_number()

        with open("{}/{}.filteredSNP.withPvalue{}.summary.txt".format(self.qc_dir, output_name, pvalue), "w") as file:
            for k in snp_summary:
                file.write('%s\n' % k)
        print("summary txt generated")