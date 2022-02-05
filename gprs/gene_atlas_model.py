import os
from gprs.gprs import GPRS
import pandas as pd


class GeneAtlasModel( GPRS ):
    def filter_data(self, snp_id_header, allele_header, beta_header, se_header, pvalue_header, output_name, pvalue=1):
        # extract SNPs ID for filtering
        for nb in range(1, 23):
            chrnb = "chr{}".format(nb)
            for i in os.listdir( self.data_dir ):
                if "chr" in i and i.endswith( ".csv" ) or i.endswith( ".txt" ) and ".gz" not in i and chrnb in i:
                    print( "chromosome information are find")
                    print( "start with chr" + ":" + i )
                    # read colnames by sep, and filter by pvalue
                    df = pd.read_csv( "{}/{}".format( self.data_dir, i ), delim_whitespace=True )

                    # change column name
                    df = df.rename( columns={snp_id_header: 'SNPID',
                                        allele_header: 'Allele',
                                        beta_header: 'Beta',
                                        se_header: 'SE',
                                        pvalue_header: 'Pvalue'}, inplace=True )

                    # IMPORTANT!!! Change column name before run the script
                    # extract the SNP information, notice: column name have to check
                    data = df.loc[(df['Pvalue'] <= pvalue)]

                    # reorder column, and extract three columns only
                    filtered_df = data[['SNPID', 'Allele', 'Beta', 'SE', 'Pvalue']]

                    # drop duplicate SNPs
                    filtered_df.drop_duplicates( subset=['SNPID'] )

                    filtered_df.to_csv( "{}/{}_{}.QC.csv".format( self.qc_dir(), chrnb, output_name ), sep=' ',
                                        index=False,
                                        header=True )

                    # output new csv files with chr number and SNPs ID with header
                    snps = data['SNPID']
                    snps.to_csv( "{}/{}_{}.csv".format( self.snplists_dir(), chrnb, output_name ), sep=' ',
                                 index=False, header=True )
                    print("end the job with chr" + ":" + i)

                elif "chr" not in i and i.endswith( ".csv" ) or i.endswith( ".txt" ) and ".gz" not in i:
                    print("chromosome information are not find")
                    print("start with" + ":" + i)

                    df = pd.read_csv( "{}/{}".format( self.data_dir, i ), delim_whitespace=True )
                    # change column name
                    df = df.rename( columns={snp_id_header: 'SNPID',
                                        allele_header: 'Allele',
                                        beta_header: 'Beta',
                                        se_header: 'SE',
                                        pvalue_header: 'Pvalue'}, inplace=True)

                    # IMPORTANT!!! Change column name before run the script
                    # extract the SNP information, notice: column name have to check
                    data = df.loc[(df['Pvalue'] <= pvalue)]

                    # reorder column, and extract three columns only
                    filtered_df = data[['SNPID', 'Allele', 'Beta', 'SE', 'Pvalue']]

                    # drop duplicate SNPs
                    filtered_df.drop_duplicates( subset=['SNPID'] )

                    filtered_df.to_csv( "{}/{}.QC.csv".format( self.qc_dir(), output_name ), sep=' ',
                                                    index=False,
                                                    header=True )

                    # output new csv files with chr number and SNPs ID with header
                    snps = data['SNPID']
                    snps.to_csv( "{}/{}.csv".format( self.snplists_dir(), output_name ), sep=' ',
                                             index=False, header=True )
                    print("end the job with" + ":" + i)
        print( "all Snps lists are ready!" )

        # generate the summary file
        snp_summary = []
        visited = set()
        for nb in range(1, 23):
            chrnb = 'chr{}'.format(nb)
            for i in os.listdir(self.data_dir()):
                if "chr" in i and i.endswith(".csv") or i.endswith(".txt") and i not in visited and ".gz" not in i and chrnb in i:
                    # read colnames by sep, and filter by pvalue
                    snp_file_w_chr = "{}/{}".format(self.data_dir(), i)
                    filter = open("{}".format(snp_file_w_chr), "r")
                    Counter_filter = 0

                    # Reading from file
                    Content = filter.read()
                    CoList = Content.split("\n")

                    for k in CoList:
                        if k:
                            Counter_filter += 1
                    snp_summary.append("{}:SNP BEFORE FILTERING:{}".format(i, Counter_filter))

                elif "chr" not in i and i.endswith(".csv") or i.endswith(".txt") and i not in visited and ".gz" not in i and chrnb in i:
                    snp_file_no_chr = "{}/{}".format(self.data_dir(), i)
                    Filter = open("{}".format(snp_file_no_chr), "r")
                    counter_filter = 0

                    # Reading from file
                    content = Filter.read()
                    coList = content.split("\n")

                    for j in coList:
                        if j:
                            counter_filter += 1
                    snp_summary.append("{}:Total SNP BEFORE FILTERING:{}".format(output_name, counter_filter))
                visited.add(i)

        visited.clear()
        for i in os.listdir(self.qc_dir()):
            if "chr" in i:
                for nb in range(1, 23):
                    snp_qc_file_w_chr = "{}/chr{}_{}.QC.csv".format(self.qc_dir(), nb, output_name)
                    filter = open("{}".format(snp_qc_file_w_chr), "r")
                    Counter_filter = 0

                    # Reading from file
                    Content = filter.read()
                    CoList = Content.split("\n")

                    for k in CoList:
                        if k:
                            Counter_filter += 1
                    snp_summary.append("chr{}_{}.QC:SNP AFTER FILTERING:{}".format(nb, output_name, Counter_filter))

            elif "chr" not in i:
                Snp_qc_file_no_chr = "{}/{}.QC.csv".format(self.qc_dir(), output_name)
                Filter = open("{}".format(Snp_qc_file_no_chr), "r")
                counter_filter = 0

                # Reading from file
                content = Filter.read()
                coList = content.split("\n")

                for j in coList:
                    if j:
                        counter_filter += 1
                snp_summary.append("{}:Total SNP AFTER FILTERING:{}".format(output_name, counter_filter))
        visited.add("{}/{}.QC.csv".format(self.qc_dir(), output_name))

        with open("{}/{}.filteredSNP.withPvalue{}.summary.txt".format(self.qc_dir(),output_name,pvalue), "w" ) as file:
            for k in snp_summary:
                file.write('%s\n' % k)