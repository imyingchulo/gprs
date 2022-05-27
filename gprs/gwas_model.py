import os
import pandas as pd
from gprs.gprs_main import GPRS
import collections

class GwasModel( GPRS ):
    def filter_data(self, snp_id_header, allele_header, beta_header, se_header, pvalue_header, file_name, output_name, pvalue=1):
        # read colnames
        df = pd.read_csv( "{}/{}".format( self.data_dir, file_name ), delim_whitespace=True )

        # change column name
        df.rename( columns={snp_id_header: 'SNPID',
                            allele_header: 'Allele',
                            beta_header: 'Beta',
                            se_header: 'SE',
                            pvalue_header: 'Pvalue'}, inplace=True )

        # setting the chr number for output files name and extract SNPs
        # SNPs are select by P-value <= 0.05
        # range(1, 23) will give you the numbers from 1 to 22
        # range(22) will give you 0 to 22
        # Thus, I put range 1, 23 in this script
        for nb in range( 1, 23 ):
            # I filtered by chromosome(column 1), P-value (column 9) and RSID (column 2)
            # PLEASE make sure the column name before run this script
            filtered_df = df.loc[(df['Chr'] == nb) & (df['Pvalue'] <= pvalue)]

            # filter df (remove NaN)
            filtered_df = filtered_df.loc[filtered_df['SNPID'].notnull()]

            # reorder column, and extract three columns only
            filtered_df_2 = filtered_df[['SNPID', 'Allele', 'Beta', 'SE', 'Pvalue']]

            # replace a t c g to A T C G
            filtered_df_2.loc[:, 'Allele'] = filtered_df_2['Allele'].apply(
                {'a': 'A', 't': 'T', 'c': 'C', 'g': 'G'}.get )

            # drop duplicate SNPs
            filtered_df_2.drop_duplicates( subset=['SNPID'] )

            # output filter data
            filtered_df_2.to_csv( "{}/chr{}_{}.QC.csv".format( self.qc_dir(), nb, output_name ), sep=' ',
                                  index=False,
                                  header=True )

            # snps is a variable represents the chromosome number
            snps = filtered_df_2['SNPID']

            # output new csv files with chr number and SNPs ID with header
            snps.to_csv( "{}/chr{}_{}.csv".format( self.snplists_dir(), nb, output_name ), sep=' ', index=False, header=True )
        print("SNPs list are ready")

        df = pd.read_csv( "{}/{}".format( self.data_dir,file_name), delim_whitespace=True )
        # method 1
        # snp_summary.append(collections.Counter(df['Chr']))

        # method 2
        table = collections.defaultdict(int)
        for item in df['Chr']:
            table[item] += 1
        summary = ['{}: {}'.format(key, val) for key, val in table.items()]
        # print('\n'.join(summary))

        with open( "{}/{}.filteredSNP.withPvalue{}.summary.txt".format(self.qc_dir(),output_name,pvalue), "w" ) as file:
            file.write('\n'.join(summary))
            for nb in range( 1, 23 ):
                snp_qc_file = "{}/chr{}_{}.QC.csv".format( self.qc_dir(), nb, output_name )
                filter = open( "{}".format( snp_qc_file ), "r" )
                Counter_filter = 0

                # Reading from file
                Content = filter.read()
                CoList = Content.split( "\n" )

                for i in CoList:
                    if i:
                        Counter_filter += 1
                file.write( "\nchr{}:SNP AFTER FILTERING:{}".format( nb, Counter_filter ) )