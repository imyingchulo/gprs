import os
import pandas as pd
import gzip
import shutil
from gprs.gprs import GPRS


class Gwas( GPRS ):
    def filter_data(self, snp_id_header, allele_header, beta_header, se_header, pvalue_header,
                    output_name, pvalue=0.05):

        # # unzip the genotyped files
        for i in os.listdir( self.data_dir ):
            if "genotyped" in i and i.endswith( ".gz" ):
                gz_file = "{}/{}".format( self.data_dir, i )
                unzipped_file = gz_file.replace( ".gz", "" )
                with gzip.open( gz_file, 'rb' ) as f_in, open( unzipped_file, 'wb' ) as f_out:
                    shutil.copyfileobj( f_in, f_out )
            else:
                print( "all files are unzipped!" )
        print( "starting to extract and filter SNPs" )

        # extract SNPs ID for filtering
        for i in os.listdir( self.data_dir ):
            if i.endswith( ".csv" ):
                # read colnames by sep
                df = pd.read_csv( "{}/{}".format( self.data_dir, i ), delim_whitespace=True )

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
                    snps.to_csv( "{}/chr{}_{}.csv".format( self.snplists_dir(), nb, output_name ), sep=' ', index=False,
                                 header=True )
        print( "all Snps lists are ready!" )
