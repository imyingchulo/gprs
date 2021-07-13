import gzip
import os
import shutil
import pandas as pd
from gprs.gprs import GPRS


class GeneAtlas( GPRS ):
    def filter_data(self, snp_id_header, allele_header, beta_header, se_header, pvalue_header, output_name='geneatlas',
                    pvalue=0.05):

        # for i in os.listdir( self.data_dir ):
        #     if "genotyped" in i and i.endswith( ".gz" ):
        #         gz_file = "{}/{}".format( self.data_dir, i )
        #         unzipped_file = gz_file.replace( ".gz", "" )
        #         with gzip.open( gz_file, 'rb' ) as f_in, open( unzipped_file, 'wb' ) as f_out:
        #             shutil.copyfileobj( f_in, f_out )
        #         print( "{} are unzipped.".format( gz_file ) )
        # print( "starting to extract and filter SNPs" )

        # extract SNPs ID for filtering
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format(nb)
            for i in os.listdir( self.data_dir ):
                if i.endswith( ".csv" ) and "{}.csv".format(chrnb) in i:
                    # read colnames by sep, and filter by pvalue
                    df = pd.read_csv( "{}/{}".format( self.data_dir, i ), sep=' ' )

                    # change column name
                    df.rename( columns={snp_id_header: 'SNPID',
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
        print( "all Snps lists are ready!" )
