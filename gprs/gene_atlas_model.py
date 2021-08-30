import os
from gprs.gprs import GPRS


class GeneAtlasModel( GPRS ):
    def filter_data(self, snp_id_header, allele_header, beta_header, se_header, pvalue_header, output_name='geneatlas',
                    pvalue=0.05):
        #unzip the genotyped files
        gz_files = []
        for i in os.listdir( self.data_dir ):
            if i.endswith( ".gz" ):
                gz_files.append( i )
            else:
                print( "all files are unzipped!" )
        if gz_files:
            for gz_file in gz_files:
                gz_file = "{}/{}".format( self.data_dir, gz_file )
                unzipped_file = gz_file.replace( ".gz", "" )
                with gzip.open( gz_file, 'rb' ) as f_in, open( unzipped_file, 'wb' ) as f_out:
                    shutil.copyfileobj( f_in, f_out )
        else:
            raise ('No data files end with .gz')
        print( "starting to extract and filter SNPs" )

        # extract SNPs ID for filtering
        for i in os.listdir( self.data_dir ):
            if "chr" in i:
                print( "chromosome number were found in file name, continue to extract SNPs by chromosomes" )
                for nb in range( 1, 23 ):
                    chrnb = "chr{}".format( nb )
                    if i.endswith( ".csv" ) or i.endswith( ".txt" ) and ".gz" not in i and chrnb in i:
                        print( "with chr" + ":" + i )
                        # read colnames by sep, and filter by pvalue
                        df = pd.read_csv( "{}/{}".format( self.data_dir, i ), delim_whitespace=True )

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

            elif "chr" not in i:
                print( "chromosome number were NOT found in file name, consider a mix file, continue..." )
                if i.endswith( ".csv" ) or i.endswith( ".txt" ) and ".gz" not in i:
                    # print("without chr"+":"+i)
                    df = pd.read_csv( "{}/{}".format( self.data_dir, i ), delim_whitespace=True )
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

                    filtered_df.to_csv( "{}/{}.QC.csv".format( self.qc_dir(), output_name ), sep=' ',
                                                    index=False,
                                                    header=True )

                    # output new csv files with chr number and SNPs ID with header
                    snps = data['SNPID']
                    snps.to_csv( "{}/{}.csv".format( self.snplists_dir(), output_name ), sep=' ',
                                             index=False, header=True )
        print( "all Snps lists are ready!" )

        snp_summary=[]
        visited = set()
        for i in os.listdir( self.data_dir ):
            if  "chr" in i:
                for nb in range( 1, 23 ):
                    chrnb = 'chr{}'.format( nb )
                    if i not in visited and i.endswith( ".csv" ) or i.endswith( ".txt" ) and ".gz" not in i and chrnb in i:
                        # read colnames by sep, and filter by pvalue
                        snp_file = "{}/{}".format( self.data_dir, i )
                        filter = open( "{}".format(snp_file), "r" )
                        Counter_filter = 0

                        # Reading from file
                        Content = filter.read()
                        CoList = Content.split( "\n" )

                        for k in CoList:
                            if k:
                                Counter_filter += 1
                        snp_summary.append("{}:SNP BEFORE FILTERING:{}".format(i, Counter_filter))

            elif "chr" not in i:
                if i not in visited and i.endswith('.csv') or i.endswith('.txt') and '.gz' not in i:
                    Snp_file = "{}/{}".format( self.data_dir, i )
                    Filter = open( "{}".format( Snp_file ), "r" )
                    counter_filter = 0

                    # Reading from file
                    content = Filter.read()
                    coList = content.split( "\n" )

                    for j in coList:
                        if j:
                            counter_filter += 1
                    snp_summary.append( "{}:Total SNP BEFORE FILTERING:{}".format( output_name, counter_filter ) )
            visited.add( i )

        visited.clear()
        for i in os.listdir( self.data_dir ):
            if i not in visited:
                if  "chr" in i:
                    for nb in range( 1, 23 ):
                        snp_qc_file = "{}/chr{}_{}.QC.csv".format( self.qc_dir(), nb, output_name )
                        filter = open( "{}".format(snp_qc_file), "r" )
                        Counter_filter = 0

                        # Reading from file
                        Content = filter.read()
                        CoList = Content.split( "\n" )

                        for k in CoList:
                            if k:
                                Counter_filter += 1
                        snp_summary.append("chr{}_{}.QC:SNP AFTER FILTERING:{}".format(nb, output_name, Counter_filter))
                elif "chr" not in i:
                    Snp_qc_file = "{}/{}.QC.csv".format( self.qc_dir(), output_name )
                    Filter = open( "{}".format( Snp_qc_file ), "r" )
                    counter_filter = 0

                    # Reading from file
                    content = Filter.read()
                    coList = content.split( "\n" )

                    for j in coList:
                        if j:
                            counter_filter += 1
                    snp_summary.append( "{}:Total SNP AFTER FILTERING:{}".format( output_name, counter_filter ) )
            visited.add( i )

        with open( "{}/{}.filteredSNP.withPvalue{}.summary.txt".format(self.qc_dir(),output_name,pvalue), "w" ) as file:
            for k in snp_summary:
                file.write( '%s\n' % k )