import os
from pathlib import Path
import pandas as pd
from subprocess import call
import glob

# Built a class name GPRS
class GPRS( object ):
    def __init__(self,
                 ref='',
                 data_dir='',
                 result_dir='./result'):
                 # plink='$HOME/bin/plink'):
        """
        At the beginning of gprs package, result and sub-folders will generate automatically when the user first run the package.
        """
        self.ref = ref
        self.data_dir = data_dir
        self.result_dir = Path( '{}/{}'.format( os.getcwd(), result_dir ) ).resolve()
        # self.plink = plink
        self.setup_dir()

    def setup_dir(self):   # The setup_dir function is automatically create 10 folders
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

    def create_result_dir(self):  # A function to create result folder
        if not os.path.exists( self.result_dir ):
            os.mkdir( self.result_dir )

    def result_dir(self):  # Define result dir path
        return self.result_dir

    def create_plink_dir(self):  # A function to create plink folder
        if not os.path.exists( self.plink_dir() ):
            os.mkdir( self.plink_dir() )

    def create_pop_dir(self):  # A function to create pop folder
        if not os.path.exists( self.pop_dir() ):
            os.mkdir( self.pop_dir() )

    def create_stat_dir(self):  # A function to create stat folder
        if not os.path.exists( self.stat_dir() ):
            os.mkdir( self.stat_dir() )

    def stat_dir(self):  # Define stat dir path
        return '{}/{}'.format( self.result_dir, 'stat' )

    def pop_dir(self):  # Define pop dir path
        return '{}/{}'.format( self.result_dir, 'pop' )

    def plink_dir(self):  # Define plink dir path
        return '{}/{}'.format( self.result_dir, 'plink' )

    def create_plink_bfiles_dir(self):  # A function to create plink bfile folder
        if not os.path.exists( self.plink_bfiles_dir() ):
            os.mkdir( self.plink_bfiles_dir() )

    def plink_bfiles_dir(self):  # Define plink/bfile dir path
        return '{}/{}'.format( self.plink_dir(), 'bfiles' )

    def create_plink_clump_dir(self):  # A function to create plink clump folder
        if not os.path.exists( self.plink_clump_dir() ):
            os.mkdir( self.plink_clump_dir() )

    def plink_clump_dir(self):  # Define plink/clump dir path
        return '{}/{}'.format( self.plink_dir(), 'clump' )

    def create_prs_dir(self):  # A function to create prs folder
        if not os.path.exists( self.prs_dir() ):
            os.mkdir( self.prs_dir() )

    def prs_dir(self):  # Define prs dir path
        return '{}/{}'.format( self.plink_dir(), 'prs' )

    def create_qc_dir(self):   # A function to create qc folder
        if not os.path.exists( self.qc_dir() ):
            os.mkdir( self.qc_dir() )

    def qc_dir(self):  # Define qc dir path
        return '{}/{}'.format( self.result_dir, 'qc' )

    def create_snplists_dir(self):  # A function to create snplists folder
        if not os.path.exists( self.snplists_dir() ):
            os.mkdir( self.snplists_dir() )

    def snplists_dir(self):  # Define snplists dir path
        return '{}/{}'.format( self.result_dir, 'snplists' )

    def create_qc_clump_snpslist_dir(self):  # A function to create qc'd clump snplists folder
        if not os.path.exists( self.qc_clump_snpslist_dir() ):
            os.mkdir( self.qc_clump_snpslist_dir() )

    def qc_clump_snpslist_dir(self):  # Define qc'd clump snplists dir path
        return '{}/{}'.format( self.plink_dir(), 'qc_and_clump_snpslist' )

    # Transfer a t c g into capital A T C G
    def transfer_atcg(self, qc_file_name):
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            df = pd.read_csv( "{}/{}_{}.QC.csv".format( self.qc_dir(), chrnb, qc_file_name ), sep=' ' )
            df.loc[:, 'Allele'] = df['Allele'].apply({'a': 'A', 't': 'T', 'c': 'C', 'g': 'G'}.get)
            df.to_csv( "{}/{}_{}.QC.csv".format( self.qc_dir(), chrnb, qc_file_name ), sep=' ', index=False, header=True )
        print("transfer completed!")


    # Using plink to generate bfiles fam/bim/bed.
    def generate_plink_bfiles(self, snplist_name, output_name, symbol='.' ,extra_commands=" "):
        # Generate chr with different numbers, i.e. chr1, chr2, chr3 ... chr22
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            # Get a list of snplist files from the directory
            for snps in os.listdir( self.snplists_dir() ):
                # Get a snplist file with "chr". Take one snplist file each time and run plink command
                if "chr" and "{}_{}.csv".format( chrnb, snplist_name ) in snps:
                        # Get the list of  reference files from self.ref
                        for i in os.listdir( self.ref ):
                            # The if condition here is to exclude chr X, Y, and MT, and make sure all inputs are consistent
                            # ex: The input should be chr1 snps-list and chr1 vcf reference
                            if i.endswith('.vcf.gz') and chrnb != "chrX" and chrnb != "chrY" and chrnb != "chrMT" and "{}{}".format(chrnb, symbol) in i:
                                # Run plink command
                                os.system("plink --vcf {}/{} --extract {}/{} {} --make-bed --out {}/{}_{}".format( self.ref, i,
                                                                                                              self.snplists_dir(),snps,
                                                                                                              extra_commands,
                                                                                                              self.plink_bfiles_dir(),
                                                                                                              chrnb,
                                                                                                              output_name ) )
                            print( "with chr: {}_{} is finished!".format( chrnb, output_name ) )
                # Get a snplist file without "chr". Take one snplist file each time and run plink command
                elif "chr" not in snps and "{}.csv".format( snplist_name ) in snps:
                        for i in os.listdir( self.ref ):
                            if i.endswith('.vcf.gz' ) and chrnb != "chrX" and chrnb != "chrY" and chrnb != "chrMT" and "{}{}".format(chrnb, symbol ) in i:
                               os.system("plink --vcf {}/{} --extract {}/{} {} --make-bed --out {}/{}_{}".format( self.ref, i,
                                                                                                      self.snplists_dir(),snps,
                                                                                                      extra_commands,
                                                                                                      self.plink_bfiles_dir(),
                                                                                                      chrnb,
                                                                                                      output_name ) )
                            print( "no chr: {}_{} is finished!".format( chrnb, output_name ) )
        print( "all jobs completed!" )

    #  Sub-setting the population from phenotype data
    def subset_pop(self,input_data,column_name,pop_info,output_name):
        file = pd.read_csv("{}".format(input_data), delim_whitespace=True)
        subset_file = file.loc[(file["{}".format(column_name)] == "{}".format(pop_info))]
        subset_file.to_csv("{}".format(output_name), sep='\t', index=False, header=True)

        # file = pd.read_csv( "{}".format( input_data ), delim_whitespace=True )
        # subset_file = file.loc[(file["{}".format( column_name )] == "{}".format( pop ))]
        # popinfo = subset_file['sample'] + "\t" + subset_file['sample']
        # popinfo.to_csv( "{}/{}.txt".format( self.pop_dir(),output_name ),index=False, header=False )

    # Subset Plink bfiles with selected individuals
    def generate_plink_bfiles_w_individual_info(self, popfile_name, bfile_name, output_name):
        # Generate chr number (chr1-chr22)
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            # Find the popfile from pop directory
            for j in os.listdir(self.pop_dir()):
                if "{}.txt".format(popfile_name) in j:
                    # Get the original Plink bfiles list from plink/bfile directory
                    for i in os.listdir( self.plink_bfiles_dir() ):
                        if "{}_{}".format( chrnb, bfile_name ) in i:
                            bfile = i.split(".")[0]
                            # Run plink command
                            os.system("plink --bfile {}/{} --keep {}/{} --make-bed --out {}/{}_{}".format( self.plink_bfiles_dir(),bfile,
                                                                                                        self.pop_dir(),j,
                                                                                                        self.plink_bfiles_dir(),
                                                                                                        chrnb,output_name ) )
                        print( "with chr: {}_{} is finished!".format( chrnb, bfile_name ) )
        print( "all jobs completed!" )

    # This function is doing clumping by using Plink
    def clump(self, qc_file_name, plink_bfile_name, output_name, clump_kb, clump_p1, clump_p2, clump_r2='0.1', clump_field='Pvalue', clump_snp_field='SNPID'):
        # Create a C+T tag
        output_name_with_conditions = "{}_{}_{}_{}".format( output_name, clump_kb, clump_p1, clump_r2 )
        # Generate chr number (chr1-chr22)
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            # Get the list of .QC files (unified summary statistic from step1)
            for qc_files in os.listdir( self.qc_dir() ):
                # If file name contains "chr", take one file each time
                # If filename does not contain "chr", then jump to next "elif"
                if "{}_{}.QC.csv".format( chrnb, qc_file_name ) in qc_files:
                    # Get a list of plink bfiles
                    for plink_files in os.listdir( self.plink_bfiles_dir() ):
                        # Take one .bed file each time
                        if "{}_{}.bed".format( chrnb, plink_bfile_name ) in plink_files:
                                        plinkinput = plink_files.split( "." )[0]
                                        # Check the folder exists or not, if yes, run plink command
                                        # Check the folder exists or not, if not, jump to "else" create the output folder and then run the plink command
                                        if os.path.exists("{}/{}".format(self.plink_clump_dir(), output_name_with_conditions)):
                                            # After defined all the input files(qc and .bed files), Run plink
                                            # NOTE: the summary statistics file has to separate by space otherwise plink can not read it
                                            os.system( "plink --bfile {}/{} --clump {}/{} --clump-p1 {} --clump-p2 {} --clump-r2 {} --clump-kb {} --clump-field {} --clump-snp-field {} --out {}/{}/{}_{}".format(
                                                    self.plink_bfiles_dir(), plinkinput,
                                                    self.qc_dir(), qc_files,
                                                    clump_p1, clump_p2, clump_r2,
                                                    clump_kb, clump_field, clump_snp_field,
                                                    self.plink_clump_dir(), output_name_with_conditions, chrnb, output_name_with_conditions) )
                                        else:
                                            # Create a folder with C+T tag under plink/clump
                                            os.mkdir("{}/{}".format(self.plink_clump_dir(), output_name_with_conditions))
                                            os.system("plink --bfile {}/{} --clump {}/{} --clump-p1 {} --clump-p2 {} --clump-r2 {} --clump-kb {} --clump-field {} --clump-snp-field {} --out {}/{}/{}_{}".format(
                                                    self.plink_bfiles_dir(), plinkinput,
                                                    self.qc_dir(), qc_files,
                                                    clump_p1, clump_p2, clump_r2,
                                                    clump_kb, clump_field, clump_snp_field,
                                                    self.plink_clump_dir(), output_name_with_conditions, chrnb, output_name_with_conditions ) )
                # Take one qc file(without chr in the file name) each time
                elif "{}.QC.csv".format( qc_file_name ) in qc_files:
                    for plink_files in os.listdir( self.plink_bfiles_dir() ):
                        if "{}_{}.bed".format( chrnb, plink_bfile_name ) in plink_files:
                                        plinkinput = plink_files.split( "." )[0]
                                        # Check the folder exists or not, if yes, run plink command
                                        # Check the folder exists or not, if not, jump to "else" and create the output folder
                                        if os.path.exists("{}/{}".format( self.plink_clump_dir(), output_name_with_conditions ) ):
                                            os.system( "plink --bfile {}/{} --clump {}/{} --clump-p1 {} --clump-p2 {} --clump-r2 {} --clump-kb {} --clump-field {} --clump-snp-field {} --out {}/{}/{}_{}".format(
                                                    self.plink_bfiles_dir(), plinkinput,
                                                    self.qc_dir(), qc_files,
                                                    clump_p1, clump_p2, clump_r2,
                                                    clump_kb, clump_field, clump_snp_field,
                                                    self.plink_clump_dir(),output_name_with_conditions,  chrnb, output_name_with_conditions) )
                                        else:
                                            # Create a folder with C+T tag under plink/clump
                                            os.mkdir("{}/{}".format( self.plink_clump_dir(), output_name_with_conditions ) )
                                            os.system("plink --bfile {}/{} --clump {}/{} --clump-p1 {} --clump-p2 {} --clump-r2 {} --clump-kb {} --clump-field {} --clump-snp-field {} --out {}/{}/{}_{}".format(
                                                    self.plink_bfiles_dir(), plinkinput,
                                                    self.qc_dir(), qc_files,
                                                    clump_p1, clump_p2, clump_r2,
                                                    clump_kb, clump_field, clump_snp_field,
                                                    self.plink_clump_dir(), output_name_with_conditions, chrnb, output_name_with_conditions ) )
        print( "Finished Plink clump!" )

    # This function is to extract the snps after clumping (the snps passed C+T)
    def select_clump_snps(self, qc_file_name, clump_file_name, clumpfolder_name, output_name, clump_kb, clump_p1, clump_r2 ):
        # Create a C+T tag
        clump_conditions = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2 )
        # Generate chr number (chr1-chr22)
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format(nb)
            # I separate the function into two parts:
            # 1. generate [filename]_clumped_snplist.csv; 2. extract out the qc information for clumped snps
            # PART1: generate [filename]_clumped_snplist.csv
            # Use "try" and "except IOError" to allow the missing input files
            try:
                # Check the folder exists or not, if yes, run command
                # Check the folder exists or not, if not, jump to "else" and create the output folder
                if os.path.exists("{}/{}_{}".format(self.plink_clump_dir(),clumpfolder_name, clump_conditions)):
                    # Use pandas to read the file, and extract the snps then out put the clumped_snplist.csv
                    df = pd.read_csv( "{}/{}_{}/{}_{}_{}.clumped".format( self.plink_clump_dir(),clumpfolder_name, clump_conditions, chrnb, clump_file_name, clump_conditions ) , delim_whitespace=True )
                    new_df = df.loc[(df['CHR'] == nb)]
                    new_df_2 = new_df[['CHR', 'SNP']]
                    new_df_2.to_csv("{}/{}_{}/{}_{}_{}_clumped_snplist.csv".format( self.plink_clump_dir(),clumpfolder_name,clump_conditions, chrnb, clump_file_name, clump_conditions ), sep=' ', index=False, header=True )
                else:
                    # Create the folder with C+T tag
                    os.mkdir("{}/{}_{}".format(self.plink_clump_dir(),clumpfolder_name, clump_conditions))
                    df = pd.read_csv( "{}/{}_{}/{}_{}_{}.clumped".format( self.plink_clump_dir(),clumpfolder_name, clump_conditions, chrnb, clump_file_name, clump_conditions ) , delim_whitespace=True )
                    new_df = df.loc[(df['CHR'] == nb)]
                    new_df_2 = new_df[['CHR', 'SNP']]
                    new_df_2.to_csv("{}/{}_{}/{}_{}_{}_clumped_snplist.csv".format( self.plink_clump_dir(),clumpfolder_name,clump_conditions, chrnb, clump_file_name, clump_conditions ), sep=' ', index=False, header=True )
            except IOError:
                print("{}/{}_{}/{}_{}_{}.clumped not found. skip".format( self.plink_clump_dir(),clumpfolder_name,clump_conditions, chrnb,clump_file_name, clump_conditions ))
            print("{}_{}_{} clumped snplist created".format(chrnb,clump_file_name, clump_conditions) )
        print("All jobs finished")

        # PART2: generate snplist with qc information (Allele Beta SE Pvalue)
        # Check the folder exists or not
        if os.path.exists( "{}/{}_{}".format( self.qc_clump_snpslist_dir(), clumpfolder_name, clump_conditions ) ):
            pass
        else:
            os.mkdir( "{}/{}_{}".format( self.qc_clump_snpslist_dir(), clumpfolder_name, clump_conditions ) )

        # Use visited to make sure read the file only once
        visited = set()
        # Generate chr number (chr1-chr22)
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            # Get the full path to the clump snps file
            clump_snp_file = ("{}/{}_{}/{}_{}_{}_clumped_snplist.csv".format( self.plink_clump_dir(),
                                                                              clumpfolder_name,
                                                                              clump_conditions,
                                                                              chrnb, clump_file_name, clump_conditions ))
            # Some GWAS summary statistics does not contain "chr" and some does
            # Thus, define two type of qc file as qc_file1 and qc_file2
            qc_files1 = "{}/{}_{}.QC.csv".format( self.qc_dir(), chrnb, qc_file_name )
            qc_files2 = "{}/{}.QC.csv".format( self.qc_dir(), qc_file_name )

            # If qc_files1 not exists, jump to "elif" and read the qc_file2
            # If qc_file1 exists, then run the command below
            if os.path.exists( qc_files1 ) and  qc_files1 not in visited:
                print("start extracting {}".format(qc_files1))
                # Use "try" and "except IOError" to allow missing files
                # use "if" to check the file is exists or not. If file not exists print the error message
                try:
                    if os.path.exists(clump_snp_file):
                        clump_snp = pd.read_csv( clump_snp_file, delim_whitespace=True  )
                        clump_snp.rename( columns={'SNP': 'SNPID'}, inplace=True )
                        qc_snp = pd.read_csv( "{}/{}/{}_{}.QC.csv".format( self.qc_dir(), clump_conditions, chrnb, qc_file_name ), delim_whitespace=True )
                        newsnplist = qc_snp[qc_snp["SNPID"].isin( clump_snp["SNPID"] )]
                        newsnplist.to_csv( "{}/{}_{}/{}_{}_{}.qc_clump_snpslist.csv".format( self.qc_clump_snpslist_dir(),clumpfolder_name, clump_conditions, chrnb, output_name, clump_conditions ), sep=' ', index=False, header=True )
                        print( "{}_{} snpslist created".format( chrnb,clump_snp_file ) )
                except IOError:
                    print("{} not found. skip".format(clump_snp_file ) )
            elif os.path.exists(qc_files2) and  qc_files2 not in visited:
                print( "start extracting {}".format( qc_files2 ) )
                try:
                    if os.path.exists(clump_snp_file):
                        clump_snp = pd.read_csv( clump_snp_file, delim_whitespace=True  )
                        clump_snp.rename( columns={'SNP': 'SNPID'}, inplace=True )
                        qc_snp = pd.read_csv( "{}/{}.QC.csv".format( self.qc_dir(), qc_file_name ), delim_whitespace=True )
                        newsnplist = qc_snp[qc_snp["SNPID"].isin( clump_snp["SNPID"] )]
                        newsnplist.to_csv("{}/{}_{}/{}_{}_{}.qc_clump_snpslist.csv".format( self.qc_clump_snpslist_dir(),clumpfolder_name, clump_conditions, chrnb, output_name, clump_conditions ), sep=' ', index=False, header=True )
                        print( "{}_{} snpslist created".format( chrnb,clump_snp_file ))
                except IOError:
                    print( "{} not found. skip".format( clump_snp_file ))
        print( "All jobs are completed" )

    # Use build_prs to generate the prs model
    def build_prs(self, vcf_input, output_name, qc_clump_snplist_foldername, memory, clump_kb, clump_p1, clump_r2, symbol='.' ,columns='1 2 3', plink_modifier='no-mean-imputation'):
        # Create a C+T tag
        clump_conditions = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2 )
        # Check the folder exists or not, if not create the folder
        if os.path.exists( "{}/{}_{}".format( self.prs_dir(), output_name, clump_conditions ) ):
            pass
        else:
            os.mkdir( "{}/{}_{}".format( self.prs_dir(), output_name, clump_conditions ) )
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            # check file status
            if os.path.exists("{}/{}_{}/{}_{}_{}.sscore".format(self.prs_dir(), output_name, clump_conditions, chrnb, output_name, clump_conditions)):
                print("{}_{}_{}.sscore already finished!".format(chrnb, output_name, clump_conditions))
            else:
                for vcf_file in os.listdir( vcf_input ):
                    # Define the input files (vcf and qc files)
                    if vcf_file.endswith('.vcf.gz') and chrnb != "chrY" and chrnb != "chrX" and chrnb != "wgs" and "{}{}".format(chrnb, symbol) in vcf_file:
                        qc_file = "{}/{}_{}/{}_{}_{}.qc_clump_snpslist.csv".format( self.qc_clump_snpslist_dir(), qc_clump_snplist_foldername, clump_conditions, chrnb, qc_clump_snplist_foldername, clump_conditions )
                        try:
                            os.path.exists(qc_file)
                            os.system("plink2 --vcf {}/{} dosage=DS --score {} {} '{}' --memory {} --out {}/{}_{}/{}_{}_{}".format( vcf_input, vcf_file,
                                                                                                     qc_file, columns, plink_modifier,
                                                                                                     memory,
                                                                                                     self.prs_dir(), output_name, clump_conditions,
                                                                                                     chrnb, output_name, clump_conditions ) )
                            print( "{} GPRS model built!".format( chrnb ) )
                        except IOError:
                            print( "{}/{}/{}_{}_{}.qc_clump_snpslist.csv not found. skip".format( self.qc_clump_snpslist_dir(), clump_conditions, chrnb, qc_clump_snplist_foldername, clump_conditions ) )
                print("{}_{}_{}.sscore completed!".format(chrnb, output_name, clump_conditions))
        print( "ALL work are complete!" )

    def combine_prs(self,filename,clump_kb,clump_p1,clump_r2):
        ##created by - Soyoung Jeon
        ###calculate PRS from .sscore files (plink2 output)
        # needs PREFIX_chrX.sscore files as inputs
        # outputs PREFIX.sscore with same prefix
        # usage: python [prefix before '_chrX.sscore' separated by comma]

        POP = "{}_{}_{}_{}".format(filename,clump_kb,clump_p1,clump_r2).split( ',' )
        dir = "{}_{}_{}_{}".format(filename,clump_kb,clump_p1,clump_r2)

        if os.path.exists("{}/{}".format(self.prs_dir(),dir)):
            chr = [i for i in range( 1, 23 )]
            for p in POP:
                nm = {}  # NMISS_ALLELE_CT
                sumsc = {}  # score as SUM of dosage*score
                avg = {}
                # make dictionary for each sample
                if os.path.exists("{}/{}/chr10_{}.sscore".format(self.prs_dir(),dir, p)):
                    f = open( "{}/{}/chr10_{}.sscore".format(self.prs_dir(),dir, p) )
                    f.readline()
                    for lines in f:
                        lines = lines.split()
                        nm[lines[0]] = []
                        sumsc[lines[0]] = []
                        avg[lines[0]] = []
                else:
                    f = open( "{}/{}/chr1_{}.sscore".format( self.prs_dir(), dir, p ) )
                    f.readline()
                    for lines in f:
                        lines = lines.split()
                        nm[lines[0]] = []
                        sumsc[lines[0]] = []
                        avg[lines[0]] = []
                # sum PRS per chromosome
                for i in chr:
                    i = str( i )
                    try:
                        f = open( "{}/{}/chr{}_{}.sscore".format(self.prs_dir(),dir, i, p)  )
                        f.readline()
                        for lines in f:
                            lines = lines.split()
                            id = lines[0]
                            nm[id].append( float( lines[1] ) )
                            sumsc[id].append( float( lines[1] ) * float( lines[3] ) )
                            avg[id].append( float( lines[3] ) )
                        f.close()
                    except IOError:
                        print("'File Not Found:'{}/{}/chr{}_{}.sscore'".format(self.prs_dir(),dir, i,p))
                # sum up scores per chromosome and output
                o = open('{}/'.format(self.prs_dir()) + p + '.sscore', 'w' )
                o.write( 'id\tALLELE_CT\tSCORE_AVG\tSCORE_SUM\n' )
                for i in nm.keys():
                    o.write( i + '\t' + str( sum( nm[i] ) ) + '\t' + str( sum( avg[i] ) / 22 ) + '\t' + str(
                        sum( sumsc[i] ) ) + '\n' )
                o.close()
                # remove the .sscore file after combined sscore files
                call("rm -rf {}/{}".format(self.prs_dir(),dir),shell=True)
        else:
            print("{}/{} not exsist".format(self.prs_dir(),dir))


    # Calculate the PRS statistical results and output the statistics summary
    def prs_statistics(self, score_file, pheno_file, output_name, data_set_name, prs_stats_R, r_command, clump_kb, clump_p1, clump_r2 ):
        filter_condition = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2 )
        # Sum the SNPs number from .qc_clump_snpslist.csv file
        lines=[]
        for i in range(1,23):
            try:
                file = "{}/{}_{}/chr{}_{}_{}.qc_clump_snpslist.csv".format(self.qc_clump_snpslist_dir(),
                                                                           data_set_name,filter_condition,i,data_set_name,filter_condition)
                df = pd.read_csv("{}".format(file))
                # print("{}/{}".format(self.qc_clump_snpslist_dir(),file))
                index = df.index
                number_of_rows = len(index)
                lines.append(number_of_rows)
            except IOError:
             print( "{}/{}_{}/chr{}_{}_{}.qc_clump_snpslist.csv Not Found".format(self.qc_clump_snpslist_dir(),
                                                                                  data_set_name,filter_condition,i,data_set_name,filter_condition) )

        if os.path.exists(score_file):
            # The R script is written by Soyoung Jeon
            # USAGE Rscript --vanilla prs_stats_quantitative_phenotype.R [score file] [pheno file] [target pop for OR] [ref pop for OR] [graph pdf name]
            call("{0} --vanilla {1} {2} {3} {4} {5} {6} {8}/{7}".format(r_command, prs_stats_R, score_file, pheno_file,
                                            data_set_name, filter_condition, sum(lines),output_name, self.stat_dir()), shell=True)
            # Read the statistics result and reformat it
            # Reformat: separator = tab, change float into scientific notation
            stat_data = pd.read_csv("{}/{}_{}_stat.txt".format(self.stat_dir(),output_name, filter_condition ))
            try:
                os.path.exists("{}/{}_{}_stat.txt".format(self.stat_dir(),output_name, filter_condition ))
                stat_data.to_csv(
                    "{}/{}_{}_stat.txt".format(self.stat_dir(),output_name, filter_condition),
                    index=False,
                    header=True,
                    sep='\t',
                    float_format='%.2E' )
            except IOError:
                print( "{}/{}_{}_stat.txt Not Found".format(self.stat_dir(),output_name, filter_condition ))
        else:
            print("{} Not Found".format(score_file))

    # In combine_prs_stat function is to combine PRS statistical results as one file
    def combine_prs_stat(self,data_set_name, clump_kb, clump_p1, clump_r2):
        filter_condition = "{}_{}_{}".format(clump_kb, clump_p1, clump_r2 )
        extension = "_stat.txt"
        target_files = [i for i in glob.glob("{}/{}*{}".format(self.stat_dir(), data_set_name, extension))]
        # Combine all the csv files
        combined_csv = pd.concat( [pd.read_csv( f ) for f in target_files] )
        # Export to csv
        combined_csv.to_csv( "{}/combined_{}.txt".format(self.stat_dir(),data_set_name), index=False,header=True,sep=' ')
