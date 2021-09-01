import os
from pathlib import Path
import pandas as pd
from subprocess import call
import glob


class GPRS( object ):
    def __init__(self,
                 ref='',
                 data_dir='',
                 result_dir='./result'):
        """
        At the beginning of gprs package, result and sub-folders will generate automatically when the user first run the package.
        """
        self.ref = ref
        self.data_dir = data_dir
        self.result_dir = Path( '{}/{}'.format( os.getcwd(), result_dir ) ).resolve()
        self.setup_dir()

    def setup_dir(self):
        self.create_result_dir()
        self.create_plink_dir()
        self.create_plink_bfiles_dir()
        self.create_plink_clump_dir()
        self.create_prs_dir()
        self.create_qc_dir()
        self.create_snplists_dir()
        self.create_qc_clump_snpslist_dir()
        self.create_stat_dir()

    def create_result_dir(self):
        if not os.path.exists( self.result_dir ):
            os.mkdir( self.result_dir )

    def result_dir(self):
        return self.result_dir

    def create_plink_dir(self):
        if not os.path.exists( self.plink_dir() ):
            os.mkdir( self.plink_dir() )

    def create_stat_dir(self):
        if not os.path.exists( self.stat_dir() ):
            os.mkdir( self.stat_dir() )

    def stat_dir(self):
        return '{}/{}'.format( self.result_dir, 'stat' )

    def plink_dir(self):
        return '{}/{}'.format( self.result_dir, 'plink' )

    def create_plink_bfiles_dir(self):
        if not os.path.exists( self.plink_bfiles_dir() ):
            os.mkdir( self.plink_bfiles_dir() )

    def plink_bfiles_dir(self):
        return '{}/{}'.format( self.plink_dir(), 'bfiles' )

    def create_plink_clump_dir(self):
        if not os.path.exists( self.plink_clump_dir() ):
            os.mkdir( self.plink_clump_dir() )

    def plink_clump_dir(self):
        return '{}/{}'.format( self.plink_dir(), 'clump' )

    def create_prs_dir(self):
        if not os.path.exists( self.prs_dir() ):
            os.mkdir( self.prs_dir() )

    def prs_dir(self):
        return '{}/{}'.format( self.plink_dir(), 'prs' )

    def create_qc_dir(self):
        if not os.path.exists( self.qc_dir() ):
            os.mkdir( self.qc_dir() )

    def qc_dir(self):
        return '{}/{}'.format( self.result_dir, 'qc' )

    def create_snplists_dir(self):
        if not os.path.exists( self.snplists_dir() ):
            os.mkdir( self.snplists_dir() )

    def snplists_dir(self):
        return '{}/{}'.format( self.result_dir, 'snplists' )

    def create_qc_clump_snpslist_dir(self):
        if not os.path.exists( self.qc_clump_snpslist_dir() ):
            os.mkdir( self.qc_clump_snpslist_dir() )

    def qc_clump_snpslist_dir(self):
        return '{}/{}'.format( self.plink_dir(), 'qc_and_clump_snpslist' )


    def transfer_atcg(self, qc_file_name):
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            df = pd.read_csv( "{}/{}_{}.QC.csv".format( self.qc_dir(), chrnb, qc_file_name ), sep=' ' )
            df.loc[:, 'Allele'] = df['Allele'].apply({'a': 'A', 't': 'T', 'c': 'C', 'g': 'G'}.get)
            df.to_csv( "{}/{}_{}.QC.csv".format( self.qc_dir(), chrnb, qc_file_name ), sep=' ', index=False, header=True )
        print("transfer completed!")

    def generate_plink_bfiles(self, snplist_name, output_name, symbol='.' ):
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            for snps in os.listdir( self.snplists_dir() ):
                if "chr"  and "{}_{}.csv".format( chrnb, snplist_name ) in snps:
                        for i in os.listdir( self.ref ):
                            if i.endswith('.vcf.gz') and chrnb != "chrX" and chrnb != "chrY" and chrnb != "chrMT" and "{}{}".format(chrnb, symbol) in i:
                                    # exclude chr X, Y and MT, and setting a filter to make sure all input are consistent
                                    # ex: chr1 snps-list and chr1 vcf read at the same time
                                os.system("plink --vcf {}/{} --extract {}/{} --make-bed --out {}/{}_{}".format( self.ref, i,
                                                                                                              self.snplists_dir(),
                                                                                                              snps,
                                                                                                              self.plink_bfiles_dir(),
                                                                                                              chrnb,
                                                                                                              output_name ) )
                            print( "with chr: {}_{} is finished!".format( chrnb, output_name ) )
                elif "chr" not in snps and "{}.csv".format( snplist_name ) in snps:
                        for i in os.listdir( self.ref ):
                            if i.endswith('.vcf.gz' ) and chrnb != "chrX" and chrnb != "chrY" and chrnb != "chrMT" and "{}{}".format(chrnb, symbol ) in i:
                               os.system("plink --vcf {}/{} --extract {}/{} --make-bed --out {}/{}_{}".format( self.ref, i,
                                                                                                      self.snplists_dir(),
                                                                                                      snps,
                                                                                                      self.plink_bfiles_dir(),
                                                                                                      chrnb,
                                                                                                      output_name ) )
                            print( "no chr: {}_{} is finished!".format( chrnb, output_name ) )
        print( "all jobs completed!" )


    def clump(self, qc_file_name, plink_bfile_name, output_name, clump_kb, clump_p1, clump_p2, clump_r2='0.1', clump_field='Pvalue',clump_snp_field='SNPID'):
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            # TO perform plink clump function the following files are required:
            # .bam .bim .fam .QC(summary statistic)
            for qc_files in os.listdir( self.qc_dir() ):
                if "{}_{}.QC.csv".format( chrnb, qc_file_name ) in qc_files:
                    for plink_files in os.listdir( self.plink_bfiles_dir() ):
                        if "{}_{}.bed".format( chrnb, plink_bfile_name ) in plink_files:
                                        plinkinput = plink_files.split( "." )[0]
                                        # use plink to do the clump
                                        # NOTE: the summary statistics file has to separate by space otherwise plink can not read it
                                        os.system( "plink --bfile {}/{} --clump {}/{} --clump-p1 {} --clump-p2 {} --clump-r2 {} --clump-kb {} --clump-field {} --clump-snp-field {} --out {}/{}_{}".format(
                                                self.plink_bfiles_dir(), plinkinput,
                                                self.qc_dir(), qc_files,
                                                clump_p1, clump_p2, clump_r2,
                                                clump_kb, clump_field, clump_snp_field,
                                                self.plink_clump_dir(), chrnb, output_name ) )
                elif "{}.QC.csv".format( qc_file_name ) in qc_files:
                    for plink_files in os.listdir( self.plink_bfiles_dir() ):
                        if "{}_{}.bed".format( chrnb, plink_bfile_name ) in plink_files:
                                        plinkinput = plink_files.split( "." )[0]
                                        # use plink to do the clump
                                        # NOTE: the summary statistics file has to separate by space otherwise plink can not read it
                                        os.system( "plink --bfile {}/{} --clump {}/{} --clump-p1 {} --clump-p2 {} --clump-r2 {} --clump-kb {} --clump-field {} --clump-snp-field {} --out {}/{}_{}".format(
                                                self.plink_bfiles_dir(), plinkinput,
                                                self.qc_dir(), qc_files,
                                                clump_p1, clump_p2, clump_r2,
                                                clump_kb, clump_field, clump_snp_field,
                                                self.plink_clump_dir(), chrnb, output_name ) )
        print( "Finished Plink clump!" )

    def select_clump_snps(self, qc_file_name, clump_file_name, output_name):
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            #extract clumped snplist
            try:
                df = pd.read_csv( "{}/{}_{}.clumped".format( self.plink_clump_dir(), chrnb, output_name ), delim_whitespace=True )
                new_df = df.loc[(df['CHR'] == nb)]
                new_df_2 = new_df[['CHR', 'SNP']]
                new_df_2.to_csv("{}/{}_{}_clumped_snplist.csv".format( self.plink_clump_dir(), chrnb, output_name ), sep=' ', index=False, header=True )
            except IOError:
                print("{}/{}_{}.clumped not found. skip".format( self.plink_clump_dir(), chrnb, output_name ))
            #generate snplist by comparing clumped/qced snplists
            try:
                clump_snp = pd.read_csv("{}/{}_{}_clumped_snplist.csv".format( self.plink_clump_dir(), chrnb, clump_file_name ), sep=' ' )
                clump_snp.rename( columns={'SNP': 'SNPID'}, inplace=True )
                qc_files1 = "{}/{}_{}.QC.csv".format( self.qc_dir(), chrnb, qc_file_name )
                qc_files2 = "{}/{}.QC.csv".format( self.qc_dir(), qc_file_name )
                if os.path.exists(qc_files1):
                    print("start extracting {}".format(qc_files1))
                    qc_snp = pd.read_csv( "{}/{}_{}.QC.csv".format( self.qc_dir(), chrnb, qc_file_name ), sep=' ' )
                    newsnplist = qc_snp[qc_snp["SNPID"].isin( clump_snp["SNPID"] )]
                    newsnplist.to_csv( "{}/{}_{}.qc_clump_snpslist.csv".format( self.qc_clump_snpslist_dir(), chrnb, output_name ), sep=' ', index=False, header=True )
                    print( "{} snpslist has been created".format( chrnb ) )
                elif os.path.exists(qc_files2):
                    print( "start extracting {}".format( qc_files2 ) )
                    qc_snp = pd.read_csv( "{}/{}.QC.csv".format( self.qc_dir(), qc_file_name ), sep=' ' )
                    newsnplist = qc_snp[qc_snp["SNPID"].isin( clump_snp["SNPID"] )]
                    newsnplist.to_csv("{}/{}_{}.qc_clump_snpslist.csv".format( self.qc_clump_snpslist_dir(), chrnb, output_name ), sep=' ', index=False, header=True )
                    print( "{}x snpslist has been created".format( chrnb ) )
            except IOError:
                print("{}/{}_{}_clumped_snplist.csv not found. skip".format( self.plink_clump_dir(), chrnb, clump_file_name ))
        print( "All jobs are completed" )

    def build_prs(self, vcf_input, output_name, qc_file_name, memory, symbol='.' ,columns='1 2 3', plink_modifier='no-mean-imputation'):
        for nb in range( 1, 23 ):
            chrnb = "chr{}".format( nb )
            for vcf_file in os.listdir( vcf_input ):
                if vcf_file.endswith('.vcf.gz') and chrnb != "chrY" and chrnb != "chrX" and chrnb != "wgs" and "{}{}".format(chrnb, symbol)in vcf_file:
                    qc_file = "{}/{}_{}.qc_clump_snpslist.csv".format( self.qc_clump_snpslist_dir(), chrnb, qc_file_name )
                    try:
                        os.path.exists(qc_file)
                        os.system("plink2 --vcf {}/{} dosage=DS --score {} {} '{}' --memory {} --out {}/{}_{}".format( vcf_input, vcf_file,
                                                                                                 qc_file, columns,
                                                                                                 plink_modifier, memory,
                                                                                                 self.prs_dir(), chrnb,
                                                                                                 output_name ) )
                        print( "{} GPRS model built!".format( chrnb ) )
                    except IOError:
                        print( "{}/{}_{}.qc_clump_snpslist.csv not found. skip".format( self.qc_clump_snpslist_dir(), chrnb, qc_file_name ) )
        print( "ALL work are complete!" )

    def combine_prs(self,pop):
        ##created by - Soyoung Jeon
        ###calculate PRS from .sscore files (plink2 output)
        # needs PREFIX_chrX.sscore files as inputs
        # outputs PREFIX.sscore with same prefix
        # usage: python [prefix before '_chrX.sscore' separated by comma]

        POP = pop.split( ',' )

        chr = [i for i in range( 1, 23 )]
        for p in POP:
            nm = {}  # NMISS_ALLELE_CT
            sumsc = {}  # score as SUM of dosage*score
            avg = {}
            # make dictionary for each sample
            f = open( "{}/chr10_{}.sscore".format(self.prs_dir(), p) )
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
                    f = open( "{}/chr{}_{}.sscore".format(self.prs_dir(), i, p)  )
                    f.readline()
                    for lines in f:
                        lines = lines.split()
                        id = lines[0]
                        nm[id].append( float( lines[1] ) )
                        sumsc[id].append( float( lines[1] ) * float( lines[3] ) )
                        avg[id].append( float( lines[3] ) )
                    f.close()
                except IOError:
                    print("'File Not Found:'{}/chr{}_{}.sscore'".format(self.prs_dir(),i,p))
            # sum up scores per chromosome and output
            o = open( p + '.sscore', 'w' )
            o.write( 'id\tALLELE_CT\tSCORE_AVG\tSCORE_SUM\n' )
            for i in nm.keys():
                o.write( i + '\t' + str( sum( nm[i] ) ) + '\t' + str( sum( avg[i] ) / 22 ) + '\t' + str(
                    sum( sumsc[i] ) ) + '\n' )
            o.close()

    def prs_statistics(self, score_file, pheno_file, output_name, data_set_name, filter_pvalue, prs_stats_R, r_command):
        ##USAGE Rscript --vanilla prs_stats.R [score file] [pheno file] [target pop for OR] [ref pop for OR] [graph pdf name]
        call("{0} --vanilla {1} {2} {3} {4} {5} {7}/{6}".format(r_command, prs_stats_R, score_file, pheno_file,
                                        data_set_name, filter_pvalue, output_name, self.stat_dir()), shell=True)
        stat_data = pd.read_csv("{}/{}_{}_stat.txt".format(self.stat_dir(),output_name,filter_pvalue))
        stat_data.to_csv(
            "{}/{}_{}_stat.txt".format(self.stat_dir(),output_name,filter_pvalue),
            index=False,
            header=True,
            sep='\t',
            float_format='%.4f'
        )

    def combine_prs_stat(self,data_set_name):
        interesting_files = glob.glob( "{}/{}*.txt".format(self.stat_dir(), data_set_name) )
        combined_csv = pd.concat( [pd.read_csv( f ) for f in sorted( interesting_files )] )
        combined_csv.to_csv(
            '{}/{}_combined_stat.txt'.format( self.stat_dir(), data_set_name ),
            index=False,
            sep='\t'
        )