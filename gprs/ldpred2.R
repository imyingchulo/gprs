#!/usr/bin/Rscript

# Options
suppressMessages(library("optparse"))
option_list=list(
    make_option("--train", action="store", default=NA, type='character',
        help="Prefix to training PLINK files including full path [required]"),
    make_option("--LDref",action="store", default=NA, type='character', 
        help="Path to containing external LD reference PLINK files [optional]"),
    make_option("--LDmatrix",action="store", default='./tmp-data/LD_matrix', type='character',
        help="Path to LD matrix directory [optional]"),
    make_option("--sumstat", action="store", default=NA, type='character',
        help="Prefix to GWAS summary statistics including full path[required]"),
    make_option("--output_dir", action="store", default=NA, type='character',
        help="Path for output files [required]"),
    make_option("--h2", action="store", default=NA, type='numeric',
        help="heritability estimate [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))
if (is.na(opt$LDref)){
  opt$LDref = opt$train
}
cat('Options are:\n')
print(opt)

#############
# 0. prepare workspace
#############
#LDPred2
suppressMessages(library(bigsnpr))
suppressMessages(library(bigreadr))
suppressMessages(options(bigstatsr.check.parallel.blas = FALSE))
suppressMessages(options(default.nproc.blas = NULL))
suppressMessages(library(data.table))
suppressMessages(library(magrittr))

#Evaluation
suppressMessages(library(fmsb))
suppressMessages(library(pROC))

system(paste0('mkdir -p ',opt$output,'/'))
bigassertr::assert_dir(opt$LDmatrix)
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

info <- readRDS("./tmp-data/map_hm3_ldpred2.rds")
NCORES <- nb_cores() #as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))?

#############
# 1.  Read in GWAS summary statistics; assume pre-processed from gprs prepare-sumstat
#############

cat('Reading in Sumstat...\n')
sumstats <- data.table()

#load by chr
for (chrnb in 1:22){
  sumstat <- bigreadr::fread2(paste0(opt$sumstat,'_chr',chrnb,'.csv'), data.table=TRUE)
  cat(paste0('Processing chr',chrnb,'\n'))
  names(sumstat) <- c("snpid","chr","pos","a1","a0","beta","beta_se","p","n_eff")
  cat(paste0(nrow(sumstat),' SNPs in file... '))
  sumstat <- sumstat[sumstat$pos %in% info$pos, ]
  cat(paste0(nrow(sumstat),' SNPs in hapmap3\n')) 
  sumstats <- rbindlist( list(sumstats, sumstat))
}
cat(paste0(nrow(sumstats),' SNPs in total... \n\n'))



#############
# 2. LD
#############
# Read in LD Reference
cat('------------------------------------------------\n')
cat('Reading in LD Reference..\n')
if(file.exists(paste0(opt$LDref,'.bk'))){
    system(paste0('rm ',opt$LDref,'.bk'))
}
snp_readBed(paste0(opt$LDref,'.bed'))
obj.bigSNP <- snp_attach(paste0(opt$LDref,".rds"))
cat(paste0('LD Reference file contains ',dim(obj.bigSNP$genotypes)[1],' individuals, ',dim(obj.bigSNP$genotypes)[2],' SNPs\n\n'))
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
cat('Matching SNPs with Sumstat and HapMap3 \n')
info_snp <- snp_match(sumstats, map)
info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp))

G_imp <- snp_fastImputeSimple(obj.bigSNP$genotypes, method = "mean2", ncores=NCORES)
G <- G_imp
CHR <- map$chr
POS <- map$pos
POS2 <- snp_asGeneticPos(CHR, POS, dir = "./tmp-data", ncores = NCORES)


# Compute LD
cat('Computing LD...\n')
for(chr in 1:22){
  cat(chr, ".. ", sep = "")
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  corr <- snp_cor(
    G, 
    ind.col = ind.chr2, 
    infos.pos = POS2[ind.chr2], 
    size = 3 / 1000, ncores = NCORES)
  saveRDS(corr, file = paste0(opt$LDmatrix,'/LD_chr', chr, ".rds"), version = 2)
}
cat(paste0('LD per chromosome saved in ',opt$LDmatrix,'\n'))

# Compute LD Scores
cat('\nComputing LD Scores...\n')
map <- map[info_snp$`_NUM_ID_`,]
map$ld <- do.call('c', lapply(1:22, function(chr) {
  cat(chr, ".. ", sep = "")
  corr_chr <- readRDS(paste0(opt$LDmatrix,'/LD_chr', chr, ".rds"))
  Matrix::colSums(corr_chr^2)
}))

#add positions in different builds - haven't tested
#need to download liftover file
#liftOver <- runonce::download_file("http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver", "tmp-data")
#bigsnpr:::make_executable(liftOver)
#map$pos_hg17 <- snp_modifyBuild(map, liftOver, from = "hg19", to = "hg17")$pos
#map$pos_hg18 <- snp_modifyBuild(map, liftOver, from = "hg19", to = "hg18")$pos
#map$pos_hg38 <- snp_modifyBuild(map, liftOver, from = "hg19", to = "hg38")$pos
map_ldref <- map
info_snp <- snp_match(sumstats, map_ldref)
info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp))
saveRDS(map, paste0(opt$LDmatrix,'/map.rds'), version = 2)
cat('\nLD score saved in ',opt$LDmatrix,'\n\n\n')


#############
# 3. Read in Training data
#############
# Genotype plink file
cat('------------------------------------------------\n')
cat('Reading in Training dataset...\n')
if(opt$train == opt$LDref){
    cat('Training dataset is used as LD reference\n')
    map_test <- map
    } else { #아직안해봄,info_snp아직건드리면안되나?
    if(file.exists(paste0(opt$train,'.bk'))){
      system(paste0('rm ',opt$train,'.bk'))}
    snp_readBed(paste0(opt$train,'.bed'))
    obj.bigSNP <- snp_attach(paste0(opt$train,".rds"))
    cat(paste0('Training dataset contains ',dim(obj.bigSNP$genotypes)[1],' individuals, ',dim(obj.bigSNP$genotypes)[2],' SNPs\n'))
    G <- obj.bigSNP$genotypes
    G_imp <- snp_fastImputeSimple(G, method = "mean2")
    G <- G_imp
    map_test <- obj.bigSNP$map[-3]
    names(map_test) <- c("chr", "rsid", "pos", "a1", "a0")
    in_test <- vctrs::vec_in(info_snp[, c("chr", "pos")], map_test[, c("chr", "pos")])
    info_snp <- info_snp[in_test, ]
}

#############
# 4. Genome-wide sparse LD matrix & LDSC
#############
cat('\n\nCreating genome-wide sprase LD matrix...\n')
for (chr in 1:22) {
    cat(chr, ".. ", sep = "")
    ## indices in 'df_beta'
    ind.chr <- which(info_snp$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    ## indices in 'corr_chr'
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
    corr_chr <- readRDS(paste0(opt$LDmatrix,"/LD_chr", chr, ".rds"))[ind.chr3, ind.chr3]
    if (chr == 1) {
        corr <- as_SFBM(corr_chr, tmp)
    } else {
    corr$add_columns(corr_chr, nrow(corr))
    }
}
cat('\n\n\n')

### LDSC - estimat heritability
map_ldref <- map_ldref[info_snp$`_NUM_ID_`, ]
if (is.na(opt$h2)){
  cat('Heritability not provided by user...\n')
  ldsc <- with(info_snp, snp_ldsc( map_ldref$ld,
                    nrow(map_ldref),
                    chi2 = (beta / beta_se)^2,
                    sample_size = n_eff, 
                    blocks = NULL))
    h2_est <- ldsc[["h2"]]
    cat('Estimated SNP heritability by LDSC = ',h2_est,'\n')
} else {
  h2_est <- opt$h2
}
cat('\n\n')

#cat('file size for correlation matrix: ')
#(file.size(corr$sbk) / 1024^3)
#cat('\n\n\n\n')


#############
# 5.1 LDPred-Inf
#############
cat('-----------------------------------------\n\n')
cat('Running LDPred models\n')
cat('Infinitestimal model..\n')
beta_inf <- snp_ldpred2_inf(corr, info_snp, h2 = h2_est)
ind.test <- 1:nrow(G)
cat('Done!\n\n')

#############
# 5.2 LDPred-grid
#############
cat('Grid model..\n')
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
grid.param <-
  expand.grid(p = p_seq,
              h2 = h2_seq,
              sparse = c(FALSE, TRUE))
beta_grid <- snp_ldpred2_grid(corr, info_snp, grid.param)


beta_grid_nosp <- data.table(beta_grid[,grid.param$sparse==F])
names(beta_grid_nosp)<-gsub('-','.',paste0(grid.param$p[grid.param$sparse == F],'_',grid.param$h2[grid.param$sparse == F],'_nosparse'))
beta_grid_sp <- data.table(beta_grid[,grid.param$sparse==T])
names(beta_grid_sp)<-gsub('-','.',paste0(grid.param$p[grid.param$sparse == T],'_',grid.param$h2[grid.param$sparse == T],'_sparse'))
cat('Done!\n\n')

#############
# 5.3 LDPred-auto
#############
cat('Auto model\n')
multi_auto <- snp_ldpred2_auto(corr, info_snp, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
                               ncores = NCORES)
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
pred_auto <- big_prodMat(G, beta_auto, ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)

#filter out bad chains by comparing scale of resulting predictions(calculating sd and keeping ones without too much divergence)
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
#get mean of filtered predictions -> final prediction
final_beta_auto <- rowMeans(beta_auto[, keep])
cat('Done!\n')


### write out Beta
cat('\n\nOutputting Betas..\n')
betas <- data.table(SNPID=info_snp$rsid, CHR=info_snp$chr, POS=info_snp$pos, A1=info_snp$a1, A2=info_snp$a0, beta_inf, beta_grid_nosp, beta_grid_sp, beta_auto=final_beta_auto)
names(betas)[-1:-5]<-paste0('LDPred2_',names(betas)[-1:-5])

rem<-NULL
for(i in 6:length(names(betas))){
    #infinite/NA beta 있는 parameter는 스킵, index저장
    if(is.infinite(sum(betas[[names(betas)[i]]])) | is.na(sum(betas[[names(betas)[i]]]))){
        cat ('Skipping',names(betas)[i],'due to presence of infinite/null values.\n')
        rem<-c(rem,i)
      } else {
      system(paste0('mkdir -p ',opt$output,'/',names(betas)[i]))
      for (chr in 1:22){
        fwrite(betas[betas$CHR==chr,c('SNPID','CHR','POS','A1','A2',names(betas)[i]), with=F], 
                paste0(opt$output,'/',names(betas)[i],'/chr',chr,'_',gsub('LDPred2_','',names(betas)[i]),'.weight'), quote=F, sep='\t', col.names=F, na='NA')
      }}
}

cat('.\n.\n.\n')
cat(paste0('Betas Saved in ',opt$output,'\n'))
