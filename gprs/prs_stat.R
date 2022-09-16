#!/usr/bin/Rscript
#created by - Soyoung Jeon
suppressWarnings(library(tidyverse))
suppressMessages(library(pROC))
suppressMessages(library(DescTools))
suppressMessages(library(WebPower))

# options
args <- commandArgs(TRUE)

model_name <- args[1]
scoref <- args[2] # doesn't need header, but first column ID; second column SCORE; plink2 output works
phenof <- args[3] # need to have header, first column ID; second column PHENO; the rest covariates
family <- args[4] # binary or quantitative
pop_prev <- args[5]
plotroc <- args[6] #plotroc OR no_plot
output_name <- args[7] #full output directory + model_name

cat("Specified Options: \n")
cat(paste("model_name = ",model_name,"\n"))
cat(paste("score file = ",scoref,"\n"))
cat(paste("pheno file = ",phenof,"\n"))
cat(paste("phenotype = ",family,"\n"))
cat(paste("population prevalence = ",family,"\n"))
cat(paste("output_name = ",output_name,"\n"))
cat("\n")

# function for converting R2 into liability R2
h2l_R2 <- function(k, r2, p) {
  # K baseline disease risk
  # r2 from a linear regression model attributable to genomic profile risk score
  # P proportion of sample that are cases
  # calculates proportion of variance explained on the liability scale
  #from ABC at http://www.complextraitgenomics.com/software/
  #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
  x= qnorm(1-k)
  z= dnorm(x)
  i=z/k
  C= k*(1-k)*k*(1-k)/(z^2*p*(1-p))
  theta= i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
  h2l_R2 = C*r2 / (1 + C*theta*r2)
}


###########
# Read in prs
###########
score <- read.table(scoref)
names(score) <- c("ID","SCORE_SUM","TOTAL_ALLELE_CT")


###########
# Read in the phenotype data
###########

pheno<-read.table(phenof, header = TRUE)
# Check first two columns are named right
if (names(pheno)[1] != "ID" | names(pheno)[2] != "PHENO"){
    stop('Header names for Phenotype file are wrong. Please read manual and format accordingly.\n')}

cat('Phenotype file contains',dim(pheno)[1],'individuals and ',dim(pheno)[2]-2,'covariates.\n')
cat("\n")

# Determine whether outcome is binary or continuous and formatted correctly
if(family == 'binary'){
  if( length(unique(pheno[,2])) > 2 ){
    stop('Phenotype has more than two values.\n')
  }
  if( pop_prev == 'NA'){
    stop('Population disease prevalence for calculating liability r2 is not specified')
  } else {
    pop_prev <- as.numeric(pop_prev)
  }
}

if(length(unique(pheno[,2])) == 2 & family == 'quantitative'){
    warning('Phenotype has only two values. \n')
}


###########
# Merge the phenotype and prs
###########
prs <- inner_join(score[,c(1,2)], pheno, by="ID")

if(family == 'binary'){
    prs$PHENO<-factor(prs$PHENO, labels=c('CONTROL','CASE'))
    logit <- glm(PHENO~., data=prs[,-c(1)], family="binomial")
    prs.coef <- summary(logit)$coeff["SCORE_SUM",]
    prs.pseudor2 <- as.numeric(PseudoR2(logit,which="Nagelkerke"))
    prs.obs_r2<-cor(predict(logit), as.numeric(prs$PHENO))^2
    prs.leesr2 <- h2l_R2(pop_prev, prs.obs_r2, sum(prs$PHENO== 'CASE')/length(prs$PHENO))
    myroc <-roc(prs$PHENO, prs$SCORE_SUM, auc=TRUE, quiet=TRUE)
    if(plotroc == 'plotroc'){
      pdf(paste0(output_name,".pdf"))
      plot.roc(myroc,auc.polygon=TRUE, print.auc=TRUE)
      dev.off()
      cat(paste0('ROC plot with AUC saved: ',output_name,'.pdf\n\n'))
    }
    # OR by decile
    prs$decile <- ntile(prs$SCORE_SUM,100)
    prs$decile[prs$decile %in% c(40:60)] <- "middle_20"
    prs$decile[prs$decile == 100] <- "top_1"
    prs$decile[prs$decile  %in% c(99:100)] <- "top_2"
    prs$decile[prs$decile %in% c(96:100)] <- "top_5"
    prs$decile[prs$decile %in% c(91:100)] <- "top_10"
    
    user1 = filter(prs, prs$decile %in% c("top_1", "middle_20"))
    user2 = filter(prs, prs$decile %in% c("top_2","middle_20"))
    user3 = filter(prs, prs$decile %in% c("top_5","middle_20"))
    user4 = filter(prs, prs$decile %in% c("top_10","middle_20"))
    
    cat('Contingency table per deciles: \n')
    print(cont.table1 <- table(user1$decile,user1$PHENO))
    print(cont.table2 <- table(user2$decile,user2$PHENO))
    print(cont.table3 <- table(user3$decile,user3$PHENO))
    print(cont.table4 <- table(user4$decile,user4$PHENO))
    
    OR1<-OddsRatio(cont.table1)
    OR2<-OddsRatio(cont.table2)
    OR3<-OddsRatio(cont.table3)
    OR4<-OddsRatio(cont.table4)
    
    # result dataframe
    stat <- data.frame( Model = model_name, P=prs.coef[4],
		            Beta=prs.coef[1], SE=prs.coef[2], OR=exp(prs.coef[1]), 
                AUC=myroc$auc, PseudoR2=prs.pseudor2, LiabilityR2=prs.leesr2,
		            N=length(prs$PHENO), N_cas=sum(prs$PHENO=='CASE'), N_ctrl=sum(prs$PHENO=='CONTROL'),  
                OR_top1_to_middle20=OR1,OR_top2_to_middle20=OR2,OR_top5_to_middle20=OR3,OR_top10_to_middle20=OR4)

    
    } else {
      lnr <- glm(PHENO~., data=prs[,-c(1)], family="gaussian")
      prs.coef <- summary(lnr)$coeff["SCORE_SUM",]
      prs.obs_r2<-cor(predict(lnr), as.numeric(prs$PHENO))^2
      stat <- data.frame( Model = model_name, P=prs.coef[4],
                          Beta=prs.coef[1], SE=prs.coef[2], R2=prs.obs_r2, N=length(prs$PHENO) )
}

write.table( format(stat, digits=3), paste0(output_name, ".stat"), row.names = F, quote = F, sep=" ")
cat("\n")
cat(paste0("Statistics calculation is done. Results saved as ",output_name,".stat \n"))