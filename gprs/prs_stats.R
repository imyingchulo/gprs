##USAGE Rcript --vanilla PRS.R [score file] [pheno file] [graph pdf name]
#pheno file requires header: 1st "id" 2nd "pheno" coded 0/1, 3rd - any additional covariates to include
#score file needs header: 1st"id" ...4th "SCORE_SUM" (sscore.py output works)

# install.packages(c("tidyverse","pROC","DescTools","WebPower"))

library(tidyverse)
library(pROC)
library(DescTools)
library(WebPower)

args <- commandArgs(TRUE)

scoref <- args[1]
phenof <- args[2]
input_file_name <- args[3]
filter_pvalue <- as.numeric(args[4])
print(paste("score file= ",scoref),quote=F)
print(paste("pheno file= ",phenof),quote=F)
print(paste("input_file_name= ",input_file_name),quote=F)
print(paste("filter_pvalue= ",filter_pvalue),quote=F)
cat("\n")

##logistic regression
pheno <- read.table(phenof,header=T)
pheno$pheno <- as.factor(pheno$pheno)
score <- read.table(scoref,header = T)
prs <- inner_join(score[,c(1,4)], pheno, by="id")
logit <- glm(pheno~., data=prs[,-c(1)], family="binomial")


prs.coef <- summary(logit)$coeff["SCORE_SUM",]
prs.beta <- as.numeric(prs.coef[1])
prs.aic <- as.numeric(summary(logit)$aic)
prs.p <- as.numeric(prs.coef[4])
prs.r2 <- as.numeric(PseudoR2(logit,which="Nagelkerke"))
summary(logit)

# #odds ratio
prs$decile <- ntile(prs$SCORE_SUM,10)
user1 <- filter(prs,prs$decile %in% c(1,5))
user2 <- filter(prs,prs$decile %in% c(2,5))
user3 <- filter(prs,prs$decile %in% c(3,5))
user4 <- filter(prs,prs$decile %in% c(4,5))
user5 <- filter(prs,prs$decile %in% c(6,5))
user6 <- filter(prs,prs$decile %in% c(7,5))
user7 <- filter(prs,prs$decile %in% c(8,5))
user8 <- filter(prs,prs$decile %in% c(9,5))
user9 <- filter(prs,prs$decile %in% c(10,5))

cont.table1 <- table(user1$decile,user1$pheno)
cont.table2 <- table(user2$decile,user2$pheno)
cont.table3 <- table(user3$decile,user3$pheno)
cont.table4 <- table(user4$decile,user4$pheno)
cont.table5 <- table(user5$decile,user5$pheno)
cont.table6 <- table(user6$decile,user6$pheno)
cont.table7 <- table(user7$decile,user7$pheno)
cont.table8 <- table(user8$decile,user8$pheno)
cont.table9 <- table(user9$decile,user9$pheno)

names(dimnames(cont.table1)) <- c("decile","pheno")
names(dimnames(cont.table2)) <- c("decile","pheno")
names(dimnames(cont.table3)) <- c("decile","pheno")
names(dimnames(cont.table4)) <- c("decile","pheno")
names(dimnames(cont.table5)) <- c("decile","pheno")
names(dimnames(cont.table6)) <- c("decile","pheno")
names(dimnames(cont.table7)) <- c("decile","pheno")
names(dimnames(cont.table8)) <- c("decile","pheno")
names(dimnames(cont.table9)) <- c("decile","pheno")


OR1<-OddsRatio(cont.table1)
OR2<-OddsRatio(cont.table2)
OR3<-OddsRatio(cont.table3)
OR4<-OddsRatio(cont.table4)
OR5<-OddsRatio(cont.table5)
OR6<-OddsRatio(cont.table6)
OR7<-OddsRatio(cont.table7)
OR8<-OddsRatio(cont.table8)
OR9<-OddsRatio(cont.table9)

# head(OR1) 
# head(OR2)
# head(OR3)
# head(OR4)
# head(OR5)
# head(OR6)
# head(OR7)
# head(OR8)
# head(OR9)



#ROC curve
myroc <-roc(prs$pheno, prs$SCORE_SUM, auc=TRUE, plot=TRUE,auc.polygon=TRUE,print.auc=TRUE)
stat <- data.frame(data=input_file_name, filter_pvalue=filter_pvalue, P=prs.p,
                   BETA=prs.beta, AIC=prs.aic, AUC=myroc$auc, PseudoR2=prs.r2,
                   OR1vs5=OR1,OR2vs5=OR2,OR3vs5=OR3,OR4vs5=OR4,OR6vs5=OR5,
                   OR7vs5=OR6,OR8vs5=OR7,OR9vs5=OR8,OR10vs5=OR9)
stat <- as.data.frame(stat)
output_name <- as.character(args[5])
filename <- paste(output_name, filter_pvalue, "stat.txt", sep="_")
write.csv(stat, file = filename, row.names = FALSE)

