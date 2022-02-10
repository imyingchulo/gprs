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
filter_condition <- args[4]
snps_nb <- args[5]
print(paste("score file= ",scoref),quote=F)
print(paste("pheno file= ",phenof),quote=F)
print(paste("input_file_name= ",input_file_name),quote=F)
print(paste("filter_condition= ",filter_condition),quote=F)
print(paste("snps_nb= ",snps_nb),quote=F)
cat("\n")

##linear regression
pheno <- read.table(phenof,header=T)
pheno$pheno <- as.factor(pheno$pheno)
score <- read.table(scoref,header = T)
prs <- inner_join(score[,c(1,4)], pheno, by="id")
logit <- glm(pheno~., data=prs[,-c(1)], family="binomial")

prs.degree_of_freedom <-summary(logit)$df[2]
prs.coef <- summary(logit)$coeff["SCORE_SUM",]
prs.beta <- as.numeric(prs.coef[1])
prs.aic <- as.numeric(summary(logit)$aic)
prs.p <- as.numeric(prs.coef[4])
prs.r2 <- as.numeric(PseudoR2(logit,which="Nagelkerke"))
# summary(logit)

# print(paste("this is degree_of_freedom:",prs.degree_of_freedom),quote=F)

#odds ratio
prs$decile <- ntile(prs$SCORE_SUM,100)
prs$decile_top2 <- ntile(prs$SCORE_SUM,100)
prs$decile_top5 <- ntile(prs$SCORE_SUM,100)
prs$decile_top10 <- ntile(prs$SCORE_SUM,100)

prs$decile_top2[prs$decile_top2 %in% c(99:100)] <- "top2"
prs$decile_top5[prs$decile_top5 %in% c(96:100)] <- "top5"
prs$decile_top10[prs$decile_top10 %in% c(91:100)] <- "top10"
# prs$decile_top2[prs$decile_top2 %in% c(1:2)] <- "top2"
# prs$decile_top5[prs$decile_top5 %in% c(1:5)] <- "top5"
# prs$decile_top10[prs$decile_top10 %in% c(1:10)] <- "top10"

prs$decile[prs$decile %in% c(40:60)] <- "middle_20"
prs$decile_top2[prs$decile_top2 %in% c(40:60)] <- "middle_20"
prs$decile_top5[prs$decile_top5 %in% c(40:60)] <- "middle_20"
prs$decile_top10[prs$decile_top10 %in% c(40:60)] <- "middle_20"

# output_name <- as.character(args[6])
# filename <- paste(output_name, filter_pvalue, "stat.indv.txt", sep="_")
# test_user = c(filter(prs,prs$decile %in% c(1,,"middle_20"),
#               filter(prs,prs$decile_top2 %in% c("top2"),
#               filter(prs,prs$decile_top5 %in% c("top5"),
#               filter(prs,prs$decile_top10 %in% c("top10"))
# write.csv(prs, file = filename, row.names = FALSE)

user1 = filter(prs,prs$decile %in% c(100,"middle_20"))
# user1 = filter(prs,prs$decile %in% c(1,"middle_20"))
user2 = filter(prs,prs$decile_top2 %in% c("top2","middle_20"))
user3 = filter(prs,prs$decile_top5 %in% c("top5","middle_20"))
user4 = filter(prs,prs$decile_top10 %in% c("top10","middle_20"))

cont.table1 <- table(user1$decile,user1$pheno)
cont.table2 <- table(user2$decile_top2,user2$pheno)
cont.table3 <- table(user3$decile_top5,user3$pheno)
cont.table4 <- table(user4$decile_top10,user4$pheno)

names(dimnames(cont.table1)) <- c("decile","pheno")
names(dimnames(cont.table2)) <- c("decile","pheno")
names(dimnames(cont.table3)) <- c("decile","pheno")
names(dimnames(cont.table4)) <- c("decile","pheno")

print(cont.table1)
print(cont.table2)
print(cont.table3)
print(cont.table4)

OR1<-OddsRatio(cont.table1)
OR2<-OddsRatio(cont.table2)
OR3<-OddsRatio(cont.table3)
OR4<-OddsRatio(cont.table4)

# head(OR1)
# head(OR2)
# head(OR3)
# head(OR4)


#ROC curve
myroc <-roc(prs$pheno, prs$SCORE_SUM, auc=TRUE, plot=TRUE,auc.polygon=TRUE,print.auc=TRUE)

stat <- data.frame(data=input_file_name,filter_condition=filter_condition, snps_nb=snps_nb, P=prs.p,
                   BETA=prs.beta, Degree_of_freedom=prs.degree_of_freedom ,AIC=prs.aic, AUC=myroc$auc, PseudoR2=prs.r2,
                   OR_top1_to_middle20=OR1,OR_top2_to_middle20=OR2,OR_top5_to_middle20=OR3,OR_top10_to_middle20=OR4)

new_stat <- as.data.frame(stat)
output_name <- as.character(args[6])
filename <- paste(output_name, filter_pvalue, "stat.txt", sep="_")
write.csv(new_stat, file = filename, row.names = FALSE)



