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
output_name <- as.character(args[6])
pcs <- args[7]

print(paste("score file= ",scoref),quote=F)
print(paste("pheno file= ",phenof),quote=F)
print(paste("input_file_name= ",input_file_name),quote=F)
print(paste("filter_conditions= ",filter_condition),quote=F)
print(paste("snps_nb= ",snps_nb),quote=F)
print(paste("output_name= ",output_name),quote=F)
print(paste("pcs= ",pcs),quote=F)
cat("\n")

#read and merge the files
pheno <- read.table(phenof,header=T)
# pheno$pheno <- as.factor(pheno$pheno)
score <- read.table(scoref,header = T)
pcf <- read.table(pcs,header=F)
colnames(pcf) <- c("id", paste0("PC",1:10))
sscore_pc <- merge(score, pcf, by='id')
pheno.prs <- inner_join(sscore_pc[,c(1,3:13)], pheno, by="id")
pheno.prs$pheno <- as.factor(pheno.prs$pheno)
logit <- glm(pheno~., data=pheno.prs[,-c(1)], family="binomial", na.action = na.pass)

pheno.prs.degree_of_freedom <-summary(logit)$df[2]
pheno.prs.coef <- summary(logit)$coeff["SCORE_SUM",]
pheno.prs.beta <- as.numeric(pheno.prs.coef[1])
pheno.prs.aic <- as.numeric(summary(logit)$aic)
pheno.prs.p <- as.numeric(pheno.prs.coef[4])
pheno.prs.r2 <- as.numeric(PseudoR2(logit,which="Nagelkerke"))

#odds ratio
pheno.prs$decile<- ntile(desc(pheno.prs$SCORE_SUM),100)
pheno.prs$decile_top2 <- ntile(desc(pheno.prs$SCORE_SUM),100)
pheno.prs$decile_top5 <- ntile(desc(pheno.prs$SCORE_SUM),100)
pheno.prs$decile_top10 <- ntile(desc(pheno.prs$SCORE_SUM),100)

pheno.prs$decile_top2[pheno.prs$decile_top2 %in% c(99:100)] <- "top2"
pheno.prs$decile_top5[pheno.prs$decile_top5 %in% c(95:100)] <- "top5"
pheno.prs$decile_top10[pheno.prs$decile_top10 %in% c(90:100)] <- "top10"

pheno.prs$decile[pheno.prs$decile %in% c(40:60)] <- "middle_20"
pheno.prs$decile_top2[pheno.prs$decile_top2 %in% c(40:60)] <- "middle_20"
pheno.prs$decile_top5[pheno.prs$decile_top5 %in% c(40:60)] <- "middle_20"
pheno.prs$decile_top10[pheno.prs$decile_top10 %in% c(40:60)] <- "middle_20"

user1 = filter(pheno.prs,pheno.prs$decile %in% c(100,"middle_20"))
user2 = filter(pheno.prs,pheno.prs$decile_top2 %in% c("top2","middle_20"))
user3 = filter(pheno.prs,pheno.prs$decile_top5 %in% c("top5","middle_20"))
user4 = filter(pheno.prs,pheno.prs$decile_top10 %in% c("top10","middle_20"))

cont.table1 <- table(user1$decile,user1$pheno)
cont.table2 <- table(user2$decile_top2,user2$pheno)
cont.table3 <- table(user3$decile_top5,user3$pheno)
cont.table4 <- table(user4$decile_top10,user4$pheno)

table(user1$decile,user1$pheno)
names(dimnames(cont.table1)) <- c("decile","pheno")
names(dimnames(cont.table2)) <- c("decile","pheno")
names(dimnames(cont.table3)) <- c("decile","pheno")
names(dimnames(cont.table4)) <- c("decile","pheno")

# print(cont.table1)
# print(cont.table2)
# print(cont.table3)
# print(cont.table4)

OR1<-OddsRatio(cont.table1)
OR2<-OddsRatio(cont.table2)
OR3<-OddsRatio(cont.table3)
OR4<-OddsRatio(cont.table4)

#ROC curve
myroc <-roc(pheno.prs$pheno, pheno.prs$SCORE_SUM, auc=TRUE, plot=TRUE,auc.polygon=TRUE,print.auc=TRUE)

stat <- data.frame(data=input_file_name,filter_condition=filter_condition, snps_nb=snps_nb,
                   P=pheno.prs.p,
                   BETA=pheno.prs.beta, Degree_of_freedom=pheno.prs.degree_of_freedom ,
                   AIC=pheno.prs.aic, AUC=myroc$auc,
                   PseudoR2=pheno.prs.r2,
                   OR_top1_to_middle20=OR1,OR_top2_to_middle20=OR2,OR_top5_to_middle20=OR3,OR_top10_to_middle20=OR4)

new_stat <- as.data.frame(stat)
filename <- paste(output_name, filter_condition, "stat.txt", sep="_")
write.csv(stat, file = filename, row.names = FALSE)


