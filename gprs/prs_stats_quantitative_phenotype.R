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
remove_column <- args[6]

print(paste("score file= ",scoref),quote=F)
print(paste("pheno file= ",phenof),quote=F)
print(paste("input_file_name= ",input_file_name),quote=F)
print(paste("filter_conditions= ",filter_conditions),quote=F)
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

stat <- as.data.frame(stat)
output_name <- as.character(args[5])

stat <- data.frame(data=input_file_name,filter_condition=filter_condition, snps_nb=snps_nb, P=prs.p,
                   BETA=prs.beta, Degree_of_freedom=prs.degree_of_freedom , PseudoR2=prs.r2)
new_stat <- as.data.frame(stat)
output_name <- as.character(args[6])
filename <- paste(output_name, filter_pvalue, "stat.txt", sep="_")
write.csv(stat, file = filename, row.names = FALSE)




