# install.packages(c("tidyverse","pROC","DescTools","WebPower"))

library(tidyverse)
library(pROC)
library(DescTools)
library(WebPower)
library(rsq)

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
pheno$pheno <- as.factor(pheno$pheno)
score <- read.table(scoref,header = T)
pcf <- read.table(pcs,header=F)
colnames(pcf) <- c("id", paste0("PC",1:10))
sscore_pc <- merge(score, pcf, by='id')
pheno.prs <- inner_join(sscore_pc[,c(1,3,4:13)], pheno, by="id")

# logistic regression
prs.result <- NULL
pheno.prs$pheno <- as.numeric(pheno.prs$pheno)
pheno.prs$SCORE_SUM <- as.numeric(pheno.prs$SCORE_SUM)
lm_model <- lm(pheno~., data=pheno.prs[,-c(1)])
model <- summary(lm(pheno~., data=pheno.prs[,-c(1)]))
# model$coefficients
prs.df <- model$df[2]
model.r2 <- model$r.squared
prs.r2 <- model.r2
prs.coef <- model$coeff["SCORE_SUM",]
prs.p <- model$coeff["SCORE_SUM",4]
prs.result <- rbind(prs.result, c(R2=prs.r2,prs.coef[4], prs.coef[1], prs.coef[2], DF= prs.df))
partial.r2 <- rsq.partial(lm_model)

# output statistics result
stat <- as.data.frame(stat)
stat <- data.frame(data=input_file_name,filter_condition=filter_condition, snps_nb=snps_nb, P=prs.p,
                   Degree_of_freedom=prs.df , r2=prs.r2, Partial_r2=partial.r2$partial.rsq[1])
new_stat <- as.data.frame(stat)
filename <- paste(output_name, filter_condition, "stat.txt", sep="_")
write.csv(stat, file = filename, row.names = FALSE)





