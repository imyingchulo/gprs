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
score <- read.table(scoref,header = T)
pcf <- read.table(pcs,header=F)
colnames(pcf) <- c("id", paste0("PC",1:10))
sscore_pc <- merge(score, pcf, by='id')
pheno.prs <- inner_join(sscore_pc[,c(1,3:13)], pheno, by="id")
pheno.prs$pheno <- as.factor(pheno.prs$pheno)

# logistic models
full_logit <- glm(pheno~., data=pheno.prs[,-c(1)], family="binomial", na.action = na.pass)
partial_logit <- glm(pheno~., data=pheno.prs[,-c(1,2)], family="binomial", na.action = na.pass)

# calculate r2
pheno.prs.r2 <- as.numeric(PseudoR2(full_logit,which="Nagelkerke"))
pheno.prs.part.r2 <- as.numeric(PseudoR2(partial_logit,which="Nagelkerke"))

# get sscore r2
sscore.r2 <- (pheno.prs.r2-pheno.prs.part.r2)

# output data
stat <- data.frame(data=input_file_name,
                   Filter_condition=filter_condition,
                   Fullmodel_R2=pheno.prs.r2,
                   Partial_R2=pheno.prs.part.r2,
                   PseudoR2_sscore=sscore.r2)

new_stat <- as.data.frame(stat)
filename <- paste(output_name, filter_condition, "stat.txt", sep="_")
write.csv(stat, file = filename, row.names = FALSE)