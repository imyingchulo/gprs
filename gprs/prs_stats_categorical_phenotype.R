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

##linear regression
pheno <- read.table(phenof,header=T)
pheno$pheno <- as.factor(pheno$pheno)
score <- read.table(scoref,header = T)
prs <- inner_join(score[,c(1,4)], pheno, by="id")
logit <- glm(pheno~., data=prs[,-c(1)], family="normal")


prs.coef <- summary(logit)$coeff["SCORE_SUM",]
prs.beta <- as.numeric(prs.coef[1])
prs.aic <- as.numeric(summary(logit)$aic)
prs.p <- as.numeric(prs.coef[4])
prs.r2 <- as.numeric(PseudoR2(logit,which="Nagelkerke"))
summary(logit)

# #odds ratio
prs$decile <- ntile(prs$SCORE_SUM,100)
prs$decile[prs$decile == 40] <- "middle_20"
prs$decile[prs$decile == 41] <- "middle_20"
prs$decile[prs$decile == 42] <- "middle_20"
prs$decile[prs$decile == 43] <- "middle_20"
prs$decile[prs$decile == 44] <- "middle_20"
prs$decile[prs$decile == 45] <- "middle_20"
prs$decile[prs$decile == 46] <- "middle_20"
prs$decile[prs$decile == 47] <- "middle_20"
prs$decile[prs$decile == 48] <- "middle_20"
prs$decile[prs$decile == 49] <- "middle_20"
prs$decile[prs$decile == 50] <- "middle_20"
prs$decile[prs$decile == 51] <- "middle_20"
prs$decile[prs$decile == 52] <- "middle_20"
prs$decile[prs$decile == 53] <- "middle_20"
prs$decile[prs$decile == 54] <- "middle_20"
prs$decile[prs$decile == 55] <- "middle_20"
prs$decile[prs$decile == 56] <- "middle_20"
prs$decile[prs$decile == 57] <- "middle_20"
prs$decile[prs$decile == 58] <- "middle_20"
prs$decile[prs$decile == 59] <- "middle_20"
prs$decile[prs$decile == 60] <- "middle_20"


user1 <- filter(prs,prs$decile %in% c(1,middle20))
user2 <- filter(prs,prs$decile %in% c(2,middle20))
user3 <- filter(prs,prs$decile %in% c(5,middle20))
user4 <- filter(prs,prs$decile %in% c(10,middle20))


cont.table1 <- table(user1$decile,user1$pheno)
cont.table2 <- table(user2$decile,user2$pheno)
cont.table3 <- table(user3$decile,user3$pheno)
cont.table4 <- table(user4$decile,user4$pheno)


names(dimnames(cont.table1)) <- c("decile","pheno")
names(dimnames(cont.table2)) <- c("decile","pheno")
names(dimnames(cont.table3)) <- c("decile","pheno")
names(dimnames(cont.table4)) <- c("decile","pheno")


OR1<-OddsRatio(cont.table1)
OR2<-OddsRatio(cont.table2)
OR3<-OddsRatio(cont.table3)
OR4<-OddsRatio(cont.table4)

# head(OR1)
# head(OR2)
# head(OR3)
# head(OR4)

stat <- as.data.frame(stat)
output_name <- as.character(args[5])
filename <- paste(output_name, filter_pvalue, "stat.txt", sep="_")
write.csv(stat, file = filename, row.names = FALSE)


