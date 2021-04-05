library(dplyr)
setwd("/Users/xilu/Desktop/T2D/T2D13/final")
name1 = read.csv(file="name1.csv", header=TRUE)
name2 = read.csv(file="name2.csv", header=TRUE)

#name11 = name1[which(id1 %in% id2),2:3]
gene_name1 = merge(name1, name2, by="dbSNP_RS_ID")
#gene_name1 = inner_join(name1, name2, by="dbSNP_RS_ID")

gene_names = cbind(gene_name1[,1:3],gene_name1[,11:12])

write.csv(gene_names, file = "gene_names.csv")

### LADBLSS ###
setwd("/Users/xilu/Desktop/T2D/T2D13")
ladblss = read.csv(file="LADBLSS_prop.csv", header=TRUE)

library(stringr)
regexp <- "[[:digit:]]+"

id = str_extract(ladblss$X, regexp)
ladblss2 =cbind(id,ladblss[-1])

id = str_extract(gene_names$Probe_Set_ID, regexp)
gene_name3 = cbind(id, gene_names[,5],gene_names[,1])
colnames(gene_name3) = c("id","gene_name","rsID")
gene_name3 = as.data.frame(gene_name3)
length(which(ladblss2$id %in% gene_name3$id))

LADBLSS_update = merge(gene_name3, ladblss2, by="id")

write.csv(LADBLSS_update[-1], file = "LADBLSS_prop_update.csv")


### BLSS ###
setwd("/Users/xilu/Desktop/T2D/T2D13")
blss = read.csv(file="BLSS_prop.csv", header=TRUE)

library(stringr)
regexp <- "[[:digit:]]+"

id = str_extract(blss$X, regexp)
blss2 =cbind(id,blss[-1])

id = str_extract(gene_names$Probe_Set_ID, regexp)
gene_name3 = cbind(id, gene_names[,5],gene_names[,1])
colnames(gene_name3) = c("id","gene_name","rsID")
gene_name3 = as.data.frame(gene_name3)
length(which(blss2$id %in% gene_name3$id))

BLSS_update = merge(gene_name3, blss2, by="id")

write.csv(BLSS_update[-1], file = "BLSS_prop_update.csv")


## BL ##
setwd("/Users/xilu/Desktop/T2D/T2D13")
load("Data101.RData")
load("top_genes.RData")
load("T2D_BL.RData")
index_BL = index_BL[1:100]

n=nrow(top2000_gene); p=ncol(top2000_gene);q2=5
y=as.matrix(clinic$bmi)
y[which(is.na(y))]=mean(y[-which(is.na(y))])

e=as.matrix(cbind(clinic$age,clinic$act,clinic$trans,clinic$ceraf,clinic$chol))
for (i in 1:5) 
{
  e1=e[,i]
  e1[which(is.na(e1))]=mean(e1[-which(is.na(e1))])
  e[,i]=e1
}

colnames(e)=c("age","act","trans","ceraf","chol")
e=scale(e)

g=as.matrix(top2000_gene)
g=scale(g)
gene_name = colnames(g)
e_name = colnames(e)

coef_total = rep(0,(p+p*q2))

prop = m2[1:100,1]
coef_total[index_BL] = prop

coef_main_2 = coef_total[1:p]
coef_inter_2 = coef_total[(p+1):(p+p*q2)]
matrix_inter = matrix(coef_inter_2,nrow = p, ncol = q2, byrow = TRUE)
matrix_total = cbind(coef_main_2, matrix_inter)
rownames(matrix_total) = gene_name
colnames(matrix_total) = c("main", e_name)

main_BL = which(rowSums(matrix_total)!=0)
result_total = matrix_total[rowSums(matrix_total)!=0,]
View(result_total)
write.csv(result_total, file = "BL_prop.csv")


setwd("/Users/xilu/Desktop/T2D/T2D13")
bl = read.csv(file="BL_prop.csv", header=TRUE)

library(stringr)
regexp <- "[[:digit:]]+"

id = str_extract(bl$X, regexp)
bl2 =cbind(id,bl[-1])

id = str_extract(gene_names$Probe_Set_ID, regexp)
gene_name3 = cbind(id, gene_names[,5],gene_names[,1])
colnames(gene_name3) = c("id","gene_name","rsID")
gene_name3 = as.data.frame(gene_name3)
length(which(bl2$id %in% gene_name3$id))

BL_update = merge(gene_name3, bl2, by="id")

write.csv(BL_update[-1], file = "BL_prop_update.csv")


## LADBL ##
setwd("/Users/xilu/Desktop/T2D/T2D13")
load("Data101.RData")
load("top_genes.RData")
load("T2D_LADBL.RData")
index_LADBL = index_LADBL[1:100]

n=nrow(top2000_gene); p=ncol(top2000_gene);q2=5
y=as.matrix(clinic$bmi)
y[which(is.na(y))]=mean(y[-which(is.na(y))])

e=as.matrix(cbind(clinic$age,clinic$act,clinic$trans,clinic$ceraf,clinic$chol))
for (i in 1:5) 
{
  e1=e[,i]
  e1[which(is.na(e1))]=mean(e1[-which(is.na(e1))])
  e[,i]=e1
}

colnames(e)=c("age","act","trans","ceraf","chol")
e=scale(e)

g=as.matrix(top2000_gene)
g=scale(g)
gene_name = colnames(g)
e_name = colnames(e)


coef_total = rep(0,(p+p*q2))

prop = m2[1:100,1]
coef_total[index_LADBL] = prop

coef_main_2 = coef_total[1:p]
coef_inter_2 = coef_total[(p+1):(p+p*q2)]
matrix_inter = matrix(coef_inter_2,nrow = p, ncol = q2, byrow = TRUE)
matrix_total = cbind(coef_main_2, matrix_inter)
rownames(matrix_total) = gene_name
colnames(matrix_total) = c("main", e_name)

main_LADBL = which(rowSums(matrix_total)!=0)
result_total = matrix_total[rowSums(matrix_total)!=0,]
View(result_total)
write.csv(result_total, file = "LADBL_prop.csv")


setwd("/Users/xilu/Desktop/T2D/T2D13")
ladbl = read.csv(file="LADBL_prop.csv", header=TRUE)

library(stringr)
regexp <- "[[:digit:]]+"

id = str_extract(ladbl$X, regexp)
ladbl2 =cbind(id,ladbl[-1])

id = str_extract(gene_names$Probe_Set_ID, regexp)
gene_name3 = cbind(id, gene_names[,5],gene_names[,1])
colnames(gene_name3) = c("id","gene_name","rsID")
gene_name3 = as.data.frame(gene_name3)
length(which(ladbl2$id %in% gene_name3$id))

LADBL_update = merge(gene_name3, ladbl2, by="id")

write.csv(LADBL_update[-1], file = "LADBL_prop_update.csv")









