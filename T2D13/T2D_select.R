library(roben)
library(mnormt)
library(MASS)
library(MCMCpack)
library(expm)
library(rmutil)
library(stringr)

setwd("/Users/xilu/Desktop/T2D/T2D13")
load("Data101.RData")

n=nrow(gene)
p=ncol(gene)
q2=5

y=as.matrix(clinic$bmi)
e=as.matrix(cbind(clinic$age,clinic$act,clinic$trans,clinic$ceraf,clinic$chol))
colnames(e)=c("age","act","trans","ceraf","chol")
e=scale(e)


g=as.matrix(gene)
g=scale(g)
gene_name = colnames(gene)
e_name = colnames(e)

p_number = c()

for (i in 1:p) {
  
  x=g[,i]
  fit = lm(y~x + x:e)
  pp = summary(fit)$coefficients[,4][2:7] 
  min_pp = sum(pp<=0.02)
  #p_value = c(p_value, pp)
  p_number = c(p_number,min_pp)
}

index = c(1:p)

m1 = cbind(index, p_number)
m1 = na.omit(m1)
m2 = m1[order(m1[,2], decreasing = TRUE),]

index_gene = m2[,1]   

top2000_index = index_gene[1:2000]   #take first 2000
top2000_gene = g[,top2000_index]  # first 2000 genes

save(top2000_gene, top2000_index, index_gene, file = "top_genes.RData")  

save(top2000_gene, top2000_index, index_gene,m2,file = "top_genes_m2.RData") 



e_rest = clinic[,-c(1:4,6,11,15,26)]
e_rest_name = colnames(e_rest)

pvalue_e = c()
for (j in 1:ncol(e_rest)) {
  
  x=e_rest[,j]
  fit = lm(y~x)
  pp_e = summary(fit)$coefficients[,4][2] 
  
  pvalue_e = c(pvalue_e, pp_e)
}

m_e = cbind(e_rest_name, as.numeric(pvalue_e))
m_e = na.omit(m_e)
m_e = m_e[order(m_e[,2]),]

