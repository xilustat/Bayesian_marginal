setwd("/Users/xilu/Desktop/T2D/T2D13/final")
BL = read.csv(file="BL_update.csv", header=TRUE)
BL_gene = BL$gene_name
BLSS = read.csv(file="BLSS_update.csv", header=TRUE)
BLSS_gene = BLSS$gene_name
LADBL = read.csv(file="LADBL_update.csv", header=TRUE)
LADBL_gene = LADBL$gene_name
LADBLSS = read.csv(file="LADBLSS_update.csv", header=TRUE)
LADBLSS_gene = LADBLSS$gene_name

#unique gene in LADBLSS
#method 1
gene_41 = setdiff(LADBLSS_gene, BL_gene)
gene_42 = setdiff(LADBLSS_gene, BLSS_gene)
gene_43 = setdiff(LADBLSS_gene, LADBL_gene)

inter_412 = intersect(gene_41, gene_42)
inter_4123 = intersect(inter_412, gene_43)

#method 2

gene_123 = unique(c(BL_gene, BLSS_gene, LADBL_gene))
gene_4123 = setdiff(LADBLSS_gene, gene_123)


#[1] "ATP5C1"        "MIR5100"       "ARHGAP22"      "CHAT"          "MALRD1"       
#[6] "RP11-478H13.1" "RP11-192P3.5"  "PLXDC2"        "RP11-810B23.1"


#gene_name T2D google scholar

# A = c("Dog", "Cat", "Mouse")
# B = c("Tiger","Lion","Cat")
# A %in% B
# FALSE  TRUE FALSE
# intersect(A,B)
# "Cat"
# setdiff(A,B)
# "Dog"   "Mouse"
# setdiff(B,A)
# "Tiger" "Lion" 


## CHIAP2 SFTA1P CEACAM8 
