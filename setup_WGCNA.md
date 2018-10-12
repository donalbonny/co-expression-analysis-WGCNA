### Settup WGCNA libray in R
#I have  trouble when I first install WGCNA package in R since R required a lot of packages

source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 

library(MASS)
library(cluster)
library(impute) # install it for imputing missing value 
library(WGCNA)
options(stringsAsFactors = F)



require(data.table)
data<-as.data.frame(fread("transcript_tpms_all_samples.tsv"))


head(data)

heatmap(cor(data), trace='none', main='Sample correlations (raw)')

set.seed(2018)

#Manipulate file so it matches the format WGCNA needs


row.names(data) = data$target_id
data$target_id = NULL
data = as.data.frame(t(data))   # WGCNA requires genes as columns 

dim(data)


# 59 samples and 48359 genes 

#Run this to check if there are gene outliers

gsg = goodSamplesGenes(data, verbose = 3)
gsg$allOK

# if not, can check this 

if (!gsg$allOK)
  {if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
    datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

gsg$allOK


powers = c(c(1:10), seq(from =10, to=30, by=1))
sft = pickSoftThreshold(data, powerVector=powers, verbose =5, networkType="signed")
