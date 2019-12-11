rm(list = ls())
fileName1 <- '/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/summary_major_clade_fixed_alleles.txt'
fileName2 <- '/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/summary_basal_clade_fixed_alleles.txt'


myDataFrame <- read.table(fileName,header=T)
mySubDF <- myDataFrame


#data cleaning by project
mySubDF <- subset(mySubDF, pop_type != "del_srlD" & pop_type != "del_mutL_srlD")
mySubDF <- cbind(mySubDF, pop_100 = floor(mySubDF$pop/100))

#categorization of types of population
mySubDF <- cbind(mySubDF, pop_type2 = mySubDF$pop_type)
mySubDF[mySubDF$pop_type2 == "del_araBAD", "pop_type2"] <- "wt"
mySubDF[mySubDF$pop_type2 == "del_mutL_araBAD", "pop_type2"] <- "del_mutL"

mySubDF <- cbind(mySubDF, pop_type3 = mySubDF$pop_type2)
mySubDF <- arrange(transform(mySubDF,pop_type=factor(pop_type3,levels=c("del_mutL","wt"))),pop_type3)
mySubDF$pop_type3 <- factor(mySubDF$pop_type3, labels = c("MMR-{}", "WT"))


mySubDF <- subset(mySubDF, mut_type=="nonsynonymous" & pop_100==5)

mySelectGeneList <- as.data.frame(unique(mySubDF$gene_name))
write.table(mySelectGeneList,'/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/select_gene.txt',col.names=FALSE,row.names=FALSE,sep="\t", quote = FALSE)

