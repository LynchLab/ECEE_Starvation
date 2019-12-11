rm(list = ls())

library(dplyr) #arrange

myDataFrame1 <- read.table('/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/list_basal_fixed_alleles.txt', header = T, stringsAsFactors=F)
myDataFrame2 <- read.table('/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/list_major_minor_clade_fixed_alleles.txt', header = T, stringsAsFactors=F)
myDataFrame3 <- read.table('/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/list_basal_poly_alleles.txt', header = T, row.names=NULL)
myDataFrame <- rbind(myDataFrame1, myDataFrame2, myDataFrame3)
rm(myDataFrame2)
rm(myDataFrame3)

GeneList1 <- read.csv("/Users/wei-chinho/Documents/MURI2/01_ref_genome/GO/GO0006260.txt", sep='\t', stringsAsFactors=F)
GeneList2 <- read.csv("/Users/wei-chinho/Documents/MURI2/01_ref_genome/GO/GO0006281.txt", sep='\t', stringsAsFactors=F)
colnames(GeneList1) <- c("gene_name", "des");
colnames(GeneList2) <- c("gene_name", "des");
GeneList <- rbind(GeneList1,GeneList2)
GeneList$des <- NULL
GeneList <- unique(GeneList)

mySubDF <- myDataFrame
#data cleaning by project.
mySubDF <- subset(mySubDF, pop_type != "del_srlD" & pop_type != "del_mutL_srlD")
mySubDF <- cbind(mySubDF, pop_100 = floor(mySubDF$pop/100))

#categorization of types of population
mySubDF <- cbind(mySubDF, pop_type2 = mySubDF$pop_type)
mySubDF[mySubDF$pop_type2 == "del_araBAD", "pop_type2"] <- "wt"
mySubDF[mySubDF$pop_type2 == "del_mutL_araBAD", "pop_type2"] <- "del_mutL"

mySubDF <- cbind(mySubDF, pop_type3 = mySubDF$pop_type2)
mySubDF <- arrange(transform(mySubDF,pop_type=factor(pop_type3,levels=c("del_mutL","wt"))),pop_type3)
mySubDF$pop_type3 <- factor(mySubDF$pop_type3, labels = c("MMR-{}", "WT"))

mySubDF <- subset(mySubDF, (mut_type=="nonsynonymous" | mut_type=="indel" | mut_type=="sv"))
#mySubDF <- subset(mySubDF, pop_100==1 | pop_100==4 | pop_100==5)

mySelectList <- data.frame(matrix(nrow=0,ncol=ncol(mySubDF)))
colnames(mySelectList) <- colnames(mySubDF)

for (i in 1:nrow(GeneList)){
  mySelectList <- rbind(mySelectList, subset(mySubDF, gene_name==GeneList[i,1]))
}

mySelectList <- mySelectList[order(mySelectList$pop),]
write.table(mySelectList,'/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/select_gene_mutation_rate_related.txt',col.names=FALSE,row.names=FALSE,sep="\t", quote = FALSE)
write.table(subset(mySelectList, pop_100==1|pop_100==4 | pop_100==5),'/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/select_gene_mutation_rate_related_145.txt',col.names=TRUE,row.names=FALSE,sep="\t", quote = FALSE)
write.table(subset(mySelectList, pop_100==1|pop_100==2 | pop_100==3),'/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/select_gene_mutation_rate_related_123.txt',col.names=TRUE,row.names=FALSE,sep="\t", quote = FALSE)


myPosCount <- aggregate(pop~pos+gene_name+pop_type3+pop_100+mut_type, data=mySelectList[,c("pos","gene_name","pop","pop_type3","pop_100","mut_type")],FUN=length)
myGeneCount <- aggregate(pop~gene_name+pop_type3+pop_100+mut_type, data=mySelectList[,c("gene_name","pop","pop_type3","pop_100","mut_type")],FUN=length)
myPopCount <- aggregate(pop~pop_type3+pop_100, data=unique(mySelectList[,c("pop","pop_type3","pop_100")]),FUN=length)
myMutCount <- aggregate(pop~pop_type3+pop_100, data=mySelectList[,c("pop","pop_type3","pop_100")],FUN=length)