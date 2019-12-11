rm(list = ls())

library(plyr)  #round_any
library(dplyr)  #full_join
library(ggplot2)
library(RVAideMemoire)
cbPalette <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#CC79A7", "#999999")

basalDF <- read.table("/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/list_basal_fixed_alleles.txt", header=T)
major_minor_DF  <- read.table("/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/list_major_minor_clade_fixed_alleles.txt", header=T)

basalDF <- cbind(basalDF, pop_100=floor(basalDF$pop/100))
major_minor_DF <- cbind(major_minor_DF, pop_100=floor(major_minor_DF$pop/100))

basalDF <- cbind(basalDF, day=round_any(basalDF$fixedTP/((basalDF$pop_100==1)*3.3+(basalDF$pop_100==2)*3.3*4+(basalDF$pop_100==3)*3.3*7+(basalDF$pop_100==4)*0.33+(basalDF$pop_100==5)*0.22),10))
major_minor_DF <- cbind(major_minor_DF, day=round_any(major_minor_DF$fixedTP/((major_minor_DF$pop_100==1)*3.3+(major_minor_DF$pop_100==2)*3.3*4+(major_minor_DF$pop_100==3)*3.3*7+(major_minor_DF$pop_100==4)*0.33+(major_minor_DF$pop_100==5)*0.22),10))

basalCount <- aggregate(pos~pop+pop_type+day, data=basalDF[,c("pop","pop_type","day","pos")], FUN=length)
major_minor_Count <- aggregate(pos~pop+pop_type+day, data=major_minor_DF[,c("pop","pop_type","day","pos")], FUN=length)

colnames(basalCount)[which(colnames(basalCount)=="pos")] <- "num_fixed_basal"
colnames(major_minor_Count)[which(colnames(major_minor_Count)=="pos")] <- "num_fixed_major_minor"

allCount <- full_join(basalCount,major_minor_Count, by=c("pop","pop_type","day"))
allCount <- allCount[order(allCount$pop, allCount$day),]

allCount$num_fixed_basal[which(is.na(allCount$num_fixed_basal))]<-0
allCount$num_fixed_major_minor[which(is.na(allCount$num_fixed_major_minor))]<-0

allCount <- cbind(allCount, ratio=allCount$num_fixed_basal/(allCount$num_fixed_basal+allCount$num_fixed_major_minor))

write.table(allCount,"/Users/wei-chinho/Documents/MURI2/201_hmm/hmm_all_summary/calculate_num_fixed_basal_byTP.out")