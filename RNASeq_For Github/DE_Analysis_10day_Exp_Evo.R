rm(list=ls())

library(readr)
library(DESeq2)
library(apeglm)
library(tximport)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)


setwd("~/10_day_ExpEvo")

###Call in Transcript Sample names and subset based on treatment
RNA_sample_id <- dir(file.path("~/10_day_ExpEvo/TranscriptData/"))
RNA_sample_id
Mid_RNA_sample_id<-RNA_sample_id[c(10,12,27,28)]
Mid_RNA_sample_id
Stat_RNA_sample_id<-RNA_sample_id[c(9,11,25,26)]
Stat_RNA_sample_id

###Call in Transcript Sample t_data files for both treatment types
Mid_RNA_Files <- file.path("~/10_day_ExpEvo/TranscriptData/", Mid_RNA_sample_id, "t_data.ctab")
Mid_RNA_Files
Stat_RNA_Files <- file.path("~/10_day_ExpEvo/TranscriptData/", Stat_RNA_sample_id, "t_data.ctab")
Stat_RNA_Files

###Translate ctab files for use by DeSeq for both treatment types
Mid_RNAtmp<- read_tsv(Mid_RNA_Files[1])
Mid_tx2gene<- Mid_RNAtmp[, c("t_name", "gene_name")]
Mid_txi <- tximport(Mid_RNA_Files, type = "stringtie", tx2gene = Mid_tx2gene, countsFromAbundance="lengthScaledTPM")
Mid_txi

saveRDS(Mid_txi, file = "~/10_day_ExpEvo/TranscriptAnalysis/Mid_txi.rds")

Stat_RNAtmp<- read_tsv(Stat_RNA_Files[1])
Stat_tx2gene<- Stat_RNAtmp[, c("t_name", "gene_name")]

Stat_txi <- tximport(Stat_RNA_Files, type = "stringtie", tx2gene = Stat_tx2gene,countsFromAbundance="lengthScaledTPM")
Stat_txi
saveRDS(Stat_txi, file = "~/10_day_ExpEvo/TranscriptAnalysis/Stat_txi.rds")


##Rename columns for both treatments based on the sample represented and associate with metadata information
RNA_SampleTable<-read.table("~/10_day_ExpEvo/TranscriptAnalysis/ClonePdata.test.txt", sep="\t", header=TRUE)
RNA_SampleTable
Mid_RNA_SampleTable<-RNA_SampleTable[which(RNA_SampleTable$treatment=='mid'),]
Stat_RNA_SampleTable<-RNA_SampleTable[which(RNA_SampleTable$treatment=='stat'),]

colnames(Mid_txi$counts)<-Mid_RNA_sample_id
rownames(Mid_RNA_SampleTable) <- colnames(Mid_txi$counts)

colnames(Stat_txi$counts)<-Stat_RNA_sample_id
rownames(Stat_RNA_SampleTable) <- colnames(Stat_txi$counts)

Mid_RNA_SampleTable
Stat_RNA_SampleTable



###Run DESeq Analysis
Mid_RNA_dds <- DESeqDataSetFromTximport(Mid_txi, Mid_RNA_SampleTable, ~condition)
Mid_RNA_dds<-DESeq(Mid_RNA_dds)
results(Mid_RNA_dds)
counts(Mid_RNA_dds, normalized=T)


Stat_RNA_dds <- DESeqDataSetFromTximport(Stat_txi, Stat_RNA_SampleTable, ~condition)
Stat_RNA_dds<-DESeq(Stat_RNA_dds)
resultsNames(Stat_RNA_dds)


###Get Significant Genes
Mid_CloneA_Res<-results(Mid_RNA_dds, name="condition_CloneA_vs_ancestor")
Mid_CloneA_ResLFC <- lfcShrink(Mid_RNA_dds, coef="condition_CloneA_vs_ancestor", type="apeglm")
summary(Mid_CloneA_Res)
plotMA(Mid_CloneA_ResLFC)
Mid_CloneA_Sig_RNA<-Mid_CloneA_Res[which(Mid_CloneA_Res$padj<0.05),]
Mid_CloneA_Sig_RNA<-as.data.frame(Mid_CloneA_Sig_RNA)
setDT(Mid_CloneA_Sig_RNA, keep.rownames = "Genes")[]

Mid_CloneB_Res<-results(Mid_RNA_dds, name="condition_CloneB_vs_ancestor")
Mid_CloneB_ResLFC <- lfcShrink(Mid_RNA_dds, coef="condition_CloneB_vs_ancestor", type="apeglm")
summary(Mid_CloneB_Res)
Mid_CloneB_Sig_RNA<-Mid_CloneB_Res[which(Mid_CloneBRes$padj<0.05),]
Mid_CloneB_Sig_RNA<-as.data.frame(Mid_CloneB_Sig_RNA)
setDT(Mid_CloneB_Sig_RNA, keep.rownames = "Genes")[]


Stat_CloneARes<-results(Stat_RNA_dds, name="condition_CloneA_vs_ancestor")
summary(Stat_CloneARes)
Stat_CloneA_Sig_RNA<-Stat_CloneARes[which(Stat_CloneARes$padj<0.05),]
Stat_CloneA_Sig_RNA<-as.data.frame(Stat_CloneA_Sig_RNA)
setDT(Stat_CloneA_Sig_RNA, keep.rownames = "Genes")[]

Stat_CloneBRes<-results(Stat_RNA_dds, name="condition_CloneB_vs_ancestor")
summary(Stat_CloneB_Res)
Stat_CloneB_Sig_RNA<-Stat_CloneBRes[which(Stat_CloneBRes$padj<0.05),]
Stat_CloneB_Sig_RNA<-as.data.frame(Stat_CloneB_Sig_RNA)
setDT(Stat_CloneB_Sig_RNA, keep.rownames = "Genes")[]
Stat_CloneB_Sig_RNA

###Get Overexpressed Genes
Mid_CloneA_Sig_RNA_Over<-Mid_CloneA_Sig_RNA[which(Mid_CloneA_Sig_RNA$log2FoldChange>0),]
Mid_CloneA_Sig_RNA_Over
Mid_CloneB_Sig_RNA_Over<-Mid_CloneB_Sig_RNA[which(Mid_CloneB_Sig_RNA$log2FoldChange>0),]
Mid_CloneB_Sig_RNA_Over
Stat_CloneA_Sig_RNA_Over<-Stat_CloneA_Sig_RNA[which(Stat_CloneA_Sig_RNA$log2FoldChange>0),]
Stat_CloneA_Sig_RNA_Over
Stat_CloneB_Sig_RNA_Over<-Stat_CloneB_Sig_RNA[which(Stat_CloneB_Sig_RNA$log2FoldChange>0),]
Stat_CloneB_Sig_RNA_Over

###Get Underexpressed Genes
Mid_CloneA_Sig_RNA_Under<-Mid_CloneA_Sig_RNA[which(Mid_CloneA_Sig_RNA$log2FoldChange<(0)),]
Mid_CloneA_Sig_RNA_Under
Mid_CloneB_Sig_RNA_Under<-Mid_CloneB_Sig_RNA[which(Mid_CloneB_Sig_RNA$log2FoldChange<(0)),]
Mid_CloneB_Sig_RNA_Under
Stat_CloneA_Sig_RNA_Under<-Stat_CloneA_Sig_RNA[which(Stat_CloneA_Sig_RNA$log2FoldChange<(0)),]
Stat_CloneA_Sig_RNA_Under
Stat_CloneB_Sig_RNA_Under<-Stat_CloneB_Sig_RNA[which(Stat_CloneB_Sig_RNA$log2FoldChange<(0)),]
Stat_CloneB_Sig_RNA_Under


#Write File of All DE Genes
Mid_CloneA_Sig_RNA_Over$Treatment<-"CloneA"
Mid_CloneA_Sig_RNA_Over$Condition<-"Mid-log"
Mid_CloneB_Sig_RNA_Over$Treatment<-"CloneB"
Mid_CloneB_Sig_RNA_Over$Condition<-"Mid-log"

Stat_CloneA_Sig_RNA_Over$Treatment<-"CloneA"
Stat_CloneA_Sig_RNA_Over$Condition<-"36h"
Stat_CloneB_Sig_RNA_Over$Treatment<-"CloneB"
Stat_CloneB_Sig_RNA_Over$Condition<-"36h"

Mid_CloneA_Sig_RNA_Under$Treatment<-"CloneA"
Mid_CloneA_Sig_RNA_Under$Condition<-"Mid-log"
Mid_CloneB_Sig_RNA_Under$Treatment<-"CloneB"
Mid_CloneB_Sig_RNA_Under$Condition<-"Mid-log"

Stat_CloneA_Sig_RNA_Under$Treatment<-"CloneA"
Stat_CloneA_Sig_RNA_Under$Condition<-"36h"
Stat_CloneB_Sig_RNA_Under$Treatment<-"CloneB"
Stat_CloneB_Sig_RNA_Under$Condition<-"36h"


All_DE_Trans<-rbind(Mid_CloneA_Sig_RNA_Over,Mid_CloneB_Sig_RNA_Over,Stat_CloneA_Sig_RNA_Over,Stat_CloneB_Sig_RNA_Over,Mid_CloneA_Sig_RNA_Under,Mid_CloneB_Sig_RNA_Under,Stat_CloneA_Sig_RNA_Under,Stat_CloneB_Sig_RNA_Under)
write.table(as.data.table(All_DE_Trans),paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Dataset_Clones403DE.txt'))





#Make Over and Under Expressed Genes Data Set
write.table(as.character(Mid_CloneA_Sig_RNA_Over$Genes),paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Mid_CloneA_Sig_RNA_Over.txt'))
write.table(as.character(Mid_CloneB_Sig_RNA_Over$Genes),paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Mid_CloneB_Sig_RNA_Over.txt'))

write.table(as.character(Stat_CloneA_Sig_RNA_Over$Genes),paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Stat_CloneA_Sig_RNA_Over.txt'))
write.table(as.character(Stat_CloneB_Sig_RNA_Over$Genes),paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Stat_CloneB_Sig_RNA_Over.txt'))


write.table(as.character(Mid_CloneA_Sig_RNA_Under$Genes),paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Mid_CloneA_Sig_RNA_Under.txt'))
write.table(as.character(Mid_CloneB_Sig_RNA_Under$Genes),paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Mid_CloneB_Sig_RNA_Under.txt'))

write.table(as.character(Stat_CloneA_Sig_RNA_Under$Genes),paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Stat_CloneA_Sig_RNA_Under.txt'))
write.table(as.character(Stat_CloneB_Sig_RNA_Under$Genes),paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Stat_CloneB_Sig_RNA_Under.txt'))


library(ggVennDiagram)

Mid_List<-list(A=Mid_CloneA_Sig_RNA$Genes, B=Mid_CloneB_Sig_RNA$Genes)
Stat_List<-list(A=Stat_CloneA_Sig_RNA$Genes, B=Stat_CloneB_Sig_RNA$Genes)



MidLog<-counts(Mid_RNA_dds, normalized=T)
MidLog<-as.data.frame(MidLog)
MidLog
MidLog$Genes<-rownames(MidLog)
MidLog$MeanExp <- (MidLog$`403_A_Mid`+MidLog$`403_B_Mid`)/2
MidLog$MeanAns <- (MidLog$`WT_Mid`+MidLog$`WT_Mid_BGI`)/2




MidLog_Clones_RMLow.melt<-melt(MidLog, id=c("Genes","MeanExp","MeanAns"))
MidLog_Clones_RMLow.melt
names(MidLog_Clones_RMLow.melt)<- c("Genes","MeanExp","MeanAns","Sample","Abundance")
MidCloneA_Trans <- MidLog_Clones_RMLow.melt$Sample %in% c("403_A_Mid")
MidCloneB_Trans <- MidLog_Clones_RMLow.melt$Sample %in% c("403_B_Mid")
MidAncestor_Trans <- MidLog_Clones_RMLow.melt$Sample %in% c("WT_Mid","WT_Mid_BGI")
MidLog_Clones_RMLow.melt$Treatment[MidCloneA_Trans] <- "Ecotype A"
MidLog_Clones_RMLow.melt$Treatment[MidCloneB_Trans] <- "Ecotype B"
MidLog_Clones_RMLow.melt$Treatment[MidAncestor_Trans] <- "Ancestor"

MidLog_Clones_RMLow<-MidLog[which(MidLog$MeanExp>=25),]
MidLog_Clones_RMLow$FoldChange<-MidLog_Clones_RMLow$`403_A_Mid`/MidLog_Clones_RMLow$`403_B_Mid`
MidLog_Clones_RMLow$FoldChange_Log2<-log(MidLog_Clones_RMLow$FoldChange,2)

Over403A<-MidLog_Clones_RMLow[which(MidLog_Clones_RMLow$FoldChange_Log2>=1.5),]
Under403A<-MidLog_Clones_RMLow[which(MidLog_Clones_RMLow$FoldChange_Log2<=-1.5),]
Over403A$IncExp<-"Ecotype A"
write.table(Over403A,paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Mid_403AvB_RNA_Over.txt'))
Under403A$IncExp<-"Ecotype B"
MidDEG_Ecotype<-rbind(Over403A,Under403A)

plotting_df<-MidDEG_Ecotype %>% group_by(IncExp) %>% summarise(Freq=n()) %>% mutate (Freq=if_else(IncExp=="Ecotype A",-Freq,Freq))
plotting_df
plotting_df$Type<-"Genes"
plotting_df$Label<- abs(plotting_df$Freq)

Mid403<-plotting_df %>% ggplot(aes(x=Type,y=Freq, group=IncExp, fill=IncExp))+
  geom_bar(stat="identity", width=0.75)+
  geom_hline(yintercept = 0, color="black")+
  geom_text(aes(label=Label,y=Freq*1.07),fontface="bold", size=4)+
  coord_flip()+
  ylab("DEGs Between Ecotypes at Midlog")+
  scale_y_continuous(breaks=seq(-100,100,25),limits=c(-100,100), labels=abs(seq(-100,100,25)))+
  scale_fill_manual(values=c("#FF0000","#0000FF"))+
  theme_minimal()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"), axis.title.y = element_blank(), axis.text.y = element_blank(),legend.title=element_blank(),legend.position = "top")



write.table(Under403A,paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Mid_403AvB_RNA_Under.txt'))



Mid_fab.Melt<-MidLog_Clones_RMLow.melt[MidLog_Clones_RMLow.melt$Genes=="fabD"|MidLog_Clones_RMLow.melt$Genes=="fabH"|MidLog_Clones_RMLow.melt$Genes=="micF"|MidLog_Clones_RMLow.melt$Genes=="ompC",]
Mid_fab.Melt
Mid_fab.Melt$Ratio<-Mid_fab.Melt$Abundance/Mid_fab.Melt$MeanAns
Mid_fab.Melt$Log2<-log(Mid_fab.Melt$Ratio,2)
Mid_fab.Melt.NoAns<-Mid_fab.Melt[which(Mid_fab.Melt$Treatment!="Ancestor"),]

MidFabPlot<-ggplot(data=Mid_fab.Melt.NoAns, aes(x=Treatment,y=Ratio, fill=Treatment))+
  stat_summary(fun.y = "mean", geom = "point",pch=22, size=3, color="black")+
  geom_hline(yintercept=1,linetype="dotted")+
  xlab("")+
  ylab("Abundance")+
  scale_fill_manual(values=c("#FF0000", "#0000FF"))+
  theme_bw()+theme(legend.position="none", text=element_text(size=12), axis.text =element_text(size=10),axis.text.x =element_text(angle=35,hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(~Genes, scales="free")




Stat<-counts(Stat_RNA_dds, normalized=T)
Stat<-as.data.frame(Stat)
Stat
Stat$Genes<-rownames(Stat)
Stat$MeanExp <- (Stat$`403_A_36h`+Stat$`403_B_36h`)/2
Stat$MeanAns<-(Stat$`WT_36h`+Stat$`WT_36h_BGI`)/2
Stat_Clones_RMLow<-Stat[which(Stat$MeanExp>25),]
Stat_Clones_RMLow.melt<-melt(Stat_Clones_RMLow, id=c("Genes","MeanExp","MeanAns"))
Stat_Clones_RMLow.melt
names(Stat_Clones_RMLow.melt)<- c("Genes","MeanExp","MeanAns","Sample","Abundance")
StatCloneA_Trans <- Stat_Clones_RMLow.melt$Sample %in% c("403_A_36h")
StatCloneB_Trans <- Stat_Clones_RMLow.melt$Sample %in% c("403_B_36h")
StatAncestor_Trans <- Stat_Clones_RMLow.melt$Sample %in% c("WT_36h","WT_36h_BGI")
Stat_Clones_RMLow.melt$Treatment[StatCloneA_Trans] <- "Ecotype A"
Stat_Clones_RMLow.melt$Treatment[StatCloneB_Trans] <- "Ecotype B"
Stat_Clones_RMLow.melt$Treatment[StatAncestor_Trans] <- "Ancestor"
Stat_Clones_RMLow$FoldChange<-Stat_Clones_RMLow$`403_A_36h`/Stat_Clones_RMLow$`403_B_36h`
Stat_Clones_RMLow$FoldChange_Log2<-log(Stat_Clones_RMLow$FoldChange,2)





Over403A_Stat<-Stat_Clones_RMLow[which(Stat_Clones_RMLow$FoldChange_Log2>=1.5),]
Over403A_Stat
Under403A_Stat<-Stat_Clones_RMLow[which(Stat_Clones_RMLow$FoldChange_Log2<=-1.5),]
Over403A_Stat$IncExp<-"Ecotype A"
write.table(Over403A_Stat,paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Stat_403AvB_RNA_Over.txt'))
Under403A_Stat$IncExp<-"Ecotype B"
write.table(Under403A_Stat,paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Stat_403AvB_RNA_Under.txt'))

StatDEG_Ecotype<-rbind(Over403A_Stat,Under403A_Stat)
plotting_df_Stat<-StatDEG_Ecotype %>% group_by(IncExp) %>% summarise(Freq=n()) %>% mutate (Freq=if_else(IncExp=="Ecotype A",-Freq,Freq))
plotting_df_Stat$Type<-"Genes"
plotting_df_Stat$Label<- abs(plotting_df_Stat$Freq)

Stat403<-plotting_df_Stat %>% ggplot(aes(x=Type,y=Freq, group=IncExp, fill=IncExp))+
  geom_bar(stat="identity", width=0.75)+
  geom_text(aes(label=Label,y=((Freq*1.092)+5)),fontface="bold", size=4)+
  geom_hline(yintercept = 0, color="black")+
  coord_flip()+
  ylab("DEGs Between Ecotypes at 36h")+
  scale_y_continuous(breaks=seq(-400,400,100),limits=c(-420,420), labels=abs(seq(-400,400,100)))+
  scale_fill_manual(values=c("#FF0000","#0000FF"))+
  theme_minimal()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"), axis.title.y = element_blank(), axis.text.y = element_blank(),legend.position="top",legend.title = element_blank())



Stat_DEOver403A
write.table(Stat_DEOver403A,paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Stat_403AvB_RNA_Over.txt'))
Stat_DEUnder403A
write.table(Stat_DEUnder403A,paste0('~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Stat_403AvB_RNA_Under.txt'))

StatPlot.Melt<-Stat_Clones_RMLow.melt[Stat_Clones_RMLow.melt$Genes=="entA"|Stat_Clones_RMLow.melt$Genes=="entS"|Stat_Clones_RMLow.melt$Genes=="ftnA"|Stat_Clones_RMLow.melt$Genes=="phoB"|Stat_Clones_RMLow.melt$Genes=="pstS",]
StatPlot.Melt$Ratio<-StatPlot.Melt$Abundance/StatPlot.Melt$MeanAns
StatPlot.Melt$Log2<-log(StatPlot.Melt$Ratio,2)
StatPlot.Melt.NoAns<-StatPlot.Melt[which(StatPlot.Melt$Treatment!="Ancestor"),]
StatPlot.Melt.NoAns
StatFepPlot<-ggplot(data=StatPlot.Melt.NoAns, aes(x=Treatment,y=Ratio, fill=Treatment))+
  stat_summary(fun.y = "mean", geom = "point",pch=22, size=3, color="black")+
  geom_hline(yintercept=1,linetype="dotted")+
  xlab("")+
  ylab("Abundance")+
  scale_fill_manual(values=c("#FF0000", "#0000FF"))+
  theme_bw()+theme(legend.position="none", text=element_text(size=12), axis.text =element_text(size=10),axis.text.x =element_text(angle=35,hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(~Genes, scales="free")


Model<-ggdraw()+
  draw_image("~/10_day_ExpEvo/TranscriptAnalysis/EcotypeAvB/Figure4_E.png")


GeneExpr<-plot_grid(Mid403,MidFabPlot,StatFepPlot,nrow=3, labels=c("A","C","D"),rel_heights = c(.6,1,1))


Fig4_E<-plot_grid(Stat403,Model,nrow=2, labels=c("B","E"), rel_heights = c(.6,2))

plot_grid(GeneExpr,Fig4_E, rel_widths = c(1,1))





############10-day clones vs Ancestor################
###Call in Transcript Sample names and subset based on treatment
RNA_sample_id <- dir(file.path("~/10_day_ExpEvo/TranscriptData/"))
RNA_sample_id
Mid_RNA_sample_id<-RNA_sample_id[c(10,12,27,28)]
Mid_RNA_sample_id
Stat_RNA_sample_id<-RNA_sample_id[c(9,11,25,26)]
Stat_RNA_sample_id


###Call in Transcript Sample t_data files for both treatment types
Mid_RNA_Files <- file.path("~/10_day_ExpEvo/TranscriptData/", Mid_RNA_sample_id, "t_data.ctab")
Mid_RNA_Files
Stat_RNA_Files <- file.path("~/10_day_ExpEvo/TranscriptData/", Stat_RNA_sample_id, "t_data.ctab")
Stat_RNA_Files

###Translate ctab files for use by DeSeq for both treatment types
Mid_RNAtmp<- read_tsv(Mid_RNA_Files[1])
Mid_tx2gene<- Mid_RNAtmp[, c("t_name", "gene_name")]
Mid_txi <- tximport(Mid_RNA_Files, type = "stringtie", tx2gene = Mid_tx2gene, countsFromAbundance="lengthScaledTPM")
Mid_txi

Stat_RNAtmp<- read_tsv(Stat_RNA_Files[1])
Stat_tx2gene<- Stat_RNAtmp[, c("t_name", "gene_name")]
Stat_txi <- tximport(Stat_RNA_Files, type = "stringtie", tx2gene = Stat_tx2gene, countsFromAbundance="lengthScaledTPM")
Stat_txi


##Rename columns for both treatments based on the sample represented and associate with metadata information
RNA_SampleTable<-read.table("~/10_day_ExpEvo/TranscriptAnalysis/Pdata.test.txt", sep="\t", header=TRUE)
RNA_SampleTable
Mid_RNA_SampleTable<-RNA_SampleTable[which(RNA_SampleTable$treatment=='mid'),]
Stat_RNA_SampleTable<-RNA_SampleTable[which(RNA_SampleTable$treatment=='stat'),]

colnames(Mid_txi$counts)<-Mid_RNA_sample_id
rownames(Mid_RNA_SampleTable) <- colnames(Mid_txi$counts)

colnames(Stat_txi$counts)<-Stat_RNA_sample_id
rownames(Stat_RNA_SampleTable) <- colnames(Stat_txi$counts)

Mid_RNA_SampleTable
Stat_RNA_SampleTable



###Run DESeq Analysis
Mid_RNA_dds <- DESeqDataSetFromTximport(Mid_txi, Mid_RNA_SampleTable, ~condition)
Mid_RNA_dds<-DESeq(Mid_RNA_dds)
resultsNames(Mid_RNA_dds)
Mid_RNA_dds


Stat_RNA_dds <- DESeqDataSetFromTximport(Stat_txi, Stat_RNA_SampleTable, ~condition)
Stat_RNA_dds<-DESeq(Stat_RNA_dds)
resultsNames(Stat_RNA_dds)

Mid_TenRes<-results(Mid_RNA_dds, name="condition_Ten_vs_ancestor")
Mid_Ten_Sig_RNA<-Mid_TenRes[which(Mid_TenRes$padj<0.05),]
Mid_Ten_Sig_RNA<-as.data.frame(Mid_Ten_Sig_RNA)
setDT(Mid_Ten_Sig_RNA, keep.rownames = "Genes")[]


Stat_TenRes<-results(Stat_RNA_dds, name="condition_Ten_vs_ancestor")
Stat_Ten_Sig_RNA<-Stat_TenRes[which(Stat_TenRes$padj<0.05),]
Stat_Ten_Sig_RNA<-as.data.frame(Stat_Ten_Sig_RNA)
setDT(Stat_Ten_Sig_RNA, keep.rownames = "Genes")[]

Mid_Ten_Sig_RNA_Over<-Mid_Ten_Sig_RNA[which(Mid_Ten_Sig_RNA$log2FoldChange>0),]
Stat_Ten_Sig_RNA_Over<-Stat_Ten_Sig_RNA[which(Stat_Ten_Sig_RNA$log2FoldChange>0),]

Mid_Ten_Sig_RNA_Under<-Mid_Ten_Sig_RNA[which(Mid_Ten_Sig_RNA$log2FoldChange<0),]
Stat_Ten_Sig_RNA_Under<-Stat_Ten_Sig_RNA[which(Stat_Ten_Sig_RNA$log2FoldChange<0),]
Mid_Ten_Sig_RNA_Over$Condition<-"Mid-log"
Mid_Ten_Sig_RNA_Under$Condition<-"Mid-log"
Stat_Ten_Sig_RNA_Over$Condition<-"36h"
Stat_Ten_Sig_RNA_Under$Condition<-"36h"


All_DE_Trans<-rbind(Mid_Ten_Sig_RNA_Over,Stat_Ten_Sig_RNA_Over,Mid_Ten_Sig_RNA_Under,Stat_Ten_Sig_RNA_Under)
write.table(as.data.table(All_DE_Trans),paste0('~/10_day_ExpEvo/TranscriptAnalysis/403vWT/Dataset_403vAns.txt'))
Stat_Ten_Sig_RNA_Over$Genes
Stat_Ten_Sig_RNA_Under[which(Stat_Ten_Sig_RNA_Under$Genes=="gadX"),]

