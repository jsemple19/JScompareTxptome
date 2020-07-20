library(readxl)
library(ggVennDiagram)
library(ggplot2)
source("functions.R")

outPath="."

salmon<-readRDS(paste0(outPath,"/rds/salmon_DESeq2_fullResults.rds"))
star<-readRDS(paste0(outPath,"/rds/star_DESeq2_fullResults.rds"))
quasr<-readRDS(paste0(outPath,"/rds/quasr_DESeq2_fullResults.rds"))

padjVal=0.05
lfcVal=0.5


#####################################################
## compare to public data
#####################################################

kramerURL<-"https://doi.org/10.1371/journal.pgen.1005698.s011"
kramerFileName="Kramer_2015_PlotGen_S3file.xlsx"
#download.file(url=kramerURL,destfile=kramerFileName)
library("readxl")

kramer<-read_excel(kramerFileName,col_types=c(rep("text",3),rep("numeric",30)))
names(kramer)
kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]

lfcVal<-0.5
padjVal<-0.05
idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


salmonSig<-salmon$wormbase[!is.na(salmon$padj) &
                             salmon$padj<padjVal &
                             abs(salmon$log2FoldChange)>lfcVal]


x<-list(salmon=salmonSig, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p1<-ggVennDiagram(x) + ggtitle(label=paste0("salmon vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))


starSig<-star$wormbase[!is.na(star$padj) &
                         star$padj<padjVal &
                         star$log2FoldChange< -lfcVal &
                         star$chr!="chrX"]

x<-list(star=starSig, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p2<-ggVennDiagram(x) + ggtitle(label=paste0("star vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))


quasrSig<-quasr$wormbase[!is.na(quasr$padj) &
                           quasr$padj<padjVal &
                           quasr$log2FoldChange< -lfcVal &
                           quasr$chr!="chrX"]

x<-list(quasr=quasrSig, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p3<-ggVennDiagram(x) + ggtitle(label=paste0("quasr vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))

p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_Kramer2015_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", formatC(lfcVal,format="e",digits=0),
                                ".pdf"),
                plot=p, device="pdf",width=29,height=10,units="cm")




###############################
## different cutoffs
###############################

padjVal<-0.05
lfcVal<-0

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


salmonSig<-salmon$wormbase[!is.na(salmon$padj) &
                             salmon$padj<padjVal &
                             abs(salmon$log2FoldChange)>lfcVal]


x<-list(salmon=salmonSig, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p1<-ggVennDiagram(x) + ggtitle(label=paste0("salmon vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))

lfcVal<-0.5

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


salmonSig<-salmon$wormbase[!is.na(salmon$padj) &
                             salmon$padj<padjVal &
                             abs(salmon$log2FoldChange)>lfcVal]


x<-list(salmon=salmonSig, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p2<-ggVennDiagram(x) + ggtitle(label=paste0("salmon vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))

lfcVal<-0.75

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


salmonSig<-salmon$wormbase[!is.na(salmon$padj) &
                             salmon$padj<padjVal &
                             abs(salmon$log2FoldChange)>lfcVal]


x<-list(salmon=salmonSig, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p3<-ggVennDiagram(x) + ggtitle(label=paste0("salmon vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))


lfcVal<-1

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


salmonSig<-salmon$wormbase[!is.na(salmon$padj) &
                             salmon$padj<padjVal &
                             abs(salmon$log2FoldChange)>lfcVal]


x<-list(salmon=salmonSig, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p4<-ggVennDiagram(x) + ggtitle(label=paste0("salmon vs Kramer(2015): |lfc|>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_bestCutoff_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")



###############################
## chrX upregulated
###############################

dim(kramer)
dim(salmon)
idx<-match(kramer$Gene_WB_ID, salmon$wormbase)
kramer$chr<-salmon$chr[idx]
kramer<-kramer[!is.na(kramer$chr),]


padjVal=0.01

lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
      kramer$dpy27_RNAi_L3_padj < padjVal &
      kramer$dpy27_RNAi_L3_log2_fold_change > lfcVal &
      kramer$chr=="chrX"
kramerdpy27dcc<-kramer[idx,]


idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
  kramer$dpy21_mutant_L3_padj < padjVal &
  kramer$dpy21_mutant_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"
kramerdpy21dcc<-kramer[idx,]

x<-list(salmon=salmondcc$wormbase, dpy27=kramerdpy27dcc$Gene_WB_ID,
        dpy21=kramerdpy21dcc$Gene_WB_ID)

p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))



lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
  kramer$dpy27_RNAi_L3_padj < padjVal &
  kramer$dpy27_RNAi_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"
kramerdpy27dcc<-kramer[idx,]


idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
  kramer$dpy21_mutant_L3_padj < padjVal &
  kramer$dpy21_mutant_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"

kramerdpy21dcc<-kramer[idx,]

x<-list(salmon=salmondcc$wormbase, dpy27=kramerdpy27dcc$Gene_WB_ID,
        dpy21=kramerdpy21dcc$Gene_WB_ID)

p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))



lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
  kramer$dpy27_RNAi_L3_padj < padjVal &
  kramer$dpy27_RNAi_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"
kramerdpy27dcc<-kramer[idx,]


idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
  kramer$dpy21_mutant_L3_padj < padjVal &
  kramer$dpy21_mutant_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"

kramerdpy21dcc<-kramer[idx,]

x<-list(salmon=salmondcc$wormbase, dpy27=kramerdpy27dcc$Gene_WB_ID,
        dpy21=kramerdpy21dcc$Gene_WB_ID)

p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))



lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")

idx<-!is.na(kramer$dpy27_RNAi_L3_padj) &
  kramer$dpy27_RNAi_L3_padj < padjVal &
  kramer$dpy27_RNAi_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"
kramerdpy27dcc<-kramer[idx,]


idx<-!is.na(kramer$dpy21_mutant_L3_padj) &
  kramer$dpy21_mutant_L3_padj < padjVal &
  kramer$dpy21_mutant_L3_log2_fold_change > lfcVal &
  kramer$chr=="chrX"

kramerdpy21dcc<-kramer[idx,]

x<-list(salmon=salmondcc$wormbase, dpy27=kramerdpy27dcc$Gene_WB_ID,
        dpy21=kramerdpy21dcc$Gene_WB_ID)

p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_XchrUp_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")



###############################
## classical DCC genes
###############################

pubDCC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/published_DCCgr.rds")
pubNDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/published_NDCgr.rds")


padjVal=0.01
lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, DC=pubDCC$wormbase,
        nonDC=pubNDC$wormbase)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, DC=pubDCC$wormbase,
        nonDC=pubNDC$wormbase)
p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, DC=pubDCC$wormbase,
        nonDC=pubNDC$wormbase)
p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))

lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, DC=pubDCC$wormbase,
        nonDC=pubNDC$wormbase)
p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_papers_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")


padjVal=0.05
lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, DC=pubDCC$wormbase,
        nonDC=pubNDC$wormbase)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, DC=pubDCC$wormbase,
        nonDC=pubNDC$wormbase)
p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, DC=pubDCC$wormbase,
        nonDC=pubNDC$wormbase)
p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))

lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, DC=pubDCC$wormbase,
        nonDC=pubNDC$wormbase)
p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_papers_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")




###############################
## Jans 2009
###############################

JansDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/Jans2009_DCgr.rds")
JansNDC<-readRDS("/Users/semple/Documents/MeisterLab/dSMF/DCgenes/Jans2009_NDCgr.rds")


padjVal=0.01
lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, JansDC=JansDC$wormbase,
        JansNDC=JansNDC$wormbase)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, JansDC=JansDC$wormbase,
        JansNDC=JansNDC$wormbase)
p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, JansDC=JansDC$wormbase,
        JansNDC=JansNDC$wormbase)
p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))

lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, JansDC=JansDC$wormbase,
        JansNDC=JansNDC$wormbase)
p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_Jans2009_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")


padjVal=0.05
lfcVal=0
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, JansDC=JansDC$wormbase,
        JansNDC=JansNDC$wormbase)
p1<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.5
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, JansDC=JansDC$wormbase,
        JansNDC=JansNDC$wormbase)
p2<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


lfcVal=0.75
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, JansDC=JansDC$wormbase,
        JansNDC=JansNDC$wormbase)
p3<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))

lfcVal=1
salmondcc<-filterResults(salmon,padjVal,lfcVal,direction="gt",chr="chrX")
x<-list(salmon=salmondcc$wormbase, JansDC=JansDC$wormbase,
        JansNDC=JansNDC$wormbase)
p4<-ggVennDiagram(x) + ggplot2::ggtitle(label=paste0("chrX DCC genes: lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_Jans2009_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-1",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")


###############################
## Jans 2009 vs Kramer
###############################

kramer<-read_excel(kramerFileName,col_types=c(rep("text",3),rep("numeric",30)))
names(kramer)
kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]

padjVal<-0.05
lfcVal<-0
idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


x<-list(JansDC=JansDC$wormbase, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p1<-ggVennDiagram(x) + ggtitle(label=paste0("Jans DC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal))


x<-list(JansNDC=JansNDC$wormbase, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p2<-ggVennDiagram(x) + ggtitle(label=paste0("Jans NDC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal))


lfcVal<-0.5
idx<-!is.na(kramer$dpy27_RNAi_L3_padj) & abs(kramer$dpy27_RNAi_L3_log2_fold_change)>lfcVal & kramer$dpy27_RNAi_L3_padj<padjVal
kramerDpy27<-kramer[idx,]

idx<-!is.na(kramer$dpy21_mutant_L3_padj) & abs(kramer$dpy21_mutant_L3_log2_fold_change)>lfcVal & kramer$dpy21_mutant_L3_padj<padjVal
kramerDpy21<-kramer[idx,]


x<-list(JansDC=JansDC$wormbase, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p3<-ggVennDiagram(x) + ggtitle(label=paste0("Jans DC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal))



x<-list(JansNDC=JansNDC$wormbase, dpy27=kramerDpy27$Gene_WB_ID,
        dpy21=kramerDpy21$Gene_WB_ID)

p4<-ggVennDiagram(x) + ggtitle(label=paste0("Jans NDC vs Kramer(2015): lfc>", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
ggplot2::ggsave(filename=paste0(outPath, "/plots/venn_Jans2009vKramer2015_padj",
                                formatC(padjVal,format="e",digits=0),
                                "_lfc", "0-0.5",".pdf"),
                plot=p, device="pdf",width=19,height=20,units="cm")

