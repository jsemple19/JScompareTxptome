outPath="."

salmon<-readRDS(paste0(outPath,"/rds/salmon_DESeq2_fullResults.rds"))
star<-readRDS(paste0(outPath,"/rds/star_DESeq2_fullResults.rds"))
quasr<-readRDS(paste0(outPath,"/rds/quasr_DESeq2_fullResults.rds"))

padjVal=0.05
lfcVal=0.5



##############
## compare salmon to star
##############


ddff<-merge(as.data.frame(salmon),as.data.frame(star),by="wormbase")


pdf(file=paste0(outPath, "/plots/compare_salmonVstar_cor_padj",
                formatC(padjVal,format="e",digits=0), "_lfc",
                formatC(lfcVal,format="e",digits=0),".pdf"),
    width=11,height=5.5,paper="a4r")
par(mfrow=c(1,2))
plot(ddff$log2FoldChange.x, ddff$log2FoldChange.y,
     xlab="salmon log2FC",
     ylab="star log2FC",pch=16, col="#BB110077",
     main=paste0("statistically significant in salmon (lfc>",
                 lfcVal,", p<",padjVal,")"))
abline(v=c(-lfcVal,lfcVal),h=c(-lfcVal,lfcVal),lty=2,col="light grey")
idx<-which(!is.na(ddff$padj.x) & abs(ddff$log2FoldChange.x)>lfcVal & ddff$padj.x<padjVal)
points(ddff$log2FoldChange.x[idx], ddff$log2FoldChange.y[idx], pch=3)



plot(ddff$log2FoldChange.x, ddff$log2FoldChange.y,
     xlab="salmon log2FC",
     ylab="star log2FC",pch=16, col="#BB110077",
     main=paste0("statistically significant in star (lfc>",
                 lfcVal, ", p<",padjVal,")"))
#abline(v=c(-1,1),h=c(-1,1),lty=2,col="dark grey")
abline(v=c(-lfcVal,lfcVal),h=c(-lfcVal,lfcVal),lty=2,col="light grey")
idx<-which(!is.na(ddff$padj.y) & abs(ddff$log2FoldChange.y)>lfcVal & ddff$padj.x<padjVal)
points(ddff$log2FoldChange.x[idx], ddff$log2FoldChange.y[idx], pch=2)

dev.off()

# look manually at which are different in salmon but not star
idx<-which(abs(ddff$log2FoldChange.x)>lfcVal & abs(ddff$log2FoldChange.y)<lfcVal)
ddff[idx,]


##############
## compare salmon to quasr
##############
df1Name="salmon"
df2Name="quasr"
ddff<-merge(as.data.frame(salmon),as.data.frame(quasr),by="wormbase")


pdf(file=paste0(outPath, "/plots/compare_",df1Name,"V",df2Name,"_cor_padj",
                formatC(padjVal,format="e",digits=0), "_lfc",
                formatC(lfcVal,format="e",digits=0),".pdf"),
    width=11, height=5.5, paper="a4r")
par(mfrow=c(1,2))
plot(ddff$log2FoldChange.x, ddff$log2FoldChange.y,
     xlab=paste0(df1Name," log2FC"),
     ylab=paste0(df2Name," log2FC"),pch=16, col="#BB110077",
     main=paste0("statistically significant in ",df1Name,
                 " (lfc>",lfcVal,", p<",padjVal,")"))
abline(v=c(-lfcVal,lfcVal),h=c(-lfcVal,lfcVal),lty=2,col="light grey")
idx<-which(!is.na(ddff$padj.x) & abs(ddff$log2FoldChange.x)>lfcVal & ddff$padj.x<padjVal)
points(ddff$log2FoldChange.x[idx], ddff$log2FoldChange.y[idx], pch=3)



plot(ddff$log2FoldChange.x, ddff$log2FoldChange.y,
     xlab=paste0(df1Name," log2FC"),
     ylab=paste0(df2Name," log2FC"),pch=16, col="#BB110077",
     main=paste0("statistically significant in ",df2Name,
                 " (lfc>",lfcVal,", p<",padjVal,")"))
#abline(v=c(-1,1),h=c(-1,1),lty=2,col="dark grey")
abline(v=c(-lfcVal,lfcVal),h=c(-lfcVal,lfcVal),lty=2,col="light grey")
idx<-which(!is.na(ddff$padj.y) & abs(ddff$log2FoldChange.y)>lfcVal & ddff$padj.x<padjVal)
points(ddff$log2FoldChange.x[idx], ddff$log2FoldChange.y[idx], pch=2)

dev.off()

# look manually at which are different in salmon but not star
idx<-which(abs(ddff$log2FoldChange.x)>lfcVal & abs(ddff$log2FoldChange.y)<lfcVal)
ddff[idx,]


##############
#############
library(ggVennDiagram)

par(mfrow=c(3,1))

salmonSig<-salmon$wormbase[!is.na(salmon$padj) &
                             salmon$padj<padjVal &
                             abs(salmon$log2FoldChange)>lfcVal]
starSig<-star$wormbase[!is.na(star$padj) &
                             star$padj<padjVal &
                             abs(star$log2FoldChange)>lfcVal]
quasrSig<-quasr$wormbase[!is.na(quasr$padj) &
                             quasr$padj<padjVal &
                             abs(quasr$log2FoldChange)>lfcVal]


x<-list(salmon=salmonSig, star=starSig, quasr=quasrSig)

p1<-ggVennDiagram(x) + ggtitle(label=paste0("All significant genes: |lfc|>", lfcVal, ", padj<",padjVal))


## X chr genes

salmonSig<-salmon$wormbase[!is.na(salmon$padj) &
                             salmon$padj<padjVal &
                             salmon$log2FoldChange>lfcVal &
                             salmon$chr=="chrX"]
starSig<-star$wormbase[!is.na(star$padj) &
                         star$padj<padjVal &
                         star$log2FoldChange>lfcVal &
                         star$chr=="chrX"]
quasrSig<-quasr$wormbase[!is.na(quasr$padj) &
                           quasr$padj<padjVal &
                           quasr$log2FoldChange>lfcVal &
                           quasr$chr=="chrX"]


x<-list(salmon=salmonSig, star=starSig, quasr=quasrSig)

p2<-ggVennDiagram(x) + ggtitle(label=paste0("Up on chrX: lfc>", lfcVal, ", padj<",padjVal))




## autosomal genes

salmonSig<-salmon$wormbase[!is.na(salmon$padj) &
                             salmon$padj<padjVal &
                             salmon$log2FoldChange< -lfcVal &
                             salmon$chr!="chrX"]
starSig<-star$wormbase[!is.na(star$padj) &
                         star$padj<padjVal &
                         star$log2FoldChange< -lfcVal &
                         star$chr!="chrX"]
quasrSig<-quasr$wormbase[!is.na(quasr$padj) &
                           quasr$padj<padjVal &
                           quasr$log2FoldChange< -lfcVal &
                           quasr$chr!="chrX"]


x<-list(salmon=salmonSig, star=starSig, quasr=quasrSig)

p3<-ggVennDiagram(x) + ggtitle(label=paste0("Down on autosomes: lfc< -", lfcVal, ", padj<",padjVal))


p<-ggpubr::ggarrange(p1,p2,p3,ncol=3,nrow=1)
ggplot2::ggsave(filename=paste0(outPath, "/plots/", "venn_", "padj",
                                formatC(padjVal,format="e",digits=0),
                       "_lfc", formatC(lfcVal,format="e",digits=0), ".pdf"),
       plot=p, device="pdf",width=29,height=10,units="cm")
