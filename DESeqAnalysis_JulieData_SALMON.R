library(DESeq2)
library(Organism.dplyr)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library("TxDb.Celegans.UCSC.ce11.refGene")
library("TxDb.Celegans.UCSC.ce11.ensGene")
library(tximport)
library(GenomicFeatures)
library(ggplot2)
library("RColorBrewer")
library("PoiClaClu")
library("pheatmap")
library(tidyr)
library(EnhancedVolcano)
library(affy)
library("gplots")
#library(REBayes)
# get funciton for converting gene names from WBID to publicID
source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")

###############################################################
### some variables
###############################################################
fileNamePrefix="salmon_"
outPath="."
genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

dir.create(paste0(outPath,"/rds/"),recursive=T)
dir.create(paste0(outPath,"/plots/"),recursive=T)
dir.create(paste0(outPath,"/txt/"),recursive=T)
dir.create(paste0(outPath,"/tracks/"),recursive=T)

#bamList<-list.files(path="bamSTAR",pattern="Aligned.out.bam$",
#                       full.names=TRUE)

#sampleNames<-gsub("bam/","",bamList)
#sampleNames<-gsub("Aligned.out.bam","",sampleNames)

fileList<-read.table(paste0(outPath,"/fastqList.txt"),stringsAsFactors=F,header=T)


sampleNames<-paste(fileList$sampleName,fileList$repeatNum,sep="_")

fileNames<-paste0(outPath,"/salmon/mRNA/",sampleNames,"/quant.sf")

sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,stringsAsFactors=F)

# extract the technical replicate variable
#sampleTable$replicate<-as.factor(gsub("^[0-9]{3}_", "",sampleTable$sampleName))
sampleTable$replicate=fileList$repeatNum

# extract the strain variable
#sampleTable$strain<-factor(gsub("_HS[1-3]$", "",sampleTable$sampleName),levels=c("500","493"))
sampleTable$strain<-factor(as.character(fileList$sampleName),levels=c("500","493"))
sampleTable$dpy26<-sampleTable$strain
levels(sampleTable$dpy26)[levels(sampleTable$dpy26)=="493"]<-"TEVcs"
levels(sampleTable$dpy26)[levels(sampleTable$dpy26)=="500"]<-"wt"

###############################################################
### create metadata
###############################################################

if(!file.exists(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.",
                       genomeVer, ".annotations.sqlite"))){
  dir.create(paste0(genomeDir,"/annotations"),recursive=T)
  system(paste0("curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.",
                genomeVer,".annotations.gff3.gz -o ",genomeDir,
                "/annotations/c_elegans.PRJNA13758.",genomeVer,
                ".annotations.gff3.gz"))

  system(paste0("gunzip ",genomeDir,"/annotations/c_elegans.PRJNA13758.",
                genomeVer,".annotations.gff3.gz"))
  si<-seqinfo(Celegans)
  genome(si)<-genomeVer
  seqnames(si)<-gsub("M","MtDNA",gsub("chr","",seqnames(si)))
  wstxdb<-makeTxDbFromGFF(file=paste0(genomeDir,
                                      "/annotations/c_elegans.PRJNA13758.",
                                      genomeVer,".annotations.gff3"),
                          format="gff3",organism="Caenorhabditis elegans",
                          chrominfo=si)

  saveDb(wstxdb,paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
                     ".annotations.sqlite"))
  file.remove(paste0(genomeDir,"/annotations/c_elegans.PRJNA13758.",
                     genomeVer, ".annotations.gff3"))
}

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

srcref <- src_organism("TxDb.Celegans.UCSC.ce11.refGene")
srcens <- src_organism("TxDb.Celegans.UCSC.ce11.ensGene")

#src_tbls(srcref)
#tbl(srcref, "id")
#columns(srcref)
#supportedFilters(srcref)
#wbgenes<-WormbaseFilter("WBGene","startsWith")

metadata<-inner_join(tbl(srcref, "id"), tbl(srcref, "ranges_gene")) %>%
  dplyr::select(wormbase,alias,genename,gene_chrom,gene_start, gene_end, gene_strand) %>%
  collect %>% GenomicRanges::GRanges()



###############################################################
### get samples into DESeq2
###############################################################

# import the count matrices
txi<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene)

# read samples into DESeq2
dds <- DESeqDataSetFromTximport(txi=txi, colData=sampleTable,
                                design=~replicate+dpy26)



###############################################################
### DESeq2 differential expression analysis (using negative binomial distribution)
###############################################################

#dds<-collapseReplicates(dds,groupby=samples$sampleID,renameCols=T)
#dds <- DESeq(dds)
# This function performs a default analysis through the steps:
#   1. estimation of size factors: estimateSizeFactors
#   2. estimation of dispersion: estimateDispersions
#   3. Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
#   returns a DESeqDataSet object

idx<-match(rownames(dds),metadata$wormbase)
# add gene and chormosome names as metadata
featureData <- data.frame(gene=rownames(dds),
                          chr=as.vector(seqnames(metadata))[idx]) #,

rowData(dds) <- DataFrame(mcols(dds), featureData)

#only take expressed genes
#dds <- dds[ rowSums(counts(dds)) > 1, ]
#print("Number of expressed genes:")
#dim(dds)[1]
#16606 genes from 20127

dds<-DESeq(dds)




###############################################################
### get significant genes
###############################################################

Threshold=0.05

res<-results(dds)
sink(file=paste0(outPath,"/txt/",fileNamePrefix,"logfile.txt"), append=TRUE,
     type="output")
cat("Number of genes that change expression at different padj cutoffs:\n")
print(summary(res))
print(summary(res,alpha=0.05))
print(summary(res,alpha=0.01))
sink()

### add metadata
res$wormbase<-rownames(res)
idx<-match(rownames(res),metadata$wormbase)
res$chr<-factor(seqnames(metadata),levels=paste0("chr",c("I","II","III","IV","V","X")))[idx]
res$start<-as.vector(start(metadata))[idx]
res$end<-as.vector(end(metadata))[idx]
res$strand<-as.vector(strand(metadata))[idx]

# shrink LFC estimates
#resultsNames(dds) # to get names of coefficients
resLFC<-lfcShrink(dds,coef="dpy26_TEVcs_vs_wt",type="apeglm",res=res)
class(resLFC)
### add metadata
resLFC$wormbase<-rownames(resLFC)
idx<-match(rownames(resLFC),metadata$wormbase)
resLFC$chr<-factor(seqnames(metadata),levels=paste0("chr",c("I","II","III","IV","V","X")))[idx]
resLFC$start<-as.vector(start(metadata))[idx]
resLFC$end<-as.vector(end(metadata))[idx]
resLFC$strand<-as.vector(strand(metadata))[idx]
saveRDS(resLFC,file=paste0(outPath,"/rds/", fileNamePrefix,
                           "DESeq2_fullResults.rds"))

#export csv with ordered results
write.csv(resLFC[order(resLFC$padj),],
          file=paste0(outPath,"/txt/", fileNamePrefix,"DESeq2_resultsTable.csv"),
          quote=F,row.names=F)

# remove NAs
res<-na.omit(res)
resLFC<-na.omit(resLFC)


#################
#### Basic sample stats
#################

## basic sample stats
sink(file=paste0(outPath,"/txt/", fileNamePrefix,"logfile.txt"),
     append=FALSE, type="output")
statsPerSample<-data.frame(t(apply(counts(dds),2,summary)))
rownames(statsPerSample)<-colData(dds)$sampleName
colnames(statsPerSample)<-c("min", "Q1", "median", "mean", "Q3", "max")
statsPerSample$zeros <- apply(counts(dds)==0, 2, sum)
statsPerSample$percZeros <- round(100*statsPerSample$zeros/nrow(counts(dds)),1)
print(statsPerSample)
sink()

pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"sampleQC.pdf"), width=8,height=8,paper="a4")

#########
## sample counts summary: boxplots and density plots
########
## box plots
epsilon <- 1 # pseudo-count to avoid problems with log(0)
boxplot(log2(counts(dds) + epsilon), col=c("blue","red")[as.factor(colData(dds)$strain)], pch=".",
        horizontal=TRUE, cex.axis=1,
        las=1, ylab=NULL, names=colData(dds)$sampleName, xlab="log2(counts +1)")


## density plots
plotDensity(log2(counts(dds) + epsilon), lty=1, col=as.factor(colData(dds)$sampleName), lwd=2)
grid()
legend("topright", legend=colData(dds)$sampleName, col=as.factor(colData(dds)$sampleName), lwd=2)

##########
## pairwise correlation between genes in different samples
##########
#Define a function to draw a scatter plot for a pair of variables (samples) with density colors
plotFun <- function(x,y){
  dns <- densCols(x,y);
  points(x,y, col=dns, pch=".", panel.first=grid());
  #  abline(a=0, b=1, col="brown")
}

corFun <- function(x,y){
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 1.5*cor(x, y))
}

pairs(log2(counts(dds) + epsilon),
      panel=plotFun, lower.panel=corFun, labels=colData(dds)$sampleName)

#########
## plot dispersion estimates
#########
plotDispEsts(dds)

#######
## plot results filtering threshold
#######
plot(metadata(resLFC)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter",main="Threshold for independant filtering of results")
lines(metadata(resLFC)$lo.fit, col="red")
abline(v=metadata(resLFC)$filterTheta)
legend("topright",legend=paste0("Read count \nthreshold: ",
                                round(metadata(resLFC)$filterThreshold,2)))

#########
### sample to sample heatmap
#########
vsd <- vst(dds, blind=TRUE)
colnames(vsd)<-colData(dds)$sampleName
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colData(dds)$sampleName
colnames(sampleDistMatrix) <- colData(dds)$sampleName
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#########
## Heatmap of most highly expressed genes
#########
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:200]
df <- as.data.frame(colData(dds)[,c("strain","replicate")])
rownames(df)<-colData(dds)$sampleName
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df, main="Top 200 expressed genes")

##########
## heirarchical clustering of most significantly changed genes
##########
# select gene names based on FDR (5%)
gene.kept <- rownames(resLFC)[resLFC$padj <= Threshold & !is.na(resLFC$padj)]

# Retrieve the normalized counts for gene of interest
countTable.kept <- log2(counts(dds) + epsilon)[gene.kept, ]
dim(countTable.kept)
colnames(countTable.kept)<-colData(dds)$sampleName

# Perform the hierarchical clustering with
# A distance based on Pearson-correlation coefficient
# and average linkage clustering as agglomeration criteria
heatmap.2(as.matrix(countTable.kept),
          scale="row",
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          trace="none",
          density="none",
          labRow="",
          #labCol = names(countTable.kept),
          cexCol=1,
          main=paste0("Significantly changed genes (p<",Threshold,")"))


###########
## pca
###########
plotPCA(vsd, intgroup=c("strain", "replicate"))

dev.off()


##########
## plot individual genes
##########

pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"topGenes_normCounts.pdf"), width=8,height=11,paper="a4")
par(mfrow=c(5,2))
selectedGenes <- rownames(resLFC)[order(resLFC$padj)][1:20]
for (g in selectedGenes) {
  barplot(counts(dds, normalized=TRUE)[g,],
          col=colData(dds)$dpy26,
          main=g, las=2, cex.names=1,names.arg=colData(dds)$sampleName)
  legend("topleft",legend=levels(colData(dds)$dpy26),fill=c(1,2), cex=0.7)
}
dev.off()


##########
# make GRanges for LogFoldChanges for bigwig and bed
##########
#remove nas
resGR<-GenomicRanges::GRanges(seqnames=resLFC$chr,
                              IRanges::IRanges(start=resLFC$start,
                                               end=resLFC$end),
                              strand=resLFC$strand)
seqlengths(resGR)<-seqlengths(Celegans)[1:6]
mcols(resGR)<-resLFC[,c("wormbase","log2FoldChange","padj")]

names(mcols(resGR))[names(mcols(resGR))=="log2FoldChange"]<-"score"
resGR<-sort(resGR,ignore.strand=TRUE)

#https://github.com/hochwagenlab/hwglabr2/wiki/Apply-function-to-GRanges-scores-in-genomic-tiles
#######
### make bed file for significant genes
#######
idx<-which(resGR$padj<0.05)
forBed<-resGR[idx]
mcols(forBed)<-forBed$score
names(mcols(forBed))<-"score"
export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,"TEVcs_wt_lfc_p",gsub("^0.","",Threshold),".bed"),format="bed")


###########
# MAplot ALL genes
############

pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"MAplots_results.pdf"), width=5,height=5,paper="a4")


plotMA(res, main="DESeq2", ylim=c(-3,3),alpha=Threshold)
#plotCounts(dds, gene=which.min(res$padj), intgroup="sampleType")
plotMA(resLFC, main="DESeq2", ylim=c(-3,3),alpha=Threshold)
#plotCounts(dds, gene=which.min(res$padj), intgroup="sampleType")


#############
# MAplot X chr genes
#############
sink(file=paste0(outPath,"/txt/",fileNamePrefix,"logfile.txt"),append=TRUE, type="output")

chrXgenes<-mcols(dds)$gene[mcols(dds)$chr=="chrX"]
chrXres<-resLFC[rownames(resLFC) %in% chrXgenes,]

chrXres05<-chrXres[chrXres$padj<Threshold,]

cat("Using any log2FoldChange and padj smaller than",Threshold,":\n")

cat("Number of chrX genes with a significant increase in expression: ")
cat(sum(chrXres05$log2FoldChange>0),"\n")
#274
cat("Number of chrX genes with a significant decrease in expression: ")
cat(sum(chrXres05$log2FoldChange<0),"\n")
#40

upOnX<-chrXres05[chrXres05$log2FoldChange>0,]
write.table(rownames(upOnX), file=paste0(outPath,"/txt/",fileNamePrefix,"DCCgenes_JulieRNAseq.xls"),
            row.names=FALSE,col.names=FALSE)

plotMA(chrXres,main="chrX genes",ylim=c(-4,4),alpha=Threshold)
#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]


#############
# MAplotautosomal genes
#############

autosomalGenes<-mcols(dds)$gene[mcols(dds)$chr %in% c("chrI","chrII","chrIII","chrIV","chrV")]
autosomalRes<-resLFC[rownames(resLFC) %in% autosomalGenes,]

autosomalRes05<- autosomalRes[autosomalRes$padj<Threshold,]

cat("Number of autosomal genes with a significant increase in expression: ")
cat(sum(autosomalRes05$log2FoldChange>0),"\n")
#309
cat("Number of autosomal genes with a significant decrease in expression: ")
cat(sum(autosomalRes05$log2FoldChange<0), "\n")
#695

plotMA(autosomalRes,main="Autosomal genes",ylim=c(-4,4),alpha=Threshold)
#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]
dev.off()

#############
# Fisher test of number of up and down genes on X v autosomes
#############

upVdownXvA<-matrix(data=c(sum(chrXres05$log2FoldChange>0),
                          sum(chrXres05$log2FoldChange<0),
                          sum(autosomalRes05$log2FoldChange>0),
                          sum(autosomalRes05$log2FoldChange<0)),nrow=2,dimnames=list(group=c("Up","Down"),chr=c("chrX","chrA")))

cat("\nFisher Test, up v down:\n")
print(upVdownXvA)
print(fisher.test(upVdownXvA))


#############
# Fisher test of number of differentially expressed genes on X v autosomes
#############

testEnrich<-matrix(c(dim(chrXres)[1],dim(chrXres05)[1],
                     dim(autosomalRes)[1],
                     dim(autosomalRes05)[1]),
                   nrow=2,dimnames=list(group=c("NumTotal","NumSig"),chr=c("chrX","chrA")))
cat("\nFisher Test, enrichment of differentially expressed genes:\n")
print(testEnrich)
print(fisher.test(testEnrich))
sink()

#############
# Box plot by X v autosomes
#############
pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"boxPlots_expnByChrType.pdf"), width=5,height=5,paper="a4")

chrType<-factor(rownames(resLFC) %in% chrXgenes)
levels(chrType)<-c("Autosomal","X chr")
geneCounts<-table(chrType)

boxplot(log2FoldChange~chrType, data=resLFC, varwidth=TRUE, outline=FALSE, notch=TRUE,
        main="Expression changes after cleavage of dpy-26", col="grey", ylab="log2 Fold Change",
        xlab="chromosome type (number of genes)", names=paste(names(geneCounts)," \n(",geneCounts,")",sep=""))
#stripchart(log2FoldChange~chrType,data=res,method="jitter",vertical=TRUE,pch=20,col="#11115511",cex=0.5,add=TRUE)
abline(h=0,lty=2,col="blue")
dev.off()


#############
# Box plot by chromosome
#############
pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"boxPlots_expnByChr.pdf"), width=8,height=5,paper="a4")
chrName<-factor(resLFC$chr)
geneCounts<-table(chrName)

boxplot(log2FoldChange~chrName,data=resLFC,varwidth=TRUE,outline=FALSE,notch=TRUE,
        main="Expression changes after cleavage of dpy-26", ylab="log2 Fold Change",
        col=c(rep("grey",5),"purple"),xlab="chromosome (number of genes)",
        names=paste(names(geneCounts)," \n(",geneCounts,")",sep=""))
#stripchart(log2FoldChange~chrType,data=res,method="jitter",vertical=TRUE,pch=20,col="#11115511",cex=0.5,add=TRUE)
abline(h=0,lty=2,col="blue")
dev.off()


#############
# Volcano plot
#############
# https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"volcanoPlot_allBlack.pdf"), width=8,height=6,paper="a4")
EnhancedVolcano(resLFC,
                lab=rownames(resLFC),
                x="log2FoldChange",y="padj",
                selectLab=rownames(resLFC)[12366],
                xlim=c(-5.5,5.5),
                ylim=c(0,65),
                title= "TEVcs vs wildtype dpy-26",
                subtitle=NULL,
                caption = paste0(nrow(resLFC), ' expressed genes'),
                captionLabSize = 12,
                pCutoff=Threshold,
                pLabellingCutoff=10e-6,
                FCcutoff=1.0,
                xlab=bquote(~Log[2]~'fold change TEV+DPY-26cs / TEV+DPY-26wt'),
                ylab=bquote(~-Log[10]~adjusted~italic(P)),
                transcriptPointSize = 1.5,
                transcriptLabSize = 0,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 5.0,
                axisLabSize=14,
                col = c("grey20", "grey20", "grey20", "grey20"),
                colAlpha=0.8)
dev.off()





resByChr<-resLFC[order(resLFC$chr),]
# create custom key-value pairs for 'low', 'chrX', 'autosome' expression by fold-change
# set the base colour as 'black'
keyvals <- rep('black', nrow(resByChr))
# set the base name/label as 'NS'
names(keyvals) <- rep('NS', nrow(resByChr))
# modify keyvals for variables with fold change > 2.5
keyvals[which(resByChr$chr=="chrX")] <- 'red2'
names(keyvals)[which(resByChr$chr=="chrX")] <- 'chrX'

# modify keyvals for variables with fold change < -2.5
keyvals[which(resByChr$chr!="chrX")] <- 'royalblue'
names(keyvals)[which(resByChr$chr!="chrX")] <- 'autosomes'

# modify keyvals for variables with fold change < -2.5
#keyvals[which(abs(resByChr$log2FoldChange)<1 | resByChr$padj>5*10e-3)] <- 'grey10'
#names(keyvals)[which(abs(resByChr$log2FoldChange)<1 | resByChr$padj>5*10e-3)] <- 'NS'



pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"volcanoPlot_expnByChr.pdf"), width=8,height=6,paper="a4")
EnhancedVolcano(resByChr,
                lab=rownames(resByChr),
                x="log2FoldChange",y="padj",
                selectLab=rownames(resByChr)[12366],
                xlim=c(-5.5,5.5),
                ylim=c(0,65),
                title= "TEVcs vs wildtype dpy-26",
                subtitle=NULL,
                caption = paste0(nrow(resByChr), ' expressed genes'),
                captionLabSize = 12,
                pCutoff=Threshold,
                pLabellingCutoff=10e-6,
                FCcutoff=1.0,
                xlab=bquote(~Log[2]~'fold change TEV+DPY-26cs / TEV+DPY-26wt'),
                ylab=bquote(~-Log[10]~adjusted~italic(P)),
                transcriptPointSize = 1.5,
                transcriptLabSize = 0,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 5.0,
                axisLabSize=14,
                colCustom=keyvals,
                colAlpha=0.8)
dev.off()

idx<-resByChr$chr=="chrX"
pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"volcanoPlot_chrX.pdf"), width=8,height=6,paper="a4")
EnhancedVolcano(resByChr[idx,],
                lab=rownames(resByChr[idx,]),
                x="log2FoldChange",y="padj",
                selectLab=rownames(resByChr)[12366],
                xlim=c(-5.5,5.5),
                ylim=c(0,65),
                title= "TEVcs vs wildtype dpy-26 - chrX genes",
                subtitle=NULL,
                caption = paste0(nrow(resByChr[idx,]), ' expressed genes'),
                captionLabSize = 12,
                pCutoff=Threshold,
                pLabellingCutoff=10e-6,
                FCcutoff=1.0,
                xlab=bquote(~Log[2]~'fold change TEV+DPY-26cs / TEV+DPY-26wt'),
                ylab=bquote(~-Log[10]~adjusted~italic(P)),
                transcriptPointSize = 1.5,
                transcriptLabSize = 0,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 5.0,
                axisLabSize=14,
                colCustom=keyvals[idx],
                colAlpha=0.8)
dev.off()



idx<-resByChr$chr!="chrX"
pdf(file=paste0(outPath,"/plots/",fileNamePrefix,"volcanoPlot_autosomes.pdf"), width=8,height=6,paper="a4")
EnhancedVolcano(resByChr[idx,],
                lab=rownames(resByChr[idx,]),
                x="log2FoldChange",y="padj",
                selectLab=rownames(resByChr)[12366],
                xlim=c(-5.5,5.5),
                ylim=c(0,65),
                title= "TEVcs vs wildtype dpy-26 - chrX genes",
                subtitle=NULL,
                caption = paste0(nrow(resByChr[idx,]), ' expressed genes'),
                captionLabSize = 12,
                pCutoff=Threshold,
                pLabellingCutoff=10e-6,
                FCcutoff=1.0,
                xlab=bquote(~Log[2]~'fold change TEV+DPY-26cs / TEV+DPY-26wt'),
                ylab=bquote(~-Log[10]~adjusted~italic(P)),
                transcriptPointSize = 1.5,
                transcriptLabSize = 0,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 5.0,
                axisLabSize=14,
                colCustom=keyvals[idx],
                colAlpha=0.8)
dev.off()


pThresh=0.05
LFCthresh=0
summaryByChr<-function(resLFC,pThresh,LFCthresh) {
  up<-resLFC[resLFC$padj < pThresh & resLFC$log2FoldChange > LFCthresh,]
  down<-resLFC[resLFC$padj < pThresh & resLFC$log2FoldChange < -LFCthresh, ]
  allChr<-as.data.frame(rbind(up=table(up$chr),down=table(down$chr)))
  allChr$autosomes<-rowSums(allChr[,1:5])
  allChr$total<-rowSums(allChr[,1:6])
  rownames(allChr)<-paste0(rownames(allChr),"_p",pThresh,"_lfc",LFCthresh)
  return(allChr)
}


sink(file=paste0(outPath,"/txt/", fileNamePrefix,"logfile.txt"),append=TRUE, type="output")
cat("Summary by Chr: \n")
cat("\np=0.05, LFC=0: \n")
print(summaryByChr(resLFC,pThres=0.05,LFCthresh=0))
cat("\np=0.05, LFC=0.5: \n")
print(summaryByChr(resLFC,pThres=0.05,LFCthresh=0.5))
cat("\np=0.05, LFC=1: \n")
print(summaryByChr(resLFC,pThres=0.05,LFCthresh=1))

cat("\np=0.01, LFC=0: \n")
print(summaryByChr(resLFC,pThres=0.01,LFCthresh=0))
cat("\np=0.01, LFC=0.5: \n")
print(summaryByChr(resLFC,pThres=0.01,LFCthresh=0.5))
cat("\np=0.01, LFC=1: \n")
print(summaryByChr(resLFC,pThres=0.01,LFCthresh=1))

sink()
#summary(resLFC,alpha=0.05)


