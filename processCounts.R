# ff <- list.files( path = "./counts", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
# counts.files <- lapply( ff, read.table, skip = 4 )
# counts <- as.data.frame( sapply( counts.files, function(x) x[ , number ] ) )
# ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff )
# ff <- gsub( "[.]/counts/", "", ff )
# colnames(counts) <- ff
# row.names(counts) <- counts.files[[1]]$V1


library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(ggplot2)
library("RColorBrewer")
library("PoiClaClu")
library("pheatmap")

# get funciton for converting gene names from WBID to publicID
source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")

###############################################################
### some variables
###############################################################
genomeVer="ws260"

sampleList<-read.csv("../sampleList.csv",stringsAsFactors=F)


###############################################################
### create tx2gene object
###############################################################

# create a txdb from wormbase data
# gffFile=paste0("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/",genomeVer,"/c_elegans.PRJNA13758.",toupper(genomeVer),".annotations.gff3.gz")
# if (grepl(".gz",gffFile)) {
#   system(paste("gunzip ",gffFile))
#   gffFile<-gsub(".gz","",gffFile)
# }
# txdb<-makeTxDbFromGFF(gffFile,format="gff3",dataSource="wormbase",organism="Caenorhabditis elegans")
# saveDb(txdb, file=paste0("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/",
#                                genomeVer,"/",genomeVer,".sqlite"))

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/",
              genomeVer,"/",genomeVer,".sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")


###############################################################
### get samples into DESeq2
###############################################################
# modify sampleList
samples<-sampleList
#samples$starCountFiles<-gsub(".fastq.gz",".ReadsPerGene.out.tab",basename(as.vector(sampleList$fastqFile)))
samples$salmonCountFiles<-gsub(".fastq.gz","/quant.sf",basename(as.vector(sampleList$fastqFile)))
samples<-samples[,-which(names(samples)=="fastqFile")]

# add path to salmon count files
files<-paste0("salmon/mRNA/",samples$salmonCountFiles)
# get the sample name (dirname in salmon file structure)
names(files)<-dirname(samples$salmonCountFiles)

# import the count matrices
txi<-tximport(files,type="salmon",tx2gene=tx2gene)

# make factors from the variables of interest for comparisons
samples$hlh1exp<-factor(ifelse(samples$hs*samples$hsHLH1==1,"HLH1exp","noHLH1"),levels=c("noHLH1","HLH1exp"))
samples$mes2<-factor(ifelse(samples$mes2==1,"mes2","wt"),levels=c("wt","mes2"))
samples$fed<-factor(ifelse(samples$fed==1,"fed","starved"),levels=c("fed","starved"))
samples$hs<-factor(ifelse(samples$hs==1,"hs","noHs"),levels=c("noHs","hs"))
#batch<-apply( samples[ , c("lane","sampleID") ] , 1 , paste , collapse = "-" )
samples$lane<-factor(samples$lane)
samples$sampleID<-factor(samples$sampleID)
samples$strain<-factor(samples$strain)
samples$hsHLH1<-factor(ifelse(samples$hsHLH1==1,"hsHLH1","wt"),levels=c("wt","hsHLH1"))

# read samples into DESeq2
dds <- DESeqDataSetFromTximport(txi, samples,~lane+hs+fed+hlh1exp+mes2+
                                  fed:mes2+hlh1exp:mes2+
                                  fed:mes2:hlh1exp)


###############################################################
### DESeq2 differential expression analysis (using negative binomial distribution)
###############################################################

#dds<-collapseReplicates(dds,groupby=samples$sampleID,renameCols=T)
dds <- DESeq(dds)
# This function performs a default analysis through the steps:
#   1. estimation of size factors: estimateSizeFactors
#   2. estimation of dispersion: estimateDispersions
#   3. Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
#   returns a DESeqDataSet object


#####################
### heatmap sampleVgene
#####################
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:200]
df <- as.data.frame(colData(dds)[,c("hsHLH1","hs","hlh1exp","mes2","fed")])
ntd <- normTransform(dds)  # does simple log transform + 1 pseudocounts
if(!dir.exists("plots")) {
  dir.create("plots")
}

pdf(paste0("plots/heatmap_sampleVgene.pdf"),paper="a4r",height=8, width=11)
pheatmap(assay(ntd)[select,], cluster_rows=F, show_rownames=F,
         cluster_cols=TRUE, annotation_col=df, fontsize_col=8,
         main="sample clustering with top 200 expressed genes")
dev.off()

###############################################################
### variance stabilizing transformation for clustering and PCA
###############################################################
vsd <- vst(dds, blind=TRUE) # unbiased transformation wihtout prior knowledge of samples

#####################
### hclust
#####################
## Euclidian distance between sample (after vsd)
pdf(paste0("plots/heatmap_sampleVsample",".pdf"),paper="a4r",height=8, width=11)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(samples$sampleID,samples$mes2,samples$hlh1exp,samples$fed,samples$hs,samples$hsHLH1,sep="_")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,main="Euclidian sample-sample distances")

## Poisson distance between samples
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste(samples$sampleID,samples$mes2,samples$hlh1exp,samples$fed,samples$hs,samples$hsHLH1,sep="_")
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,main="Poisson sample-sample distances")
dev.off()


#####################
### pca
#####################
par(mfrow=c(4,2))
pdf(paste0("plots/pca_samples.pdf"),paper="a4",height=11, width=8)
# do pca plots colouring by different factors
plotPCA(vsd, intgroup=c("lane"))
plotPCA(vsd, intgroup=c("sampleID"))
plotPCA(vsd, intgroup=c("strain"))
plotPCA(vsd, intgroup=c("hsHLH1"))
plotPCA(vsd, intgroup=c("hs"))
plotPCA(vsd, intgroup=c("hlh1exp"))
plotPCA(vsd, intgroup=c("mes2"))
plotPCA(vsd, intgroup=c("fed"))
dev.off()
par(mfrow=c(1,1))
### for plotting more than 1 factor
# pcaData <- plotPCA(vsd, intgroup=c("hlh1exp", "fed"), returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(PC1, PC2, color=hlh1exp, shape=fed)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   coord_fixed()


#####################
### dispersion
#####################
pdf(paste0("plots/dispersion.pdf"),paper="a4",height=11, width=8)
plotDispEsts(dds,main="dispersion")
dev.off()

###############################################################
### look at specific contrasts
###############################################################

if(!dir.exists("txt")) {
  dir.create("txt")
}
if(!dir.exists("csv")) {
  dir.create("csv")
}

allCoefs<-resultsNames(dds) # lists the coefficients
idx<-match(c("hlh1exp_HLH1exp_vs_noHLH1","mes2_mes2_vs_wt",
             "fedstarved.mes2mes2","hlh1expHLH1exp.mes2mes2",
             "fedstarved.hlh1expHLH1exp.mes2wt",
             "fedstarved.hlh1expHLH1exp.mes2mes2"),allCoefs)

pdf(paste0("plots/contrasts.pdf"),paper="a4",height=11, width=8)
par(mfrow=c(3,1))

# loop through all the different contrasts of interest
for (i in idx) {
  pthresh=0.05
  lfcthresh=0.5
  i=idx[1]
  # get results table
  res <- results(dds,name=allCoefs[i])

  # print summary to file
  cat(allCoefs[i],file=paste0("txt/summary_",allCoefs[i],".txt"),sep="\n")
  capture.output(summary(res,alpha=thresh),file=paste0("txt/summary_",allCoefs[i],".txt"),append=T)

  # get shrunken logFC estimates
  resAsh <- lfcShrink(dds, coef=i, type="ashr")
  # print summary to file
  capture.output(summary(resAsh,alpha=thresh),file=paste0("txt/summary_",allCoefs[i],".txt"),append=T)

  #####################
  ### plots
  #####################
  # make MA plot
  plotMA(resAsh, ylim=c(-3,3),alpha=thresh,main=allCoefs[i])

  # plot histograms of p values
  hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
       col = "grey50", border = "white")
  hist(res$padj[res$baseMean > 1], breaks = 0:20/20,
       col = "grey50", border = "white")


  #####################
  ### data tables
  #####################
  # order file by pvalue and save
  resAshOrdered <- as.data.frame(resAsh[order(resAsh$padj),])
  IDS<-convertGeneNames(rownames(resAshOrdered),inputType="WBgeneID",
                     outputType=c("seqID","WBgeneID","publicID"))
  resAshOrdered<-cbind(resAshOrdered,IDS)
  # write to file
  write.csv(resAshOrdered, file = paste0("csv/results_",allCoefs[i],"_all.csv"))

  # extract data for genes of interest (GOI) and write to file
  GOI<-c("lin-12","lag-2")
  j<-match(GOI,resAshOrdered$publicID)
  cat(allCoefs[i],file=paste0("txt/GOI.txt"),sep="\n",append=T)
  cat(paste(colnames(resAshOrdered),collapse="\t"),file=paste0("txt/GOI.txt"),sep="\n",append=T)
  cat(paste(apply(resAshOrdered[j,],1,paste,collapse="\t"),collapse="\n"),file=paste0("txt/GOI.txt"),sep="\n",append=T)

  # filter for best results
  resAshOrdered <- resAshOrdered[resAshOrdered$padj<pthresh & abs(resAshOrdered$log2FoldChange)>lfcthresh,]
  write.csv(resAshOrdered, file = paste0("csv/results_",allCoefs[i],"_p",pthresh,"_lfc",lfcthresh,".csv"))
}


