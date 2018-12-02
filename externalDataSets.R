
library(tidyr)

# get funciton for converting gene names from WBID to publicID
source("~/Documents/MeisterLab/GenomeVer/geneNameConversion/convertingGeneNamesFunction1.R")



###############################
## get stage specific expression from Boeck-Waterston_GR2016
###############################

if(!dir.exists("externalData")) {
  dir.create("externalData")
}
tcFile="expressionTC_Boeck-Waterston_GR2016"
if(!file.exists(paste0("externalData/",tcFile))) {
  link="https://genome.cshlp.org/content/suppl/2016/09/20/gr.202663.115.DC1/Supplemental_Table_S2.gz"
  download.file(link,paste0("externalData/",tcFile,".gz"))
  system(paste0("gunzip -S txt externalData/",tcFile,".gz"))
}
# Exact times for growth after plating of starved L1s (all at 25°C) were as follows:
# (1) L1: worms were grown 4.0 h post-L1 plating;
# (2) L2: for 17.75 h (Pn.p cells visible but not divided, gonad just starting to proliferate);
# (3) L3: for 26.75 h (Pn.p cells divided once or twice; gonad
# just starting to turn up);
# (4) L4: for 34.25 h (vulvae are in Christmas tree stage, gonad has passed bend, sperm are present);
# (5) young adult: for 46 h (vulvae fully formed and oocytes present in gonad, but no embryos);
# (6) dauer entry: daf-2(e1370) 48 h post-L1 stage larvae;
# (7) dauer: daf-2(e1370) 91 h post-L1 stage larvae;
# (8) dauer exit: daf-2(e1370) at 25°C for 91 h and at 15°C for 12 h;
# male L4: him-8(e1480) mid-L4 30 h post-L1 stage larvae (filtered through mesh to purify males);
# soma: JK1107(glp-1(q224)) mid-L4 30 h post-L1 stage larvae.
# Dissected gonads were from N2 (wild type) animals grown for 48 h at 20°C post-L1 stage larvae; approximately
# 200 gonads dissected and isolated from carcasses.
tcData<-read.table(paste0("externalData/",tcFile),stringsAsFactors=F,header=T)
countCols<-grep("_counts",names(tcData))
dataSets<-list(EMB=c("N2_EE_50.600_counts","X20120223_EMB.600_counts"),
              L1=c("L1_counts","N2_L1.1_counts"),
              L2=c("L2_counts","N2_L2.4_counts"),
              L3=c("L3_counts","N2_L3.1_counts"),
              L4=c("L4_counts","L4b_counts"),
              YA=c("YA_counts","N2_Yad.1_counts"),
              Dentry=c("DauerEntryDAF2_counts","DauerEntryDAF2.1.1_counts",
                        "DauerEntryDAF2.2_counts","DauerEntryDAF2.4.1_counts"),
              DA=c("DauerDAF2.2_counts","DauerDAF2.2.1_counts","DauerDAF2.5.1_counts"))

processDataset<-function(sourceDF,sampleName){
  p<-ggpairs(log2(sourceDF+1))
  ggsave(paste0("externalData/",tcFile,"_plots/",sampleName,"_Scatter.pdf"),plot=p)
  avgDF<-data.frame(rowMeans(sourceDF))
  names(avgDF)<-sampleName
  return(avgDF)
}

if(!dir.exists(paste0("externalData/",tcFile,"_plots"))) {
  dir.create(paste0("externalData/",tcFile,"_plots"))
}
tcAvgDF<-data.frame(matrix(ncol=length(dataSets),nrow=nrow(tcData)))
names(tcAvgDF)<-names(dataSets)
for (d in 1:length(dataSets)) {
  tcAvgDF[,names(dataSets)[d]]<-processDataset(sourceDF=tcData[,dataSets[[d]]],
                                sampleName=names(dataSets)[d])
}
tcAvgDF$WormbaseName<-tcData$WormbaseName
write.csv(tcAvgDF,file=paste0("externalData/",tcFile,"_plots/stageCounts.csv"),row.names=F)


###############################
## get germline-soma genes from Boeck-Waterston_GR2016
###############################

glCols<-c("L4JK1107soma_counts","L4JK1107soma.2_counts","N2_Ad_gonad.1.RZLI_counts")
glCounts<-as.matrix(tcData[,glCols])
IDS<-convertGeneNames(tcData$WormbaseName,inputType="seqID",
                      outputType=c("seqID","WBgeneID","publicID"))
row.names(glCounts)<-IDS$WBgeneID
coldata<-data.frame(sampleNames=glCols,germline=factor(c("soma","soma","germline")))

glCounts<-glCounts[!is.na(row.names(glCounts)),]


dds <- DESeqDataSetFromMatrix(countData = glCounts,
                              colData = coldata,
                              design = ~ germline)

featureData <- data.frame(gene=rownames(glCounts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

#remove genes with few reads
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]

# estimate parameters
dds <- DESeq(dds)
res <- results(dds,alpha=0.05)
summary(res)
resLFC<- lfcShrink(dds,contrast=c("germline","germline","soma"),type="ashr")
plotMA(resLFC)

resUp_gl<-resLFC[resLFC$padj<0.05 & resLFC$log2FoldChange>1,]
resDown_soma<-resLFC[resLFC$padj<0.05 & resLFC$log2FoldChange<(-1),]

glvSoma<-data.frame(WBgeneID=c(rownames(resUp_gl),rownames(resDown_soma)),
           germline=c(rep("germline",length(rownames(resUp_gl))), rep("soma",length(rownames(resDown_soma)))))
write.csv(glvSoma,paste0("externalData/",tcFile,"_plots/germlineSomaGenes.csv"))

###############################
## get germline-soma genes from Reinke_DEV2004
###############################

glFile="germlineVsoma_Reinke_Dev2004"
if (!file.exists(paste0("externalData/",glFile,"/Fig1\ I\ wt\ vs\ glp4\ enriched\ genes.txt"))) {
  link2="http://dev.biologists.org/highwire/filestream/1201187/field_highwire_adjunct_files/0/Data_S1.zip"
  download.file(link2,paste0("externalData/",glFile,".zip"))
  system(paste0("unzip externalData/",glFile,".zip -d externalData/",glFile))
}
glData<-read.delim(paste0("externalData/",glFile,"/Fig1\ I\ wt\ vs\ glp4\ enriched\ genes.txt"),stringsAsFactors=F,header=T,sep="\t")
glData<-glData[,c("WormbaseID","exclusive.category")]
