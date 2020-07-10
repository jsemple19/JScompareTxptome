library(QuasR)
library(BSgenome)
library(Rsamtools)
library(GenomicFeatures)
library(Gviz)
library(parallel)

#setwd("/data1/projects/p025/RNA_seq_Julie/")
getwd()

#TxDb.Celegans.UCSC.ce10.ensGene <- makeTxDbFromUCSC(genome = "ce10", tablename = "ensGene")
#txdb.ce10  <- makeTxDbFromUCSC(genome = "ce10", tablename = "ensGene")
#txdb.ws240 <- loadDb("/data1/projects/p025/TxDb/txdb_ws240.sqlite")
#txdb.ws250 <- loadDb("/data1/projects/p025/TxDb/txdb_ws250.sqlite")


txdb <- c(txdb.ce10, txdb.ws240, txdb.ws250)

sampleFile <- "samples.txt"
#sampleFile <- "1_sample.txt"
genomeFile <- "../genome_ce10/genome.fa"
exonFile <- c("exons-ce10-3.tab", "exons-ws240-3.tab", "exons-ws250-3.tab")
geneFile <- c("gene-ce10-3.tab", "gene-ws240-3.tab", "gene-ws250-3.tab")


# define cores
cl <- makeCluster(8)

proj <- qAlign(sampleFile, genomeFile, splicedAlignment=TRUE, paired="no", clObj=cl)

save(proj, file="proj-3.R")
qQCReport(proj, pdfFilename="qc_report-3.pdf")

#qExportWig(proj, binsize=10L, scaling=TRUE, collapseBySample=TRUE)


  exonLevels <- qCount(proj, txdb.ce10, reportLevel="exon", orientation="opposite", selectReadPosition="end")
  write.table(exonLevels, file=exonFile[1], sep="\t")
  geneLevels <- qCount(proj, txdb.ce10, reportLevel="gene", orientation="opposite", selectReadPosition="end")
  write.table(geneLevels, file=geneFile[1], sep="\t")

  exonLevels <- qCount(proj, txdb.ws240, reportLevel="exon", orientation="opposite", selectReadPosition="end")
  write.table(exonLevels, file=exonFile[2], sep="\t")
  geneLevels <- qCount(proj, txdb.ws240, reportLevel="gene", orientation="opposite", selectReadPosition="end")
  write.table(geneLevels, file=geneFile[2], sep="\t")

  exonLevels <- qCount(proj, txdb.ws250, reportLevel="exon", orientation="opposite", selectReadPosition="end")
  write.table(exonLevels, file=exonFile[3], sep="\t")
  geneLevels <- qCount(proj, txdb.ws250, reportLevel="gene", orientation="opposite", selectReadPosition="end")
  write.table(geneLevels, file=geneFile[3], sep="\t")

####### realign to newer annotation to compare to star and salmon

  genomeVer="WS275"
  genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)
  genomeFile=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer,"/sequence/c_elegans.PRJNA13758.",genomeVer,".genomic.fa")
  txdb<-loadDb(paste0(genomeDir, "/annotations/c_elegans.PRJNA13758.",
                      genomeVer, ".annotations.sqlite"))


  library(QuasR)
  cl <- makeCluster(2)
  sampleFile <- "samples.txt"
  proj <- qAlign(sampleFile, genomeFile, splicedAlignment=TRUE, paired="no",
                 clObj=cl)

  save(proj, file="proj-4.R")
  qQCReport(proj, pdfFilename="qc_report-4.pdf")

  geneLevels <- qCount(proj, txdb, reportLevel="gene", orientation="opposite",
                       selectReadPosition="end")
  write.table(geneLevels, file="quasr-gene-WS275.tab", sep="\t")
