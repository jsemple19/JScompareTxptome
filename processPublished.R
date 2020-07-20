## 2020-07-16
## Some published data sets were collected. A big list of dosage
## compensated and non dosage compenstated genes was manually copied
## from PDF of supplementary data (supl tables 4 & 5).
## A shorter list of "classical" dosage compenstaed and escaper genes was
## curated manually from figures or text in papers, assuming that if they
## are used as examples they are higher confidence (BUT they are all from
## embryos...):
## Wheeler 2016 just mentions some control genes in text. pcr in supl fig
## Csankovszki 2004 has a few genes in a figure first identified in Meyer 1986
## Kruesi 2016 in supl fig 1 to fig 7 took DC genes from Jans and tested
## 27 (?) of them for GROseq changes +- SDC-2, i took those with FC>1.5
## Jans 2009 has a few genes it confirmed by qPCR in supl table 6

library(magrittr)


genomeVer="WS275"
genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)

# txdb<-AnnotationDbi::loadDb(paste0(genomeDir,
#                                    "/annotations/c_elegans.PRJNA13758.",
#                                    genomeVer, ".annotations.sqlite"))
# #columns(txdb) # what kind of data is retrievable
# #keytypes(txdb)
# k <- keys(txdb, keytype = "GENEID")
# geneChr <- AnnotationDbi::select(txdb, k, columns=c("CDSCHROM"),
#                                  keytype="GENEID")


srcref <- Organism.dplyr::src_organism("TxDb.Celegans.UCSC.ce11.refGene")
metadata<-dplyr::inner_join(dplyr::tbl(srcref, "id"),
                                   dplyr::tbl(srcref, "ranges_gene")) %>%
  dplyr::select(wormbase, alias, genename, gene_chrom,
                gene_start, gene_end, gene_strand) %>%
  dplyr::collect() %>% GenomicRanges::GRanges()



#######################
## manually curated from papers
#######################

pubDC<-data.table::fread(input="published_DC.txt")
pubNDC<-data.table::fread(input="published_Xescapers.txt")


pubDCgr<-metadata[metadata$wormbase %in% pubDC$wbid]

mcols(pubDCgr)<-cbind(mcols(pubDCgr),pubDC[match(pubDCgr$wormbase,pubDC$wbid),])
pubDCgr$wbid<-NULL
pubDCgr$alias<-NULL
pubDCgr

saveRDS(pubDCgr,file="published_DCgr.rds")


pubNDCgr<-metadata[metadata$wormbase %in% pubNDC$wbid]

mcols(pubNDCgr)<-cbind(mcols(pubNDCgr),pubNDC[match(pubNDCgr$wormbase,pubNDC$wbid),])
pubNDCgr$wbid<-NULL
pubNDCgr$alias<-NULL
pubNDCgr

saveRDS(pubNDCgr,file="published_NDCgr.rds")


#######################
## Jans 2009
#######################

JansDC<-data.table::fread(input="Jans2009_DC_suplTable4.txt")
JansNDC<-data.table::fread(input="Jans2009_notDC_suplTable5.txt")


JansDCgr<-metadata[metadata$wormbase %in% JansDC$WormBaseId]

mcols(JansDCgr)<-cbind(mcols(JansDCgr),JansDC[match(JansDCgr$wormbase,JansDC$WormBaseId),])
JansDCgr$WormBaseId<-NULL
JansDCgr$alias<-NULL
JansDCgr

saveRDS(JansDCgr,file="Jans2009_DCgr.rds")


JansNDCgr<-metadata[metadata$wormbase %in% JansNDC$WormBaseId]

mcols(JansNDCgr)<-cbind(mcols(JansNDCgr),JansNDC[match(JansNDCgr$wormbase, JansNDC$WormBaseId),])
JansNDCgr$WormBaseId<-NULL
JansNDCgr$alias<-NULL
JansNDCgr

saveRDS(JansNDCgr,file="Jans2009_NDCgr.rds")




#######################
## Kramer 2015
#######################
kramerURL<-"https://doi.org/10.1371/journal.pgen.1005698.s011"
kramerFileName="Kramer_2015_PlotGen_S3file.xlsx"
#download.file(url=kramerURL,destfile=kramerFileName)

kramer<-readxl::read_excel(kramerFileName,col_types=c(rep("text",3),rep("numeric",30)))
names(kramer)
kramer<-kramer[,c(1:3,grep("_L3_",names(kramer)))]

kramergr<-metadata[metadata$wormbase %in% kramer$Gene_WB_ID]

mcols(kramergr)<-cbind(mcols(kramergr),kramer[match(kramergr$wormbase, kramer$Gene_WB_ID),])
kramergr$Gene_WB_ID<-NULL
kramergr$alias<-NULL

saveRDS(kramergr,file="kramer2015_gr.rds")
kramergr<-readRDS(file="kramer2015_gr.rds")


lfcVal<-0.5
padjVal<-0.05
idx<-!is.na(kramergr$dpy27_RNAi_L3_padj) &
  kramergr$dpy27_RNAi_L3_padj < padjVal &
  kramergr$dpy27_RNAi_L3_log2_fold_change > lfcVal &
  seqnames(kramergr)=="chrX"
kramerdpy27dc<-kramergr[idx,]
colIdx<-grep("(set)|(dpy21)|(mixed_sex)",colnames(mcols(kramerdpy27dc)))
mcols(kramerdpy27dc)[,colIdx]<-NULL
saveRDS(kramerdpy27dc,file=paste0("kramer2015_chrXup_dpy27_lfc",
                             formatC(lfcVal,format="e",digits=0),"_p",
                             formatC(padjVal,format="e",digits=0),
                             "_gr.rds"))


idx<-!is.na(kramergr$dpy21_mutant_L3_padj) &
  kramergr$dpy21_mutant_L3_padj < padjVal &
  kramergr$dpy21_mutant_L3_log2_fold_change > lfcVal &
  seqnames(kramergr)=="chrX"
kramerdpy21dc<-kramergr[idx,]
colIdx<-grep("(set)|(dpy27)",colnames(mcols(kramerdpy21dc)))
mcols(kramerdpy21dc)[,c(colIdx)]<-NULL
saveRDS(kramerdpy21dc,file=paste0("kramer2015_chrXup_dpy21_lfc",
                                  formatC(lfcVal,format="e",digits=0),"_p",
                                  formatC(padjVal,format="e",digits=0),
                                  "_gr.rds"))




