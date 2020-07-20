
#' Filter DESeq2 table results
#'
#' @param resultsTable Table of DESeq results
#' @param padj Adjusted p value threshold
#' @param lfc Log fold change threshold
#' @param direction Whether to find genes that are less than (lt), or greater than (gt) the log fold change threshold, or both extreme tails ("both")
#' @param chr Include all genes in genome ("all") only those on the X chromosome ("chrX"), or only autosomes ("autosomes")
#' @param outPath Path to working directory
#' @return filtered table of results which is also automatically written to disk
#' @export
filterResults<-function(resultsTable, padj=0.05, lfc=0, direction="both",
                        chr="all",outPath=".") {
  if(direction=="both") {
    idx<-!is.na(resultsTable$padj) & resultsTable$padj<padj & abs(resultsTable$log2FoldChange)>lfc
  } else if(direction=="gt") {
    idx<-!is.na(resultsTable$padj) & resultsTable$padj<padj & resultsTable$log2FoldChange>lfc
  } else if(direction=="lt") {
    idx<-!is.na(resultsTable$padj) & resultsTable$padj<padj & resultsTable$log2FoldChange<lfc
  } else {
    print("direction must be 'both' to get both tails, \n'gt' to get lfc larger than a specific value, \nor 'lt' to get lfc less than a certain value")
  }
  if(chr=="all"){
    idx<-idx
  } else if(chr=="chrX"){
    idx<-idx & !is.na(resultsTable$chr) & resultsTable$chr=="chrX"
  } else if(chr=="autosomes"){
    idx<-idx & !is.na(resultsTable$chr) & resultsTable$chr!="chrX"
  } else {
    print("chr must be one of 'all', 'chrX' or 'autosomes'")
  }
  filtTable<-resultsTable[idx,c("baseMean","log2FoldChange","padj",
                                "wormbase","chr","start","end","strand")]
  if(!dir.exists(paste0(outPath,"/txt"))){
    dir.create(paste0(outPath,"/txt"))
  }
  write.csv(filtTable,file=paste0(outPath,"/txt/filtResults_p",
                                  padj,"_",direction,"-","lfc",
                                  lfc,"_",chr,".csv"), row.names=F,
            quote=F)
  return(filtTable)
}
