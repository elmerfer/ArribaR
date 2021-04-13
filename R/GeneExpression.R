#Gene Expression
#' GetCounts
#' it estimate the gene counts wrapping the Rsubread::featureCounts function
#' @usage ft <- GetCounts("my.aligned_bam", feature = "exon")
#' @param sbjBAM the BAM file name
#' @param feature string indicating if exons or genes are returned
#' @param nthreads integer, total of treads to use, default parallel::detectCores() -1
#' @export
#' 
GetCounts <- function(sbjBAM, feature = c("exon","gene", "transcript"), nthreads ){
  software <- .OpenConfigFile()
  if(!file.exists(sbjBAM)){
    stop("ERROR: sbjt not found")
  }
  if(missing(nthreads)){
    nthreads <- parallel::detectCores()-4
  }
  
  geneAnnot <- function(gID){
    gInfo <- SAF[which(SAF$GeneID == gID),]
    geneName <- unique(gInfo$geneName)
    if(length(geneName) >1){
      stop(paste("Error geneName length",geneName))
    }
    transcript <- paste(paste(gInfo$transcript,gInfo$exonNumber,sep="_"),collapse = ";")
    return(c(geneName,transcript))
  }
  
  stopifnot(nthreads>0)
  stopifnot(file.exists(sbjBAM))
  feature <- match.arg(feature[1], choices = c("exon","gene"))
  
  # SAF <- subset(ArribaR:::.ReadExonAnnotationFile(),source == "ENSEMBL" & type == "exon")
  software <- .OpenConfigFile()
  message("Reading gen annotation")
  SAF <- rtracklayer::readGFF(file.path(software$annotation,"GENCODE19.gtf"))
  
  SAF <- subset(SAF, type == feature & gene_type =="protein_coding")
  if(feature == "exon"){
    message("Getting protein coding and experimentally confirmed exon annotation")
    SAF <- subset(SAF, tag %in% c("CCDS","exp_conf"))
    SAF$GeneID <- SAF$gene_id
  }else{
    
    SAF$GeneID <- SAF$gene_id
  }
  
  colnames(SAF)[colnames(SAF) %in% c("seqid","start","end","strand")] <- c("Chr","Start","End","Strand")
  SAF$Chr <- paste0("chr",SAF$Chr)
  
  require(org.Hs.eg.db)
  ft <- Rsubread::featureCounts(files = sbjBAM,
                                 annot.ext = SAF[,c("Chr","Start","End","Strand","GeneID")],
                                 isGTFAnnotationFile = FALSE,
                                # annot.inbuilt = "hg19",
                                useMetaFeatures = ifelse(feature == "exon", FALSE, TRUE),
                                isPairedEnd = TRUE,
                                nthreads = nthreads)
  message("adding GeneSymbol annotation from org.Hs.eg Bioconductor library ")
  ft$annotation <- cbind(ft$annotation, SAF[,which(!c(colnames(SAF) %in% colnames(ft$annotation)))])
  # rownames(ft$annotation) <- 
  #   comb <-paste0(ft$annotation$GeneID,ft$annotation$Start,ft$annotation$End,ft$annotation$Chr)
  #   comb2 <- paste0(SAF$GeneID,SAF$Start,SAF$End,SAF$Chr)
   # ft$annotation <- merge(ft$annotation, SAF[,c("Chr","GeneID","Start","End","geneName","exonNumber","attributes")],by= "GeneID",sort=FALSE, all.x = TRUE)
   # if(feature == "gene"){
   #   geneann <- do.call(rbind,lapply(ft$annotation$GeneID, geneAnnot))
   #   colnames(geneann) <- c("GeneName","Transcripts")
   #   ft$annotation <- cbind(ft$annotation,geneann)
   # }else{
   #   ft$annotation <- merge(ft$annotation, SAF[,c("Chr","GeneID","Start","End","geneName","transcript","exonNumber","attributes")],sort=FALSE, all.x = TRUE)
   # }
  return(ft)
}

