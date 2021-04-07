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
  stopifnot(nthreads>0)
  stopifnot(file.exists(sbjBAM))
  feature <- match.arg(feature[1], choices = c("exon","gene"))
  
  SAF <- subset(ArribaR:::.ReadExonAnnotationFile(),source == "ENSEMBL" & type == "exon")
  SAF$Chr <- paste0("chr",SAF$contig)
  colnames(SAF)[colnames(SAF) %in% c("contig","start","end","strand","geneID")] <- c("Chr","Start","End","Strand","GeneID")
  
  ft2 <- Rsubread::featureCounts(files = sbjBAM,
                                annot.ext = SAF,
                                isGTFAnnotationFile = FALSE,
                                useMetaFeatures = ifelse(feature == "exon", FALSE, TRUE),
                                isPairedEnd = TRUE,
                                nthreads = 8)
  if(feature == "gene"){
    ft$annotation <- merge(ft$annotation, unique(SAF[,c("GeneID","geneName","attributes")]), by.x = "GeneID", by.y = "GeneID",sort=FALSE, all = FALSE)
  }else{
    ft$annotation <- merge(ft$annotation, SAF[,c("GeneID","geneName","attributes")], by.x = "GeneID", by.y = "GeneID",sort=FALSE, all = FALSE)
  }
  return(ft)
}

