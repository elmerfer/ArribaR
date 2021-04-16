#Gene Expression
#' GetCounts
#' it estimate the gene counts wrapping the Rsubread::featureCounts function
#' @usage ft <- GetCounts("my.aligned_bam", feature = "exon")
#' @param sbjBAM the BAM file name
#' @param feature string indicating if exons or genes are returned
#' @param annotSource string possible values are "GENCODE19" or "inbuilt". If GENECODE19, it will use GENOCDE from Arriba software, if "inbuilt" it will use as in \code{\link[Rsubread]{featureCounts}}
#' @param nthreads integer, total of treads to use, default parallel::detectCores() -1
#' @return an extended object as returned by \code{\link[Rsubread]{featureCounts}}, the extension include in the targets slots the feature extracted (exon or gene)
#' @export
#' @examples
#' \dontrun{
#' test.subject <- GetArribaRTest()
#' bam.subject <- RunSTAR(test.subject)
#' ft <- GetCounts(bam.subject)
#' }
GetCounts <- function(sbjBAM, feature = c("exon","gene"), annotSource = c("GENCODE19","inbuilt"), nthreads ){
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
  annotSource = match.arg(toupper(annotSource[1]), toupper(c("GENCODE19","inbuilt")))
  if(annotSource == "INBUILT"){
    ft <- Rsubread::featureCounts(files = sbjBAM,
                                  annot.inbuilt = "hg19",
                                  useMetaFeatures = ifelse(feature == "exon", FALSE, TRUE),
                                  isPairedEnd = TRUE,
                                  nthreads = nthreads)  
    return(invisible(ft))
  }
  # SAF <- subset(ArribaR:::.ReadExonAnnotationFile(),source == "ENSEMBL" & type == "exon")
  software <- .OpenConfigFile()
  message("Reading gen annotation")
  SAF <- rtracklayer::readGFF(file.path(software$annotation,"GENCODE19.gtf"))
  
  SAF <- subset(SAF, type == feature & gene_type %in% c("protein_coding","pseudogene","processed_transcript","IG_C_gene","antisense",
                                                        "TR_V_gene","polymorphic_pseudogene","TR_C_gene"))
  if(feature == "exon"){
    message("Getting protein coding and experimentally confirmed exon annotation")
    SAF <- subset(SAF, tag %in% c("CCDS","exp_conf"))
    SAF$GeneID <- SAF$gene_id
  }else{
    
    SAF$GeneID <- SAF$gene_id
  }
  
  colnames(SAF)[colnames(SAF) %in% c("seqid","start","end","strand")] <- c("Chr","Start","End","Strand")
  SAF$Chr <- paste0("chr",SAF$Chr)
  
  ft <- Rsubread::featureCounts(files = sbjBAM,
                                 annot.ext = SAF[,c("Chr","Start","End","Strand","GeneID")],
                                 isGTFAnnotationFile = FALSE,
                                # annot.inbuilt = "hg19",
                                useMetaFeatures = ifelse(feature == "exon", FALSE, TRUE),
                                isPairedEnd = TRUE,
                                nthreads = nthreads)
  
  ft$annotation <- cbind(ft$annotation, SAF[,which(!c(colnames(SAF) %in% colnames(ft$annotation)))])
  ft$targets <- c(ft$targets, feature)
  return(invisible(ft))
}

#' GetImmuneContent
#' 
#' Estimate the Tissue Immune Microenvironment content based on LM22 signature. It is based on TPMs calculated from the gene counts. duplicated genes will be removed
#' If MIXTURE package is not already installed, it will be installed
#' @param ftc object returned by \code{\link{GetCounts}} The ftc should be an object returned by  \code{\link{GetCounts}} called with the parameter feature = "gene"
#' this parameter is saved in ftc$targets[2] position and this function will check this, otherwise it trough an error
#' @param filter (integer) minimal count level to consider (default 0L)
#' @param plot (logical) if TRUE (default) the subject proportion plot, no plot otherwise
#' @return a MIXTURE object (see \code{\link[MIXTURE]{MIXTURE}})
#' @examples
#' \dontrun{
#' test.subject <- GetArribaRTest()
#' bam.subject <- RunSTAR(test.subject)
#' ft <- GetCounts(bam.subject)
#' ft.imm <- GetImmuneContent(ft)
#' }
GetImmuneContent <- function(ftc, filter = 0, plot = TRUE){
  software <- ArribaR:::.OpenConfigFile()
  if(c("MIXTURE" %in% rownames(installed.packages()))==FALSE){
    message("MIXTURE not installed")
    if(c("devtools" %in% rownames(installed.packages()))==FALSE){
      install.packages("devtools")
    }
    message("Installing MIXTURE")
    devtools::install_github("elmerfer/MIXTURE")
  }
  if(any(stringr::str_detect(ftc$targets,"gene"))==FALSE){
    stop("ERROR: ftc should contain gene counts")
  }

    dups <- which(duplicated(ftc$annotation$gene_name))
    counts <- ftc$counts[-dups]
    cn <- ftc$annotation$gene_name[-dups]
    effLen <- log(ftc$annotation$Length[-dups])

    # Process one column at a time.
    # tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    idx <- which(counts > filter)
    cn <- cn[idx]
    rate = log(counts[idx]) - effLen[idx]
    denom = log(sum(exp(rate)))
    tpm <- exp(rate - denom + log(1e6))
    # }))
    
    # Copy the row and column names from the original matrix.
    tpm <- as.data.frame(tpm)
    colnames(tpm) <- colnames(ftc$counts)
    rownames(tpm) <- cn
    ct <- MIXTURE(expressionMatrix = TPMs[,,drop=F], useCores = 1L)
    if(plot){
      prop <- GetMixture(ct)
      df <- data.frame(Proportions = as.vector(prop), Cells = colnames(prop))
      p <- ggplot2::ggplot(df, aes(x=Proportions, y = Cells)) + ggplot2::geom_bar(stat = "identity") +  ggplot2::scale_x_continuous(limits = c(0, 1)) + 
        ggplot2::ggtitle(paste0("Subject : ",unlist(stringr::str_split(rownames(prop),"_"))[1]))
      print(p)  
    }
    return(invisible(ct))
}

# 
# 
# 
# 
