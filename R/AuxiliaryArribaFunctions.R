# read fusions


.ArribaEnvironment <- new.env(parent = emptyenv())

# convenience functions to add/remove "chr" prefix
addChr <- function(contig) {
  ifelse(contig == "MT", "chrM", paste0("chr", contig))
}
removeChr <- function(contig) {
  sub("^chr", "", sub("^chrM", "MT", contig), perl=T)
}

# convenience function to check if a value is between two others
between <- function(value, start, end) {
  value >= start & value <= end
}

#' .ReadArribaFustionTable
#' Read and format arriba output fusion table
#' @param fusionFile string with the full path of the Fussion.tsv table
.ReadArribaFusionTable <-function(fusionsFile ){
  fusions <- read.table(fusionsFile, stringsAsFactors=F, sep="\t", header=T, comment.char="", quote="")
  if (colnames(fusions)[1] == "X.gene1") { # Arriba output
    colnames(fusions)[colnames(fusions) %in% c("X.gene1", "strand1.gene.fusion.", "strand2.gene.fusion.")] <- c("gene1", "strand1", "strand2")
    fusions$display_contig1 <- sub(":[^:]*$", "", fusions$breakpoint1, perl=T)
    fusions$display_contig2 <- sub(":[^:]*$", "", fusions$breakpoint2, perl=T)
    fusions$contig1 <- removeChr(fusions$display_contig1)
    fusions$contig2 <- removeChr(fusions$display_contig2)
    fusions$breakpoint1 <- as.numeric(sub(".*:", "", fusions$breakpoint1, perl=T))
    fusions$breakpoint2 <- as.numeric(sub(".*:", "", fusions$breakpoint2, perl=T))
    fusions$split_reads <- fusions$split_reads1 + fusions$split_reads2
  } else if (colnames(fusions)[1] == "X.FusionName") { # STAR-Fusion
    fusions$gene1 <- sub("\\^.*", "", fusions$LeftGene, perl=T)
    fusions$gene2 <- sub("\\^.*", "", fusions$RightGene, perl=T)
    fusions$strand1 <- sub(".*:(.)$", "\\1/\\1", fusions$LeftBreakpoint, perl=T)
    fusions$strand2 <- sub(".*:(.)$", "\\1/\\1", fusions$RightBreakpoint, perl=T)
    fusions$display_contig1 <- sub(":[^:]*:[^:]*$", "", fusions$LeftBreakpoint, perl=T)
    fusions$display_contig2 <- sub(":[^:]*:[^:]*$", "", fusions$RightBreakpoint, perl=T)
    fusions$contig1 <- removeChr(fusions$display_contig1)
    fusions$contig2 <- removeChr(fusions$display_contig2)
    fusions$breakpoint1 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", fusions$LeftBreakpoint, perl=T))
    fusions$breakpoint2 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", fusions$RightBreakpoint, perl=T))
    fusions$direction1 <- ifelse(grepl(":\\+$", fusions$LeftBreakpoint, perl=T), "downstream", "upstream")
    fusions$direction2 <- ifelse(grepl(":\\+$", fusions$RightBreakpoint, perl=T), "upstream", "downstream")
    fusions$transcript_id1 <- ifelse(rep(!("CDS_LEFT_ID" %in% colnames(fusions)), nrow(fusions)), ".", fusions$CDS_LEFT_ID)
    fusions$transcript_id2 <- ifelse(rep(!("CDS_RIGHT_ID" %in% colnames(fusions)), nrow(fusions)), ".", fusions$CDS_RIGHT_ID)
    fusions$fusion_transcript <- ifelse(rep(!("FUSION_CDS" %in% colnames(fusions)), nrow(fusions)), ".", toupper(sub("([a-z]*)", "\\1|", fusions$FUSION_CDS)))
    fusions$reading_frame <- ifelse(rep(!("PROT_FUSION_TYPE" %in% colnames(fusions)), nrow(fusions)), ".", ifelse(fusions$PROT_FUSION_TYPE == "INFRAME", "in-frame", ifelse(fusions$PROT_FUSION_TYPE == "FRAMESHIFT", "out-of-frame", ".")))
    fusions$split_reads <- fusions$JunctionReadCount
    fusions$discordant_mates <- fusions$SpanningFragCount
    fusions$site1 <- rep("exon", nrow(fusions))
    fusions$site2 <- rep("exon", nrow(fusions))
    fusions$confidence <- rep("high", nrow(fusions))
  } else {
    stop("Unrecognized fusion file format")
  }

return(invisible(fusions))
}

#' .ReadCytobandFile
#' Internal function. Read cytoband information
#' @param fusionFile string with the full path of the Fussion.tsv table
.ReadCytobandFile <- function(){
  software <- .OpenConfigFile()
  ##-----------
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA

  ##--------
  ##------------
  if(!exists("cytobands", envir = .ArribaEnvironment, inherits = FALSE)){
    message("Loading cytobands")
    cytobandsFile <- list.files(software$arriba$database, full.names = TRUE)
    cytobandsFile <- cytobandsFile[stringr::str_detect(cytobandsFile, "cytobands")]
    cytobandsFile <- cytobandsFile[stringr::str_detect(cytobandsFile, software$arriba$assemblyVersion)]
    cytobands <- read.table(cytobandsFile, header=T, colClasses=c("character", "numeric", "numeric", "character", "character"))
    cytobands <- cytobands[order(cytobands$contig, cytobands$start, cytobands$end),]
    .ArribaEnvironment$cytobands <- cytobands
  }else{
    message("cytobands already loaded")
    cytobands <- .ArribaEnvironment$cytobands
  }
  return(cytobands)
}



.ReadExonAnnotationFile <- function(fusionsTable ){
  software <- .OpenConfigFile()
  ##-----------
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA

  ##--------
  if(!exists("exons", envir = .ArribaEnvironment, inherits = FALSE)){
    message("Loading exons")
    exonsFile <- file.path(software$annotation,"GENCODE19.gtf")
    exons <- read.table(exonsFile, header=F, sep="\t", comment.char="#", quote="", stringsAsFactors=F)[,c(1,2, 3, 4, 5, 7, 9)]
    # exons <- data.table::fread(paste0("grep -v '^#' ",exonsFile),
    #                            nThread = 4,
    #                            header=F, sep="\t", quote="", stringsAsFactors=F)[,c(1, 3, 4, 5, 7, 9)]
    colnames(exons) <- c("contig", "source" ,"type", "start", "end", "strand", "attributes")
    exons <- exons[exons$type %in% c("exon", "CDS"),]
    exons$contig <- removeChr(exons$contig)
    
    parseGtfAttribute <- function(attribute, exons) {
      parsed <- gsub(paste0(".*", attribute, " \"?([^;\"]+)\"?;.*"), "\\1", exons$attributes)
      failedToParse <- parsed == exons$attributes
      if (any(failedToParse)) {
        warning(paste0("Failed to parse '", attribute, "' attribute of ", sum(failedToParse), " GTF record(s)."))
        parsed <- ifelse(failedToParse, "", parsed)
      }
      return(parsed)
    }
    exons$geneID <- parseGtfAttribute("gene_id", exons)
    exons$geneName <- parseGtfAttribute("gene_name", exons)
    exons$geneName <- ifelse(exons$geneName == "", exons$geneID, exons$geneName)
    exons$transcript <- parseGtfAttribute("transcript_id", exons)
    exons$exonNumber <- ifelse(rep(printExonLabels, nrow(exons)), parseGtfAttribute("exon_number", exons), "")    
    
    .ArribaEnvironment$exons <- exons
  }else{
    message("Using preloaded exons")
    exons <- .ArribaEnvironment$exons[,]
  }
  # exonsFile <- file.path(software$annotation,"GENCODE19.gtf")
  # exons <- read.table(exonsFile, header=F, sep="\t", comment.char="#", quote="", stringsAsFactors=F)[,c(1, 3, 4, 5, 7, 9)]
  # # exons <- data.table::fread(paste0("grep -v '^#' ",exonsFile),
  # #                            nThread = 4,
  # #                            header=F, sep="\t", quote="", stringsAsFactors=F)[,c(1, 3, 4, 5, 7, 9)]
  # colnames(exons) <- c("contig", "type", "start", "end", "strand", "attributes")
  # exons <- exons[exons$type %in% c("exon", "CDS"),]
  # exons$contig <- removeChr(exons$contig)
  # 
  # parseGtfAttribute <- function(attribute, exons) {
  #   parsed <- gsub(paste0(".*", attribute, " \"?([^;\"]+)\"?;.*"), "\\1", exons$attributes)
  #   failedToParse <- parsed == exons$attributes
  #   if (any(failedToParse)) {
  #     warning(paste0("Failed to parse '", attribute, "' attribute of ", sum(failedToParse), " GTF record(s)."))
  #     parsed <- ifelse(failedToParse, "", parsed)
  #   }
  #   return(parsed)
  # }
  # exons$geneID <- parseGtfAttribute("gene_id", exons)
  # exons$geneName <- parseGtfAttribute("gene_name", exons)
  # exons$geneName <- ifelse(exons$geneName == "", exons$geneID, exons$geneName)
  # exons$transcript <- parseGtfAttribute("transcript_id", exons)
  # exons$exonNumber <- ifelse(rep(printExonLabels, nrow(exons)), parseGtfAttribute("exon_number", exons), "")
  if(missing(fusionsTable)){
    return(exons)
  }

  # insert dummy annotations for intergenic breakpoints
  if (any(fusionsTable$site1 == "intergenic" | fusionsTable$site2 == "intergenic")) {
    intergenicBreakpoints <- rbind(
      setNames(fusionsTable[fusionsTable$site1 == "intergenic",c("gene1", "strand1", "contig1", "breakpoint1")], c("gene", "strand", "contig", "breakpoint")),
      setNames(fusionsTable[fusionsTable$site2 == "intergenic",c("gene2", "strand2", "contig2", "breakpoint2")], c("gene", "strand", "contig", "breakpoint"))
    )
    exons <- rbind(exons, data.frame(
      contig=intergenicBreakpoints$contig,
      type="intergenic",
      start=intergenicBreakpoints$breakpoint-1000,
      end=intergenicBreakpoints$breakpoint+1000,
      strand=".",
      attributes="",
      geneName=intergenicBreakpoints$gene,
      geneID=intergenicBreakpoints$gene,
      transcript=intergenicBreakpoints$gene,
      exonNumber="intergenic", 
      source = ""
    ))
  }
 return(exons)
}


.ReadProteinDomainsFile <- function(){
  software <- .OpenConfigFile()
  ##-----------
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA
  ##--------
  if(!exists("proteinDomains", envir = .ArribaEnvironment, inherits = FALSE)){
    proteinDomainsFile <- list.files(software$arriba$database, full.names = TRUE)
    proteinDomainsFile <- proteinDomainsFile[stringr::str_detect(proteinDomainsFile,"protein")]
    proteinDomainsFile <- proteinDomainsFile[stringr::str_detect(proteinDomainsFile,software$arriba$assemblyVersion)]

    message("Loading protein domains from file")
    proteinDomains <- read.table(proteinDomainsFile, header=F, sep="\t", comment.char="", quote="", stringsAsFactors=F)[,c(1,4,5,7,9)]
    colnames(proteinDomains) <- c("contig", "start", "end", "strand", "attributes")
    proteinDomains$color <- sub(";.*", "", sub(".*color=", "", proteinDomains$attributes, perl=T), perl=T)
    proteinDomains$proteinDomainName <- sapply(sub(";.*", "", sub(".*Name=", "", proteinDomains$attributes, perl=T), perl=T), URLdecode)
    proteinDomains$proteinDomainID <- sub(";.*", "", sub(".*protein_domain_id=", "", proteinDomains$attributes, perl=T), perl=T)
    proteinDomains <- proteinDomains[,colnames(proteinDomains) != "attributes"]
    .ArribaEnvironment$proteinDomains <- proteinDomains
    
  }else{
    message("protein domains already loaded")
    proteinDomains <- .ArribaEnvironment$proteinDomains
  }
  return(proteinDomains)
}





