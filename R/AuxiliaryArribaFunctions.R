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
  if(stringr::str_detect(fusionsFile,".xlsx")){
    fusions <- openxlsx::read.xlsx(fusionsFile)  
    first.pattern <- "gene1"
    colnames(fusions) <- make.names(colnames(fusions))
    attr(fusions,"AssemblyInfo") <- openxlsx::read.xlsx(fusionsFile,sheet = "Info")
  }else{
    fusions <- read.table(fusionsFile, stringsAsFactors=F, sep="\t", header=T, comment.char="", quote="")
    first.pattern <- "X.gene1"
    attr(fusions,"AssemblyInfo") <- "Unknown"
  }
  # fusions <- read.table(fusionsFile, stringsAsFactors=F, sep="\t", header=T, comment.char="", quote="")
  if (colnames(fusions)[1] == first.pattern) { # Arriba output
    colnames(fusions)[colnames(fusions) %in% c(first.pattern, "strand1.gene.fusion.", "strand2.gene.fusion.")] <- c("gene1", "strand1", "strand2")

    fusions$display_contig1 <- sub(":[^:]*$", "", fusions$breakpoint1, perl=T)
    fusions$display_contig2 <- sub(":[^:]*$", "", fusions$breakpoint2, perl=T)
    fusions$contig1 <- removeChr(fusions$display_contig1)
    fusions$contig2 <- removeChr(fusions$display_contig2)
    fusions$breakpoint1 <- as.numeric(sub(".*:", "", fusions$breakpoint1, perl=T))
    fusions$breakpoint2 <- as.numeric(sub(".*:", "", fusions$breakpoint2, perl=T))
    fusions$split_reads <- fusions$split_reads1 + fusions$split_reads2
    
    # fusions$contig1 <- removeChr(fusions$display_contig1)
    # fusions$contig2 <- removeChr(fusions$display_contig2)
  } else if (colnames(fusions)[1] == "X.FusionName") { # STAR-Fusion
    fusions$gene1 <- sub("\\^.*", "", fusions$LeftGene, perl=T)
    fusions$gene2 <- sub("\\^.*", "", fusions$RightGene, perl=T)
    # fusions$contig1 <- removeChr(fusions$display_contig1)
    # fusions$contig2 <- removeChr(fusions$display_contig2)
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
#' @param species see \code{\link[GenomeDB]{GetGenomes}}
#' @param version see \code{\link[GenomeDB]{GetGenomes}}
.ReadCytobandFile <- function(species,version){
  software <- GenomeDB:::.OpenConfigFile()
  ##-----------
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA

  ##--------
  ##------------
  if(!exists("cytobands", envir = .ArribaEnvironment, inherits = FALSE)){
    message("Loading cytobands")
    cytobandsFile <- list.files(software$Software$ARRIBA$database, full.names = TRUE)
    cytobandsFile <- cytobandsFile[stringr::str_detect(cytobandsFile, "cytobands")]
    cytobandsFile <- cytobandsFile[stringr::str_detect(cytobandsFile, "GRCh38")]
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
  software <- GenomeDB:::.OpenConfigFile()
  
  if(!exists("exons", envir = .ArribaEnvironment, inherits = FALSE)){
    message("Loading exons")
    assembly <- attr(fusionsTable,"AssemblyInfo")
    fasta.gtf.files <- GenomeDB::GetGenome("Human","GRCh38+GENECODE")
    # exonsFile <- software$Software$STAR$main[,assembly["ASSEMBLY"]]
    exons <- rtracklayer::readGFF(fasta.gtf.files$gtf)[, c("seqid","source","type","start","end","strand","gene_id","gene_name","transcript_id","exon_number")]
    
    colnames(exons) <- c("contig", "source" ,"type", "start", "end", "strand", "geneID","geneName","transcript","exonNumber")
    exons <- exons[exons$type %in% c("exon", "CDS"),]
    exons$contig <- removeChr(exons$contig)
    
    .ArribaEnvironment$exons <- exons
  }else{
    message("Using preloaded exons")
    exons <- .ArribaEnvironment$exons[,]
  }
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
      source = "",
      type="intergenic",
      start=intergenicBreakpoints$breakpoint-1000,
      end=intergenicBreakpoints$breakpoint+1000,
      strand=".",
      geneID=intergenicBreakpoints$gene,
      geneName=intergenicBreakpoints$gene,
      transcript=intergenicBreakpoints$gene,
      exonNumber="intergenic"
    ))
  }
  return(exons)
}


.ReadExonAnnotationFile2 <- function(fusionsTable ){
  software <- GenomeDB:::.OpenConfigFile()
  
  ##-----------
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA

  ##--------
  if(!exists("exons", envir = .ArribaEnvironment, inherits = FALSE)){
    message("Loading exons")
    fasta.gtf.files <- GenomeDB::GetGenome("Human","GRCh38+GENECODE")##OJO con esto
    # exonsFile <- file.path(software$annotation,"GENCODE19.gtf")
    exons <- read.table(fasta.gtf.files$gtf, header=F, sep="\t", comment.char="#", quote="", stringsAsFactors=F)[,c(1,2, 3, 4, 5, 7, 9)]
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
  software <- GenomeDB:::.OpenConfigFile()
  ##-----------
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA
  ##--------
  if(!exists("proteinDomains", envir = .ArribaEnvironment, inherits = FALSE)){
    proteinDomainsFile <- list.files(software$Software$ARRIBA$database, full.names = TRUE)
    proteinDomainsFile <- proteinDomainsFile[stringr::str_detect(proteinDomainsFile,"protein")]
    proteinDomainsFile <- proteinDomainsFile[stringr::str_detect(proteinDomainsFile,"GRCh38")] ##OJO con esto

    message("Loading protein domains from file versio GRCh38")
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


#' .StatsSTARfile
#' reads and format the Log-final.out alignement statistics file
#' This is an internal function, it should not be directly used
#' @param bamFile the prefix_Aligned.out.bam or prefix_SortByCoordinates.bam file
#' @return a data frame with the Log content
.StatsSTARfile <- function(bamFile){
  log.file <- list.files(dirname(bamFile), full.names = T)
  log.file <- log.file[stringr::str_detect(log.file, "Log.final.out")]
  ##get the statistics
  sd <- plyr::ldply(readLines(log.file), function(x) {
    if(stringr::str_length(x)<2) return(c(NA,NA))
    ret <- unlist(stringr::str_split(x,"\t"))
    if(length(ret)<2) return(c(ret,NA))
    ret
  })
  colnames(sd) <- c("Description","Value")
  return(sd)
}


.panel.par <- function(gene1,gene2, fileFusions) {
  cat(paste0("\nNow searching pair\t",gene1,"\t",gene2,"\t:\t"))
  ss<-system2( command="awk", args= sprintf('\'$1==\"%s" && $2==\"%s" \'  %s',gene1,gene2,fileFusions),
               stdout = TRUE)
  
  if(length(ss)>0){
    ret <- plyr::ldply(ss, function(x){
      unlist(stringr::str_split(x,"\t"))
    })  
  }else{
    ret<-NA
  }
  cat(paste0(ifelse(is.null(nrow(ret)),0,nrow(ret))," gene pairs found"))
  ret
  
}

.panel.gene1 <- function(gene1, fileFusions){
  cat(paste0("\nNow searching 5' gene\t",gene1,"\t:\t"))
  ss <- system2( command="awk", args= sprintf(' \'$1==\"%s" \' %s',gene1,fileFusions),
                 stdout = TRUE)
  if(length(ss)>0){
    ret <- plyr::ldply(ss, function(x){
      unlist(stringr::str_split(x,"\t"))
    })  
  }else{
    ret <- NA
  }
  
  cat(paste0(ifelse(is.null(nrow(ret)),0,nrow(ret))," - 5' genes found"))
  ret
}

.panel.gene2 <- function(gene2, fileFusions){
  cat(paste0("\nNow searching 3' gene\t",gene2,"\t:\t"))
  ss <- system2( command="awk", args= sprintf(' \'$2==\"%s" \' %s',gene2,fileFusions),
                 stdout = TRUE)
  if(length(ss)>0){
    ret <- plyr::ldply(ss, function(x){
      unlist(stringr::str_split(x,"\t"))
    })  
  }else{
    ret <- NA
  }
  
  cat(paste0(ifelse(is.null(nrow(ret)),0,nrow(ret))," - 3' genes found"))
  ret
}

.Mitelman <- function(fileFusions){
  cat(paste0("\nNow searching Mitelman tagged fussions\t:\t"))
  ss <- system2( command="awk", args= sprintf(' \'$17=="Mitelman" \' %s',fileFusions),
                 stdout = TRUE)
  if(length(ss)>0){
    ret <- plyr::ldply(ss, function(x){
      unlist(stringr::str_split(x,"\t"))
    })  
  }else{
    ret <- NA
  }
  
  cat(paste0(ifelse(is.null(nrow(ret)),0,nrow(ret))," fusions found"))
  ret
}


#' SearchPanel
#' This function recover those fusions that were discarded by Arriba, 
#' but they have at split_read1+spleat_read2+discordant_mates > 0. Low confident fusions
#'@param genePanel a data frame with the target genes
#'@param fileFusions the tsv discarded file fusions path, it should end in _Fusions_discarded.tsv
#'@param save (logical, default TRUE), if the results would be saved as excel file
#'@export
SearchPanel <- function(genePanel, fileFusions, save=T){
  cnames <- stringr::str_remove_all(colnames(data.table::fread(file = fileFusions,nrows=1)),"#")
  pares <- na.omit(genePanel[,c("gene1","gene2"),drop=F])
  gene1 <- genePanel[is.na(genePanel$gene2),c("gene1","gene2"),drop=F]
  gene2 <- genePanel[is.na(genePanel$gene1),c("gene1","gene2"),drop=F]
  if(nrow(pares)>0){
    ret.p<- plyr::ldply(1:nrow(pares), function(px){
      aux <- .panel.par(pares$gene1[px],pares$gene2[px], fileFusions = fileFusions)
      if(all(is.na(aux[[1]]))) return(NULL)
      aux
    })
    colnames(ret.p) <- cnames
  }else{
    ret.p<-NA
  }
  if(nrow(gene1)>0){
    ret.g1 <- plyr::ldply(1:nrow(gene1), function(px){
      aux <- .panel.gene1(pares$gene1[px], fileFusions = fileFusions)
      if(all(is.na(aux[[1]]))) return(NULL)
      aux
    })
    colnames(ret.g1) <- cnames
  }else{
    ret.g1<-NA
  } 
  if(nrow(gene2)>0){
    ret.g2 <- plyr::ldply(1:nrow(gene2), function(px){
      aux<-.panel.gene2(pares$gene2[px], fileFusions = fileFusions)
      if(all(is.na(aux[[1]]))) return(NULL)
      aux
    })
    colnames(ret.g2) <- cnames
  }else{
    ret.g2<-NA
  }
  
  ret.m <- .Mitelman(fileFusions = fileFusions)
  colnames(ret.m) <- cnames
  pan <- list(Pairs=ret.p,"gene1-any"=ret.g1,"any-gene2"=ret.g2, Mitelman = ret.m)
  
  file.name <- stringr::str_replace(fileFusions,"_discarded.tsv","_Recovered.xlsx")
  if(save) {
    openxlsx::write.xlsx(pan, file = file.name, overwrite = F)
    if(file.exists(file.name)){
      cat(paste0("\nFile ", file.name, " Saved"))
    }else{
      cat("ERROR")
    }
  }
  return(invisible(ret))
}


