#' RunSTARforARRIBA
#' This function will run STAR software
#' based on GenomeDB
#' @param sbjFile string with the full path and name of the xxxx_1_fastq or gz sequence file. This
#' function only support paired data
#' @param nThreads number of CPUs threads 
#' @export
#' @author Elmer A. Fernández
#' @details
#' In order to look for gene fusions with the ARRIBA software, the STAR aligner should be run with special parameters.
#' The RunSTARforARRIBA function performe such task.
#' @examples
#' \dontrun{
#' test.subject <- GetArribaRTest()
#' bam.subject <- RunSTAR2(test.subject)
#' ##to display the genome and annotation version
#' attr(out.file,"GenomeDB")
#' ##to display the elapsed time of the STAR run
#' attr(out.file,"ElapsedTime") 
#' }
RunSTARforARRIBA <- function(sbjFile, nThreads, version = "GRCh38+GENECODE"){
  out.file <- Aligners::RunSTAR(sbjFile = sbjFile, nThreads = nThreads , species = "Human",version=version)
  df.stat <- .StatsSTARfile(out.file)
  cat(paste0("Aligned STAR file : ",out.file))
  attr(out.file,"GenomeDB") <- version
  attr(out.file,"STAR_STATS") <- df.stat
  return(out.file)
}


  
#' RunARRIBA
#' this function will run the gene fusion detection ARRIBA software
#' @references  \url{https://genome.cshlp.org/content/early/2021/02/11/gr.257246.119}{Uhrig et al.}
#' @param sbjBamFile string with the full file name of the subject BAM file
#' @param allProteinPredictions if TRUE it will generate all the prediction of the discarded file. It may generate a very big discarded file. it is similar to set -X in the Arriba software
#' @seealso \code{\link{runSTAR()}}
#' @export
#' @return a data frame with the identified fusions
#' \url{https://arriba.readthedocs.io/en/latest/output-files/}{see details in Arriba output files}
#' @examples
#' \dontrun{
#' test.subject <- GetArribaRTest()
#' bam.subject <- RunSTAR(test.subject)
#' Fusions <- RunArriba(bam.subject)
#' View(Fusions)
#' }
RunARRIBA <- function(sbjBamFile, allProteinPredictions = FALSE){
  assemblyVersion <- attr(sbjBamFile,"GenomeDB")
  if(is.null(assemblyVersion)){
    stop("\nPlease process the subject by RunSTAR")
  }
  
  software <- GenomeDB:::.OpenConfigFile()
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA
  
  ##--------
  # call arriba
  if(!file.exists(sbjBamFile)){
    stop(paste0("\nERROR ", basename(sbjBamFile), "NOT FOUND"))
  }
  
  fusion.file <- stringr::str_replace(sbjBamFile,paste0(software$Software$STAR$alignmentPrefix,"Aligned.out"),"_Fusions.tsv")
  fusion.file <- stringr::str_remove_all(fusion.file,".bam")
  fusion.discarded.file <- stringr::str_replace(sbjBamFile,paste0(software$Software$STAR$alignmentPrefix,"Aligned.out"),"_Fusions_discarded.tsv")
  fusion.discarded.file <- stringr::str_remove_all(fusion.discarded.file,".bam")
  
  # genomeversion <- ifelse(assemblyVersion == "GRCh37", "hg19_hs37d5_","hg38_")
  Genomefiles <- list.files(software$GenomesDB$Human$version[assemblyVersion],full.names = T)
  GenomeAssembly <- Genomefiles[stringr::str_detect(Genomefiles,".fa")]
  GenomeGTF <- Genomefiles[stringr::str_detect(Genomefiles,".gtf")]
  if(!any(file.exists(c(GenomeAssembly,GenomeGTF)))){
    stop("Genome not found")
  }
  genomeversion <- ifelse(stringr::str_detect(GenomeAssembly,"GRCh38"),"hg38_GRCh38","hg19hs37d5_GRCh37")
  system2(command = software$Software$ARRIBA$command,
          args = c(paste0("-x ",sbjBamFile),
                   paste0("-o ",fusion.file),
                   paste0("-O ",fusion.discarded.file),
                   paste0("-a ",GenomeAssembly),
                   paste0("-g ", GenomeGTF),
                   paste0("-b ", file.path(software$Software$ARRIBA$database,
                                           paste0("blacklist_",genomeversion,"_v2.1.0.tsv.gz"))),
                   paste0("-k ", file.path(software$Software$ARRIBA$database,
                                           paste0("known_fusions_",genomeversion,"_v2.1.0.tsv.gz"))),
                   paste0("-t ", file.path(software$Software$ARRIBA$database,
                                           paste0("known_fusions_",genomeversion,"_v2.1.0.tsv.gz"))),
                   paste0("-p ", file.path(software$Software$ARRIBA$database,
                                           paste0("protein_domains_",genomeversion,"_v2.1.0.gff3"))),
                   ifelse(allProteinPredictions,"-X","")
          ))
  
  if(!all(file.exists(fusion.file,fusion.discarded.file))){
    cat(paste0("\nFUSION FILES NOT FOUND\n",fusion.file,"\n", fusion.discarded.file))
  }
  ##write as an excel file
  fusionsTable <- .ReadArribaFusionTable(fusion.file)
  version.info <- data.frame(  STAR = software$Software$STAR$version, 
                               ARRIBA = software$Software$ARRIBA$version, 
                               ASSEMBLY =assemblyVersion,
                               ArribaR=paste0(unlist(packageVersion("ArribaR")),collapse = "." ),
                               SAMPLE = stringr::str_remove(basename(fusion.file),"_Fusions.tsv"))
  sts <- na.omit(attr(sbjBamFile,"STAR_STATS"))
  openxlsx::write.xlsx(list(Fusions=fusionsTable,Info=version.info, 
                       stats=sts), stringr::str_replace(fusion.file,".tsv",".xlsx"),overwrite = T)
  system2( command="awk", args= sprintf(' \'$10!=0 && $11!=0 && $12 !=0\' %s',fusion.discarded.file),
           stdout = stringr::str_replace(fusion.discarded.file,".tsv","_filtered.tsv" ))
  if(file.exists(stringr::str_replace(fusion.discarded.file,".tsv","_filtered.tsv" ))){
    file.remove(fusion.discarded.file)
  }
  file.remove(fusion.file)
  #return(fusionsTable) ##for continuous processing into R
  # if(!stringr::str_detect(sbjBamFile,software$star$alignmentPrefixSorted)){
  #   Rsamtools::sortBam(file = sbjBamFile,
  #                      destination = stringr::str_replace(sbjBamFile, software$star$aligmentPrefix,software$star$alignmentPrefixSorted))
  #   test.out <- stringr::str_replace(sbjBamFile, software$star$aligmentPrefix,software$star$alignmentPrefixSorted)
  #   Rsamtools::indexBam(test.out)
  # }
  return(invisible(fusionsTable))
}

#' #' RunSortIndexBam
#' #' This function will built the sorted bam file and its index which are necesary for visualization
#' #' @param sbjBamFile string full path of aligned bam file see \code{\link{RunSTARforARRIBA}}
#' #' @param remove boolean (default  TRUE), if unsorted bam file should be removed or not.
#' #' @export
#' #' @seealso \code{\link[Rsamtools]{sortBam}}
RunSortIndexBam <- function(sbjBamFile, remove = T){
  software <- .OpenConfigFile()
  if(!file.exists(sbjBamFile)){
    stop(paste0("\nERROR ",sbjBamFile," NOT FOUND"))
  }
  if(!stringr::str_detect(sbjBamFile,".bam")){
    stop(paste0("\nERROR this is not a bam file",sbjBamFile,"\n"))
  }
  if(stringr::str_detect(sbjBamFile, software$star$alignmentPrefixSorted)){
    stop(paste0("\nThe file seems to be already sorted, pls verify"))
  }
  hd <- ArribaR:::.ReadBamHeader(sbjBamFile)
  if(hd$ProgramName!="STAR"){
    stop("It should be aligned by RunSTARforARRIBA function")
  }else{
    destination.file <- stringr::str_replace(sbjBamFile, ".bam", "SortedByCoordinates")
  }
  Rsamtools::sortBam(file = sbjBamFile,
                     destination = destination.file)
  destination.file <- paste0(destination.file,".bam")
  if(file.exists(destination.file)){
    Rsamtools::indexBam(destination.file)
    if(remove) file.remove(sbjBamFile)
  }

  return(destination.file)

}

##' .ReadBamHeader
##' @param bamFile bam file full path
##' @return a list with the following slots:
##' CHR: the list of chromosomes and sequence names in the bam file
##' ProgramName : the alignment program (PG:PN see [BAM header](https://samtools.github.io/hts-specs/SAMv1.pdf))
##' ProgramVersion : the version of the aligner
##' Code : executed source line code
##' GenomeDBversion : the genome version used (see GenomeDB)
##' GenomeDBpath : the path to genome version
##' 
.ReadBamHeader <- function(bamFile){
  hd <- Rsamtools::scanBamHeader(bamFile)
  ##program name and version
  pname <- stringr::str_remove_all(grep(pattern = "PN", unlist(hd),value=TRUE),"PN:")
  pversion <- stringr::str_remove_all(grep(pattern = "VN", unlist(hd),value=TRUE),"VN:")
  bam.chrs <- names(hd[1]$targets)
  CL <- grep(pattern = "CL", unlist(hd),value=TRUE)
  code <- unlist(str_split(CL," "))
  id<- which(str_detect(code,"--genomeDir"))
  
  gpath <- code[id+1]
  genome <- basename(gpath)
  id<- which(str_detect(code,"--sjdbGTFfile"))
  gtf <- code[id+1]
  return(list(CHR=bam.chrs,ProgramName = pname, ProgramVersion=pversion, Code =CL, GenomeDBversion=genome, GenomeDBpath = gpath, GTF=gtf))
}

#' GetArribaRTest
#' This function will retrieve the test files coming with the arriba software
#' @usage test.files <- GetArribaRTest()
#' @return an string with the full path of the firs read file of the test
GetArribaRTest <- function(){
  software <- GenomeDB:::.OpenConfigFile()
  files <- list.files(paste0(software$Software$ARRIBA$main,"/test"),full.names = T)
  id.fastq <- which(stringr::str_detect(files,"read1.fastq.gz"))
  id.fastq2 <- which(stringr::str_detect(files,"read2.fastq.gz"))
  if(length(id.fastq)>0){
    file.rename(files[id.fastq],stringr::str_replace(files[id.fastq],"read1.fastq.gz","read_1.fastq.gz"))
    file.rename(files[id.fastq2],stringr::str_replace(files[id.fastq2],"read2.fastq.gz","read_2.fastq.gz"))
  }else{
    if(file.exists(paste0(software$Software$ARRIBA$main,"/test/read_1.fastq.gz"))){
      return(paste0(software$Software$ARRIBA$main,"/test/read_1.fastq.gz"))
    }
    return(NULL)
  }
  if(file.exists(paste0(software$Software$ARRIBA$main,"/test/read_1.fastq.gz"))){
    return(paste0(software$Software$ARRIBA$main,"/test/read_1.fastq.gz"))
  }
  return(NULL)
}
