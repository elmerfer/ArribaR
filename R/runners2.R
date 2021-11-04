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
  software <- GenomeDB:::.OpenConfigFile()
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA
  cat(paste0("\nCurrent STAR version : ", software$star$version, "\n"))
  # assemblyVersion <- paste0("GRCh",match.arg(as.character(assemblyVersion), choices = c("37","38")))
  cat(paste0("\nChosen assembly version GRCh : ", version, "\n"))
  ##--------
  ##hay que revisar si es GZ
  file1 <- sbjFile
  if(stringr::str_detect(file1, "_R1.fastq")){
    file2 <- stringr::str_replace(file1,"_R1.fastq","_R2.fastq")
  }else{
    if(stringr::str_detect(file1, "_1.fastq")){
      file2 <- stringr::str_replace(file1,"_1.fastq","_2.fastq")
    }else{
      stop("unrecognizable file, it should be sbj_R1.fastq/.gz or sbj_1.fastq/.gz")
    }
  }
  
  cat(paste0("\nSubject files\n",file1,"\n",file2,"\n"))
  if(!all(file.exists(file1,file2))){
    stop("ERROR some files not found")
  }
  out.file <- file1
  out.file <- stringr::str_remove(file1,".gz")
  out.file <- stringr::str_replace(out.file,"1.fastq",software$Software$STAR$alignmentPrefix)
  out.file <- paste0(out.file,".bam")
  
  
  nThreads <- .SetThreads(nThreads)
  
  if(nThreads < 4){
    message(paste0("STAR running on ",nThreads, " may be not optimal"))
  }else{
    message(paste0("STAR running on ",nThreads, " threads"))
  }
  genomeDir <- software$Software$STAR$main[version] #paste0(software$assembly,"/STAR_",assemblyVersion)
  # genomeDir <- ifelse(assemblyVersion == "GRCh37",paste0(genomeDir,"_GENCODE19_index"),
  #                     paste0(genomeDir,"_GENCODE38_index"))
  t1 <- Sys.time()
  system2(command = software$Software$STAR$command,
          args = c(paste0("--runThreadN ",nThreads),
                   paste0("--genomeDir " ,genomeDir),
                   paste0("--readFilesIn ",file1, " ", file2),
                   paste0("--readFilesCommand ", ifelse(stringr::str_detect(file1,".gz"),"zcat","-")),
                   "--outStd BAM_Unsorted",
                   "--outSAMtype BAM Unsorted",
                   "--outSAMunmapped Within",
                   "--outBAMcompression 0",
                   "--outFilterMultimapNmax 50",
                   "--peOverlapNbasesMin 10",
                   "--alignSplicedMateMapLminOverLmate 0.5",
                   "--alignSJstitchMismatchNmax 5 -1 5 5",
                   "--chimSegmentMin 10",
                   "--chimOutType WithinBAM HardClip",
                   "--chimJunctionOverhangMin 10",
                   "--chimScoreDropMax 30",
                   "--chimScoreJunctionNonGTAG 0",
                   "--chimScoreSeparation 1",
                   "--chimSegmentReadGapMax 3",
                   "--chimMultimapNmax 50"
          ), stdout = out.file)
  if(!file.exists(out.file)){
    stop("ERROR aligment")
  }
  cat(paste0("Aligned STAR file : ",out.file))
  attr(out.file,"GenomeDB") <- version
  attr(out.file,"ElapsedTime") <- round(Sys.time() - t1,3)
  return(out.file)
}


#' RunARRIBA
#' this function will run the gene fusion detection ARRIBA software
#' @references  \url{https://genome.cshlp.org/content/early/2021/02/11/gr.257246.119}{Uhrig et al.}
#' @param sbjBamFile string with the full file name of the subject BAM file
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
RunARRIBA <- function(sbjBamFile){
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
  
  fusion.file <- stringr::str_replace(sbjBamFile,software$Software$STAR$alignmentPrefix,"_Fusions.tsv")
  fusion.file <- stringr::str_remove_all(fusion.file,".bam")
  fusion.discarded.file <- stringr::str_replace(sbjBamFile,software$Software$STAR$alignmentPrefix,"_Fusions_discarded.tsv")
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
                                           paste0("protein_domains_",genomeversion,"_v2.1.0.gff3")))
          ))
  
  if(!all(file.exists(fusion.file,fusion.discarded.file))){
    cat(paste0("\nFUSION FILES NOT FOUND\n",fusion.file,"\n", fusion.discarded.file))
  }
  ##write as an excel file
  fusionsTable <- .ReadArribaFusionTable(fusion.file)
  version.info <- data.frame(  STAR = software$Software$STAR$version, 
                               ARRIBA = software$Software$ARRIBA$version, 
                               ASSEMBLY =assemblyVersion,
                               ArribaR=packageVersion("ArribaR"),
                               SAMPLE = stringr::str_remove(basename(fusion.file),"_Fusions.tsv"))
  openxlsx::write.xlsx(list(Fusions=fusionsTable,Info=version.info), stringr::str_replace(fusion.file,".tsv",".xlsx"))
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

#' #' runSortIndexBam
#' #' This function will built the sorted bam file and its index.
#' #' @param sbjBamFile string full path of aligned bam file see \code{\link[runSTAR]{runSTAR}}
#' #' @param remove boolean (default  TRUE), if unsorted bam file should be removed or not.
#' #' @export
#' #' @seealso \code{\link[Rsamtools]{sortBam}}
#' RunSortIndexBam <- function(sbjBamFile, remove = T){
#'   software <- .OpenConfigFile()
#'   if(!file.exists(sbjBamFile)){
#'     stop(paste0("\nERROR ",sbjBamFile," NOT FOUND"))
#'   }
#'   if(!stringr::str_detect(sbjBamFile,".bam")){
#'     stop(paste0("\nERROR this is not a bam file",sbjBamFile,"\n"))
#'   }
#'   if(stringr::str_detect(sbjBamFile, software$star$alignmentPrefixSorted)){
#'     stop(paste0("\nThe file seems to be already soprted, pls verify"))
#'   }
#'   if(!stringr::str_detect(sbjBamFile,software$star$alignmentPrefix)){
#'     warning("The file seems not being aligned trhough ArribaR::runSTAR() function")
#'     destination.file <- stringr::remove(
#'       paste0(sbjBamFile, software$star$alignmentPrefixSorted), ".bam")
#'   }else{
#'     destination.file <- stringr::str_remove(
#'       stringr::str_replace(sbjBamFile, software$star$alignmentPrefix, software$star$alignmentPrefixSorted), ".bam")
#'   }
#'   Rsamtools::sortBam(file = sbjBamFile,
#'                      destination = destination.file)
#'   destination.file <- paste0(destination.file,".bam")
#'   if(file.exists(destination.file)){
#'     Rsamtools::indexBam(destination.file)
#'     if(remove) file.remove(sbjBamFile)
#'   }
#'   
#'   return(destination.file)
#'   
#' }
