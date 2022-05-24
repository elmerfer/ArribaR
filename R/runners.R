##Test Arriba pipeline (STAR + arriba + feature + MIXTURE)

.SetThreads <- function(nThreads){
  if(missing(nThreads)){
    nThreads <- parallel::detectCores()-2
  }
  if(nThreads>parallel::detectCores()){
    nThreads <- parallel::detectCores()-2
  }
  if(nThreads < 0){
    nThreads <- parallel::detectCores()
  }
  return(nThreads)
}

#'RunSTAR
#' This function will run STAR software
#' @param sbjFile string with the full path and name of the xxxx_1_fastq or gz sequence file. This
#' function only support paired data
#' @param nThreads number of CPUs threads 
#' @export
#' @author Elmer A. Fernández
#' @examples
#' \dontrun{
#' test.subject <- GetArribaRTest()
#' bam.subject <- RunSTAR(test.subject)
#' }
RunSTAR <- function(sbjFile, nThreads, assemblyVersion = c("37","38")){
  software <- .OpenConfigFile()
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA
  cat(paste0("\nCurrent STAR version : ", software$star$version, "\n"))
  assemblyVersion <- paste0("GRCh",match.arg(as.character(assemblyVersion), choices = c("37","38")))
  cat(paste0("\nChosen assembly version GRCh : ", assemblyVersion, "\n"))
  ##--------
  file1 <- sbjFile
  file2 <- stringr::str_replace(file1,"_1.fastq","_2.fastq")
  cat(paste0("\nSubject files\n",file1,"\n",file2,"\n"))
  if(!all(file.exists(file1,file2))){
    stop("ERROR some files not found")
  }
  out.file <- stringr::str_remove(file1,".gz")
  out.file <- stringr::str_replace(out.file,"_1.fastq",software$star$alignmentPrefix)
  
  
  nThreads <- .SetThreads(nThreads)
  
   if(nThreads < 4){
     message(paste0("STAR running on ",nThreads, " may be not optimal"))
   }else{
     message(paste0("STAR running on ",nThreads))
   }
  genomeDir <- paste0(software$assembly,"/STAR_",assemblyVersion)
  genomeDir <- ifelse(assemblyVersion == "GRCh37",paste0(genomeDir,"_GENCODE19_index"),
                      paste0(genomeDir,"_GENCODE38_index"))
  
  system2(command = software$star$command,
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
  attr(out.file,"assemblyVersion") <- assemblyVersion
  return(out.file)
}
#' #' GetArribaRTest
#' #' This function will retrieve the test files coming with the arriba software
#' #' @usage test.files <- GetArribaRTest()
#' #' @return an string with the full path of the firs read file of the test
#' GetArribaRTest <- function(){
#'   software <- G.OpenConfigFile()
#'   files <- list.files(paste0(software$arriba$path,"/test"),full.names = T)
#'   id.fastq <- which(stringr::str_detect(files,"read1.fastq.gz"))
#'   id.fastq2 <- which(stringr::str_detect(files,"read2.fastq.gz"))
#' if(length(id.fastq)>0){
#'   file.rename(files[id.fastq],stringr::str_replace(files[id.fastq],"read1.fastq.gz","read_1.fastq.gz"))
#'   file.rename(files[id.fastq2],stringr::str_replace(files[id.fastq2],"read2.fastq.gz","read_2.fastq.gz"))
#' }else{
#'   if(file.exists(paste0(software$arriba$path,"/test/read_1.fastq.gz"))){
#'     return(paste0(software$arriba$path,"/test/read_1.fastq.gz"))
#'   }
#'   return(NULL)
#' }
#'   if(file.exists(paste0(software$arriba$path,"/test/read_1.fastq.gz"))){
#'     return(paste0(software$arriba$path,"/test/read_1.fastq.gz"))
#'   }
#'   return(NULL)
#' }

#' RunArriba
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
RunArriba <- function(sbjBamFile){
  assemblyVersion <- attr(sbjBamFile,"assemblyVersion")
  if(is.null(assemblyVersion)){
    stop("\nPlease process the subject by RunSTAR")
  }
  
  software <- .OpenConfigFile()
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA

  ##--------
  # call arriba
  if(!file.exists(sbjBamFile)){
    stop(paste0("\nERROR ", basename(sbjBamFile), "NOT FOUND"))
  }

  fusion.file <- stringr::str_replace(sbjBamFile,software$star$alignmentPrefix,"_Fusions.tsv")
  fusion.discarded.file <- stringr::str_replace(sbjBamFile,software$star$alignmentPrefix,"_Fusions_discarded.tsv")
  
  genomeversion <- ifelse(assemblyVersion == "GRCh37", "hg19_hs37d5_","hg38_")
  
   system2(command = software$arriba$command,
           args = c(paste0("-x ",sbjBamFile),
                    paste0("-o ",fusion.file),
                    paste0("-O ",fusion.discarded.file),
                    paste0("-a ",software$assembly,ifelse(assemblyVersion == "GRCh37","/GRCh37.fa","/GRCh38.fa")),
                    paste0("-g ", software$annotation,ifelse(assemblyVersion == "GRCh37","/GENCODE19.gtf","/GENCODE38.gtf")),
                    paste0("-b ", file.path(software$arriba$database,
                                            paste0("blacklist_",genomeversion,assemblyVersion,"_v2.1.0.tsv.gz"))),
                    paste0("-k ", file.path(software$arriba$database,
                                            paste0("known_fusions_",genomeversion,assemblyVersion,"_v2.1.0.tsv.gz"))),
                    paste0("-t ", file.path(software$arriba$database,
                                            paste0("known_fusions_",genomeversion,assemblyVersion,"_v2.1.0.tsv.gz"))),
                    paste0("-p ", file.path(software$arriba$database,
                                            paste0("protein_domains_",genomeversion,assemblyVersion,"_v2.1.0.gff3")))
           ))
   
   if(!all(file.exists(fusion.file,fusion.discarded.file))){
     cat(paste0("\nFUSION FILES NOT FOUND\n",fusion.file,"\n", fusion.discarded.file))
   }
  ##write as an excel file
  fusionsTable <- .ReadArribaFusionTable(fusion.file)
  version.info <- data.frame(  STAR = c(software$star$version), 
                             ARRIBA = c(software$arriba$version, assemblyVersion,packageVersion("ArribaR")),
                             SAMPLE = c(stringr::str_remove(basename(fusion.file),"_Fusions.tsv")))
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

#' runSortIndexBam
#' This function will built the sorted bam file and its index.
#' @param sbjBamFile string full path of aligned bam file see \code{\link[runSTAR]{runSTAR}}
#' @param remove boolean (default  TRUE), if unsorted bam file should be removed or not.
#' @export
#' @seealso \code{\link[Rsamtools]{sortBam}}
RunSortIndexBam <- function(sbjBamFile, remove = T){
  software <- .OpenConfigFile()
  if(!file.exists(sbjBamFile)){
    stop(paste0("\nERROR ",sbjBamFile," NOT FOUND"))
  }
  if(!stringr::str_detect(sbjBamFile,".bam")){
    stop(paste0("\nERROR this is not a bam file",sbjBamFile,"\n"))
  }
  if(stringr::str_detect(sbjBamFile, software$star$alignmentPrefixSorted)){
    stop(paste0("\nThe file seems to be already soprted, pls verify"))
  }
  if(!stringr::str_detect(sbjBamFile,software$star$alignmentPrefix)){
    warning("The file seems not being aligned trhough ArribaR::runSTAR() function")
    destination.file <- stringr::str_remove(
      paste0(sbjBamFile, software$star$alignmentPrefixSorted), ".bam")
  }else{
    destination.file <- stringr::str_remove(
      stringr::str_replace(sbjBamFile, software$star$alignmentPrefix, software$star$alignmentPrefixSorted), ".bam")
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

#' RunAnalysis
#' 
#' An all in one line to run STAR+Arriba+MIXTURE 
#' @param sbjFile (character) full path to the first fasta file from a subject (subject_1.fastq or subject_1.fastq.gz ), the scond read file should be subject_2.fastq
#' @param feature (character)  "exon" or "gene". it will be overwritten with "gene" if immuneTME == TRUE
#' @param immuneTME (logical (default TRUE) should the \code{\link[MIXTURE]{MIXTURE}} algorithm be run?
#' @export
#' @return a list with the following slots:
#' 
#' Fusions : a Gene Fusion data (also in an excel file see \code{\link{RunArriba}})
#' 
#' Counts : a list as returned by \code{\link{GetCounts}}
#' 
#' iTME : immune tumor microenvironment objets as returned by \code{\link{GetImmuneContent}}
#' @examples 
#' \dontrun{
#' results <- RunAnalysis(sbjFile = subject_1.fastq.gz)
#' ##if a do not want to run TME and get the exon counts
#' results <- RunAnalysis(sbjFile = subject_1.fastq.gz, feature = "exon" , immuneTME = FALSE)
#' }
RunAnalysis <- function(sbjFile , feature , immuneTME = TRUE){
  bam.subject <- RunSTAR(sbjFile)
  Fusions <- RunArriba(bam.subject)
  ft <- GetCounts(bam.subject, feature = ifelse(immuneTME,"gene",feature))
  if(immuneTME){
    ct <- GetImmuneContent(ftc = ft, filter = 0, plot = FALSE)
  }
  return(invisible(list(Fusions = Fusions, Counts = ft, iTME = ct)))
}
