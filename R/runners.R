##Test Arriba pipeline (STAR + arriba + feature + MIXTURE)
#'RunSTAR
#' This function will run STAR software
#' @param sbjFile string with the full path and name of the xxxx_1_fastq or gz sequence file. This
#' function only support paired data
#' @export
#' @author Elmer A. Fernández
#' @examples
#' \dontrun{
#' test.subject <- GetArribaRTest()
#' bam.subject <- RunSTAR(test.subject)
#' }
RunSTAR <- function(sbjFile){
  software <- .OpenConfigFile()
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA
  cat(paste0("\nCurrent STAR version : ", software$star$version, "\n"))
  ##--------
  file1 <- sbjFile
  file2 <- stringr::str_replace(file1,"_1.fastq","_2.fastq")
  cat(paste0("\nSubject files\n",file1,"\n",file2,"\n"))
  if(!all(file.exists(file1,file2))){
    stop("ERROR some files not found")
  }
  out.file <- stringr::str_remove(file1,".gz")
  out.file <- stringr::str_replace(out.file,"_1.fastq",software$star$alignmentPrefix)

  system2(command = software$star$command,
          args = c("--runThreadN 8",
                   paste0("--genomeDir " ,software$arriba$genomeDir),
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
  return(out.file)
}
#' GetArribaRTest
#' This function will retrieve the test files coming with the arriba software
#' @usage test.files <- GetArribaRTest()
#' @return an string with the full path of the firs read file of the test
GetArribaRTest <- function(){
  software <- .OpenConfigFile()
  files <- list.files(paste0(software$arriba$path,"/test"),full.names = T)
  id.fastq <- which(stringr::str_detect(files,"read1.fastq.gz"))
  id.fastq2 <- which(stringr::str_detect(files,"read2.fastq.gz"))
if(length(id.fastq)>0){
  file.rename(files[id.fastq],stringr::str_replace(files[id.fastq],"read1.fastq.gz","read_1.fastq.gz"))
  file.rename(files[id.fastq2],stringr::str_replace(files[id.fastq2],"read2.fastq.gz","read_2.fastq.gz"))
}else{
  if(file.exists(paste0(software$arriba$path,"/test/read_1.fastq.gz"))){
    return(paste0(software$arriba$path,"/test/read_1.fastq.gz"))
  }
  return(NULL)
}
  if(file.exists(paste0(software$arriba$path,"/test/read_1.fastq.gz"))){
    return(paste0(software$arriba$path,"/test/read_1.fastq.gz"))
  }
  return(NULL)
}

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
  software <- .OpenConfigFile()
  ##-------- ESTO SE DEBE REEMPLAZAR CON LA LINEA DE ARRIBA

  ##--------
  # call arriba
  if(!file.exists(sbjBamFile)){
    stop(paste0("\nERROR ", basename(sbjBamFile), "NOT FOUND"))
  }

  fusion.file <- stringr::str_replace(sbjBamFile,software$star$alignmentPrefix,"_Fusions.tsv")
  fusion.discarded.file <- stringr::str_replace(sbjBamFile,software$star$alignmentPrefix,"_Fusions_discarded.tsv")
  if(!all(file.exists(fusion.file,fusion.discarded.file))){
    cat(paste0("\nFILES NOT FOUND\n",fusion.file,"\n", fusion.discarded.file,"\nPls use runSTAR(subject)\n"))
  }
   system2(command = software$arriba$command,
           args = c(paste0("-x ",sbjBamFile),
                    paste0("-o ",fusion.file),
                    paste0("-O ",fusion.discarded.file),
                    paste0("-a ",software$assembly,"/GRCh37.fa"),
                    paste0("-g ", software$annotation,"/GENCODE19.gtf"),
                    paste0("-b ", normalizePath(file.path(software$arriba$database,"blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz"))),
                    paste0("-k ", normalizePath(file.path(software$arriba$database,"known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz"))),
                    paste0("-t ", normalizePath(file.path(software$arriba$database,"known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz"))),
                    paste0("-p ", normalizePath(file.path(software$arriba$database,"protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3")))
           ))
  ##write as an excel file
  fusionsTable <- .ReadArribaFustionTable(fusion.file)
  openxlsx::write.xlsx(fusionsTable, stringr::str_replace(fusion.file,".tsv",".xlsx"))
  #return(fusionsTable) ##for continuous processing into R
  # if(!stringr::str_detect(sbjBamFile,software$star$alignmentPrefixSorted)){
  #   Rsamtools::sortBam(file = sbjBamFile,
  #                      destination = stringr::str_replace(sbjBamFile, software$star$aligmentPrefix,software$star$alignmentPrefixSorted))
  #   test.out <- stringr::str_replace(sbjBamFile, software$star$aligmentPrefix,software$star$alignmentPrefixSorted)
  #   Rsamtools::indexBam(test.out)
  # }
  return(fusionsTable)
}

#' runSortIndexBam
#' This function will built the sorted bam file and its index.
#' @param sbjBamFile string full path of aligned bam file see \code{\link[runSTAR]{runSTAR}}
#' @param remove boolean (default  TRUE), if unsorted bam file should be removed or not.
#' @export
#' @seealso \code{\link{Rsamtools::sortBam}}
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
    destination.file <- stringr::remove(
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

