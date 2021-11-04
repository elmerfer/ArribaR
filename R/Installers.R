###dependencies
library(R.utils)
#' .OpenConfigFile
#' Internal. save and restore config info
#' @param config config list structure with software$softName...
#'
.OpenConfigFile <- function(config){
  config.file <- file.path( .libPaths()[1],"GeneFusionPkg.RData")
  if(file.exists(config.file)){
    if(missing(config)){##si lo llama vacio, devuleve lo que hay en el file
      config <- readRDS(config.file)
      return(invisible(config))
    }else{
      saveRDS(config,config.file)
      return(invisible(config))
    }
  }else{
    return(list())
  }
  
}
#' InstallArribaR
#' This function will install the arri software and all its dependencies like the roght STAR version
#' It currently install version 2.1.0 of arriba and 2.7.6a of STAR, the GRCh37 genome assembly and GENCODE19.gtf annotation
#' future installer will be more felxible
#' @param extDir character string of the instalation path. It should be write enable and hold enouph space
#' if missing, it will create on the current directory the folder "Software" and everything will be installed there
#' if given, it will create the folder givenpath/Software
InstallArribaR <- function(extDir){
  if(missing(extDir)){
    extDir <- "./Software"
  }else{
    extDir <- normalizePath(file.path(extDir,"Software"))
  }
  software <- .OpenConfigFile()
  software$mainPath <- extDir
  software <- .OpenConfigFile(software)
  ##Install STAR
  .InstallSTAR()
  .DownloadAssemblies()
  .DownloadAnnotation()
  .BuildSTARindex()
  .InstallArribaGeneFusion()
}

.InstallSTAR <- function(){
  software <- .OpenConfigFile()
  if(is.null(software)){
    stop("ERROR")
  }
  cat("\nDownloading the STAR aligner software version 2.7.6a compatible with Arriba version 2.0.1")
  tmp.destfile <- tempfile()
  download.file(url = "https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz",
                method = "wget",
                destfile = tmp.destfile)
  files <- untar(tarfile = tmp.destfile,list = T )
  star.f <- files[which(stringr::str_detect(files, "Linux"))]
  star.f[stringr::str_detect(star.f, "static")]
  str.comp <- which(unlist(stringr::str_split(star.f[stringr::str_detect(star.f, "static")][2],"/") )=="STAR")-1

  software$star$path <- file.path(software$mainPath,"STAR")
  software$star$command <- normalizePath(file.path(software$star$path,"STAR"))
  software$star$version <- "2.7.6a"

  untar(tarfile = tmp.destfile,
        files = star.f[stringr::str_detect(star.f, "static")][-1],
        exdir = software$star$path,
        extras = paste0("--strip-components ",str.comp))
  file.remove(tmp.destfile)

  if(file.exists(software$star$command)){
    cat("\nSTAR 2.7.6a installed")
  }
  system2(command = software$star$command)
  software$star$alignmentPrefix <- "_STAR_Aligned.bam"
  software$star$alignmentPrefixSorted <- "_STAR_AlignedSortedByCoordinates.bam"
  .OpenConfigFile(config = software)
}

.DownloadAssemblies <- function(){
  cat("\nDownloading the STAR 2.7.6a assemblies")
  cat("\nDownloading GRCh37 Homo Sapiens...")

  software <- .OpenConfigFile()
  software$assembly <- file.path(software$mainPath,"Assemblies")
  dir.create(software$assembly)

  # assembly.destination <- "/media/respaldo4t/Assemblies"
  # annotation.destination <- "/media/respaldo4t/Annotation"
  download.file(url = "ftp://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz",
                method = "wget",
                destfile = normalizePath(file.path(software$assembly,"GRCh37.fa.gz")))
  R.utils::gunzip(filename = normalizePath(file.path(software$assembly,"GRCh37.fa.gz")),
        remove = T)

  # argum <- paste0("sed -e ", shQuote("s/^MT\\t/chrM\\t/")," -e ",shQuote("s/^\\([1-9XY]\\|[12][0-9]\\)\\t/\\1\\t/"), " ",
  #                 normalizePath(file.path(software$assembly,"GRCh37.fa")))
  # system(command = argum, ignore.stdout = T)

  if(file.exists(normalizePath(file.path(software$assembly,"GRCh37.fa")))){
    software$arriba$assemblyVersion <- "GRCh37"
    cat("\nAssembly GRCh37 installed")
  }else{
    stop("ERROR GRCh37")
  }
  .OpenConfigFile(software)

}


.DownloadAnnotation <-function(){
  software <- .OpenConfigFile()
  ###Annotation
  software$annotation <- file.path(software$mainPath,"Annotation")
  if(!dir.exists(software$annotation)){
    cat(paste0("\nCreating ",software$annotation))
    dir.create(software$annotation)
    if(dir.exists(software$annotation)){
      cat("---> OK\n")
    }
  }
  cat("\nDownloading GRCh37 Annotation Homo Sapiens...")
  arg2 <- paste0("wget --load-cookies /tmp/cookies.txt \"https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1qa9YrM-hOwzKgS8tZcU3XlTBulHnMORO' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p'",")","&id=1qa9YrM-hOwzKgS8tZcU3XlTBulHnMORO\"")
  arg2 <- paste(arg2," -O ",normalizePath(file.path(software$annotation,"GENCODE19.gtf.zip"))," && rm -rf /tmp/cookies.txt")
  system(command = arg2)
  if(file.exists(normalizePath(file.path(software$annotation,"GENCODE19.gtf.zip")))){
    cat("\nUncompressing Annotation file")
    zip::unzip(zipfile = normalizePath(file.path(software$annotation,"GENCODE19.gtf.zip")),exdir = software$annotation )
  }else{
    stop("Annotation file download failed")
  }
  if(file.exists(normalizePath(file.path(software$annotation,"GENCODE19.gtf")))){
    cat("\nAnnotation file saved as GENECODE19.gtf")
    file.remove(normalizePath(file.path(software$annotation,"GENCODE19.gtf.zip")))
    software <- .OpenConfigFile(software)
  }else{
    stop("ERROR annotation")
  }
}

.BuildSTARindex <- function(){
  software <- .OpenConfigFile()
  if(is.null(software$star$path)){
    stop("STAR not installed")
  }
  cat("\nBuilding STAR 2.7.6a - ARRIBA Genome Index...This may take some time")
  thr <- ifelse(parallel::detectCores() > 3, parallel::detectCores()-1, parallel::detectCores())
  cat(paste0("\nRunning STAR genomeGenerate with ",thr," CPU cores"))
  software$arriba$genomeDir <- file.path(software$assembly,"STAR_GRCh37_GENCODE19_index")
  
  system2(command = software$star$command,
          args = c("--runMode genomeGenerate",
                   paste0("--genomeDir ",software$arriba$genomeDir),
                   paste0("--genomeFastaFiles ",software$assembly,"/GRCh37.fa"),
                   paste0("--sjdbGTFfile ",file.path(software$annotation,"GENCODE19.gtf")),
                   paste0("--runThreadN ",thr),
                   "--sjdbOverhang 250") )
  .OpenConfigFile(software)
}

.InstallArribaGeneFusion <- function(){

  software <- .OpenConfigFile()
  if(is.null(software)){
    stop("Error ARRIBA")
  }

  software$arriba$path <- file.path(software$mainPath,"Arriba")
  software$arriba$version <- "2.1.0"
  tmp.destfile <- tempfile()
    # "/media/respaldo4t/Softwares/arriba_v2.1.0.tar.gz"
  download.file(url = "https://github.com/suhrig/arriba/releases/download/v2.1.0/arriba_v2.1.0.tar.gz",
                method = "wget",
                destfile = tmp.destfile)

  arriba.files <- untar(tmp.destfile, list =T)
  test.files <- arriba.files[stringr::str_detect(arriba.files,"test")]

  untar(tmp.destfile, files = test.files,
        exdir = software$arriba$path,
        extras = paste0("--strip-components ",1))

  database.files <- arriba.files[stringr::str_detect(arriba.files,"database")]

  software$arriba$database <- file.path(software$arriba$path, "database")

  untar(tmp.destfile, files = database.files,
        exdir = software$arriba$path,
        extras = paste0("--strip-components ",1))
  if(all(file.exists(paste0(software$arriba$database,basename(database.files[-1]))))){
    cat("Arriba database created")
  }

  arriba.soft <- arriba.files[stringr::str_detect(arriba.files,"/arriba")]
  arriba.soft <- arriba.soft[!stringr::str_detect(arriba.soft,"source")]
  untar(tmp.destfile, files = arriba.soft,
        exdir = software$arriba$path,
        extras = paste0("--strip-components ",1))
  software$arriba$command <- file.path(software$arriba$path,"arriba")
  if(file.exists(software$arriba$command)){
    system2(command = software$arriba$command, args = "-h")
    .OpenConfigFile(software)
  }else{
    stop("ERROR arriba failed")
  }

}

##########3
.InstallArribaOnGenomeDB <- function(){
  
  software <- GenomeDB:::.OpenConfigFile()
  if(is.null(software)){
    stop("Error ARRIBA")
  }
  if(!any(stringr::str_detect(names(software$Software), "STAR"))){
    stop("STAR should be previously installed, use see Aligners library at github/elmerfer/Aligners")
  }
  if(any(stringr::str_detect(names(software$Software), "ARRIBA"))){
    stop("ARRIBA already installed")
  }
  
  # software$arriba$path <- file.path(software$mainPath,"Arriba")
  software$Software$ARRIBA$main <-  file.path(software$Software$main,"ARRIBA")#to install softwar versions and index files
  stopifnot(dir.create(software$Software$ARRIBA$main))
  
  
  software$Sotware$ARRIBA$version <- "2.1.0"
  tmp.destfile <- tempfile()
  # "/media/respaldo4t/Softwares/arriba_v2.1.0.tar.gz"
  download.file(url = "https://github.com/suhrig/arriba/releases/download/v2.1.0/arriba_v2.1.0.tar.gz",
                method = "wget",
                destfile = tmp.destfile, quite=T)
  
  arriba.files <- untar(tmp.destfile, list =T)
  test.files <- arriba.files[stringr::str_detect(arriba.files,"test")]
  
  untar(tmp.destfile, files = test.files,
        exdir = software$Software$ARRIBA$main,
        extras = paste0("--strip-components ",1))
  
  database.files <- arriba.files[stringr::str_detect(arriba.files,"database")]
  
  software$Software$ARRIBA$database <- file.path(software$Software$ARRIBA$main, "database")
  
  untar(tmp.destfile, files = database.files,
        exdir = software$Software$ARRIBA$main,
        extras = paste0("--strip-components ",1))
  if(all(file.exists(paste0(software$Software$ARRIBA$database,basename(database.files[-1]))))){
    cat("Arriba database created")
  }
  
  arriba.soft <- arriba.files[stringr::str_detect(arriba.files,"/arriba")]
  arriba.soft <- arriba.soft[!stringr::str_detect(arriba.soft,"source")]
  untar(tmp.destfile, files = arriba.soft,
        exdir = software$Software$ARRIBA$main,
        extras = paste0("--strip-components ",1))
  software$Software$ARRIBA$command <- file.path(software$Software$ARRIBA$main,"arriba")
  if(file.exists(software$Software$ARRIBA$command)){
    system2(command = software$Software$ARRIBA$command, args = "-h")
    GenomeDB:::.OpenConfigFile(software)
  }else{
    stop("ERROR arriba failed")
  }
  
}




