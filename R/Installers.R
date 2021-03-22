###dependencies
library(R.utils)
#' .OpenConfigFile
#' Internal. save and restore config info
#' @param config config list structure with software$softName...
#'
.OpenConfigFile <- function(config){
  config.file <- file.path( R.home(),"GeneFusionPkg.RData")
  if(file.exists(config.file)){
    config <- readRDS(config.file)
  }else{
    config <- NULL
  }
  if(missing(config)){
    ##abrir config file
   return(config)
  }else{##guardar config
    saveRDS(congif,config.file)
  }
  return(invisible(config))
}

InstallArribaR <- function(extDir){
  if(missing(extDir)){
    extDir <- "./Software"
  }else{
    extDir <- normalizePath(file.path(extDir,"Software"))
  }
  software <- list()
  software$mainPath <- extDir
  software <- .OpenConfigFile(software)
  ##Install STAR
  .InstallSTAR()
  .DownloadAssemblies()
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
  softare$star$alignmentPrefix <- "_STAR_Aligned.bam"
  software$star$aligmentPrefixSorted <- "_STAR_AlignedSortedByCoord.bam"
  .OpenConfigFile(config = software)
}

.DownloadAssemblies <- function(){
  cat("\nDownloading the STAR 2.7.6a assemblies")
  cat("\nDownloading GRCh37 Homo Sapiens...")

  software <- .OpenConfigFile()
  software$assembly <- file.path(software$mainPath,"Assemblies")

  # assembly.destination <- "/media/respaldo4t/Assemblies"
  # annotation.destination <- "/media/respaldo4t/Annotation"
  download.file(url = "ftp://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz",
                method = "wget",
                destfile = normalizePath(file.path(software$assembly,"GRCh37.fa.gz")))
  gunzip(filename = normalizePath(file.path(software$assembly,"GRCh37.fa.gz")),
        remove = T)

  argum <- paste0("sed -e ", shQuote("s/^MT\\t/chrM\\t/")," -e ",shQuote("s/^\\([1-9XY]\\|[12][0-9]\\)\\t/\\1\\t/"), " ",
                  normalizePath(file.path(software$assembly,"GRCh37.fa")))
  system(command = argum, ignore.stdout = T)

  if(file.exists(normalizePath(file.path(software$assembly,"GRCh37.fa")))){
    software$arriba$assemblyVersion <- "GRCh37"
    cat("\nAssembly GRCh37 installed")
  }else{
    stop("ERROR GRCh37")
  }
  software <- .OpenConfigFile(software)
  ###Annotation
  software$annotation <- file.path(software$mainPath,"Annotation")
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

  if(is.null(software$star$path)){
    stop("STAR not installed")
  }
  cat("\nBuilding STAR 2.7.6a - ARRIBA Genome Index...This may take some time")
  thr <- ifelse(parallel::detectCores() > 3, parallel::detectCores()-1, parallel::detectCores())
  cat(paste0("\nRunning STAR genomeGenerate with ",thr," CPU cores"))
  software$arriba$genomeDir <- file.path(sofware$assemblies,"STAR_GRCh37_GENCODE19_index")

  system2(command = software$star$command,
         args = c("--runMode genomeGenerate",
                  paste0("--genomeDir ",software$arriba$genomeDir),
                  paste0("--genomeFastaFiles ",sofware$assemblies,"/GRCh37.fa"),
                  paste0("--sjdbGTFfile ",file.path(software$annotation,"GENCODE19.gtf")),
                  paste0("--runThreadN ",thr),
                  "--sjdbOverhang 250") )

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
    system2(command = software$arriba$command)
    .OpenConfigFile(software)
  }else{
    stop("ERROR arriba failed")
  }

}

# DetectGeneFusions <- function(sbjFQ1){{
#   sbjFQ1 <- "/media/respaldo4t/RNAseq/39738/39738_1.fastq.gz"
#
# }
#   system2(command = "/media/respaldo4t/Softwares/Arriba/arriba")
# }





