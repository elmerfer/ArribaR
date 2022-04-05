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
InstallArribaR <- function(){
  if(require(GenomeDB)==FALSE){
    stop("Pls you should install GenomeDB package See https://github.com/elmerfer/GenomeDB")
  }
  if(require(Aligners)==FALSE){
    stop("Pls you should install Aligners package See https://github.com/elmerfer/Aligners")
  }
  
  software <- GenomeDB:::.OpenConfigFile()
  
  ##Install STAR
  if(c("STAR" %in% names(software$Software)) == FALSE ){
    Aligners::InstallAligners()  
  }
  
  # 
  # .DownloadAssemblies()
  # .DownloadAnnotation()
   Aligners::UpdateSTARindex()
   .InstallArribaOnGenomeDB()
}


.InstallArribaGeneFusion <- function(){

  software <- GenomeDB:::.OpenConfigFile()
  if(is.null(software)){
    stop("Error ARRIBA")
  }
  
  software$Software$ARRIBA$main <- file.path(software$Software$main,"ARRIBA")
  
  software$Software$ARRIBA$version <- "2.1.0"
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

  software$Software$ARRIBA$database <- file.path(software$Software$ARRIBA$main, "database")

  untar(tmp.destfile, files = database.files,
        exdir = software$arriba$path,
        extras = paste0("--strip-components ",1))
  if(all(file.exists(paste0(software$Software$ARRIBA$database,basename(database.files[-1]))))){
    cat("Arriba database created")
  }

  arriba.soft <- arriba.files[stringr::str_detect(arriba.files,"/arriba")]
  arriba.soft <- arriba.soft[!stringr::str_detect(arriba.soft,"source")]
  untar(tmp.destfile, files = arriba.soft,
        exdir = software$arriba$path,
        extras = paste0("--strip-components ",1))

  software$Software$ARRIBA$command <- file.path(software$Software$ARRIBA$main,"Arriba")
  if(file.exists(software$Software$ARRIBA$command)){
    system2(command = software$Software$ARRIBA$command, args = "-h")
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
  
  
  software$Software$ARRIBA$version <- "2.1.0"
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




