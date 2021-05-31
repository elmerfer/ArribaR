library(data.table)
dataDir <- path.expand("~/dataexpo")
dataFiles <- dir(dataDir, pattern = "csv$", full.names = TRUE)
file <- "/media/respaldo4t/RNAseq/39738/39738_Fusions_discarded.tsv"
gene <- "exon"
# gene2 <- "LTBP1"

.MakeCommand <- function(gene1, gene2, file){
  if(!missing(gene1) & missing(gene2)){
    command <- sprintf(
      "grep -E -w '^%s' %s",gene1,
      paste(file, collapse = " ")
    )  
    return(command)
  }
  if(all(c(!missing(gene1), !missing(gene2)))){
    command <- sprintf(
      "egrep -E -w '%s|%s' %s",gene1,gene2,
      paste(file, collapse = " ")
    )  
    return(command)
  }
  if(missing(gene1) & !missing(gene2)){
    command <- sprintf(
      "grep -E -w '%s' %s",gene2,
      paste(file, collapse = " ")
    )  
    return(command)
  }
    
}

command <- .MakeCommand("FOXR2", file = file)


dt <- data.table::fread(file, nrows=2)
cn <- colnames(dt)
gene <- "EWSR1"
command <- sprintf(
  "grep -P -w  '(?=^((?!RP11).)*$)%s' %s",gene,
  paste(file, collapse = " ")
)
command
dt <- data.table::fread(cmd = command)
colnames(dt) <- stringr::str_remove_all(cn, "#")
View(dt)


SAF <- rtracklayer::readGFF("/media/respaldo4t/Paquetes/Software/Annotation/GENCODE19.gtf")
SAF <- subset(SAF, type == "exon" & gene_type %in% c("protein_coding","pseudogene","processed_transcript","IG_C_gene","antisense",
                                                      "TR_V_gene","polymorphic_pseudogene","TR_C_gene"))
head(SAF)
