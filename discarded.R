library(data.table)
dataDir <- path.expand("~/dataexpo")
dataFiles <- dir(dataDir, pattern = "csv$", full.names = TRUE)

GS <- openxlsx::read.xlsx("/media/respaldo4t/FLENI/Fusions/FusionPanel.xlsx")

file <- "/media/respaldo4t/RNAseq/44195/44195__Fusions_discarded.tsv"

cn <- colnames(fread(file, nrows=1,header = T))
idtags <- which(cn =="tags")
ss <- system2( command="awk", args= sprintf(' \'$17=="Mitelman" \' %s',fileFusions),
               stdout = TRUE)
ret <- plyr::ldply(ss, function(x){
  unlist(stringr::str_split(x,"\t"))
})

panel.par <- function(gene1,gene2, fileFusions) {
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
  
panel.gene1 <- function(gene1, fileFusions){
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
  
  cat(paste0(ifelse(is.null(nrow(ret)),0,nrow(ret))," 5' genes found"))
  ret
}

panel.gene2 <- function(gene2, fileFusions){
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
  
  cat(paste0(ifelse(is.null(nrow(ret)),0,nrow(ret))," 5' genes found"))
  ret
}


SearchPanel <- function(genePanel, fileFusions, save=T){
  cnames <- stringr::str_remove_all(colnames(data.table::fread(file = fileFusions,nrows=1)),"#")
  pares <- na.omit(genePanel[,c("gene1","gene2"),drop=F])
  gene1 <- genePanel[is.na(genePanel$gene2),c("gene1","gene2"),drop=F]
  gene2 <- genePanel[is.na(genePanel$gene1),c("gene1","gene2"),drop=F]
  if(nrow(pares)>0){
    ret.p<- plyr::ldply(1:nrow(pares), function(px){
      aux <- panel.par(pares$gene1[px],pares$gene2[px], fileFusions = fileFusions)
      if(all(is.na(aux[[1]]))) return(NULL)
      aux
    })
    colnames(ret.p) <- cnames
  }else{
    ret.p<-NA
  }
   if(nrow(gene1)>0){
     ret.g1 <- plyr::ldply(1:nrow(gene1), function(px){
       aux <- panel.gene1(pares$gene1[px], fileFusions = fileFusions)
       if(all(is.na(aux[[1]]))) return(NULL)
       aux
     })
     colnames(ret.g1) <- cnames
   }else{
     ret.g1<-NA
   } 
  if(nrow(gene2)>0){
     ret.g2 <- plyr::ldply(1:nrow(gene2), function(px){
       aux<-panel.gene2(pares$gene2[px], fileFusions = fileFusions)
       if(all(is.na(aux[[1]]))) return(NULL)
       aux
     })
     colnames(ret.g2) <- cnames
   }else{
     ret.g2<-NA
   }
  pan <- list(Pairs=ret.p,"gene1-any"=ret.g1,"any-gene2 '"=ret.g2)
  names(pan) <- c("Pairs","gene1-any","any-gene2")
  file.name <- stringr::str_replace(fusionsFile,"_Fusions_discarded.tsv","_Recovered.xlsx")
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

pan <- SearchPanel(GS,"/media/respaldo4t/RNAseq/44195/44195__Fusions_discarded.tsv")
names(pan) <- c("Pairs","gene1-any","any-gene2")
openxlsx::write.xlsx(pan, file = "Test.xlsx", overwrite = T)

names(pan)

panel.par(pares$gene1[12],pares$gene2[12], "/media/respaldo4t/RNAseq/44195/44195__Fusions_discarded.tsv")
head(pares)

p1 <- panel.gene1(gen1,  "/media/respaldo4t/RNAseq/44195/44195__Fusions_discarded.tsv")

View(aa)
lapply(ugenes, function(x){
  subset(GS, gene1 == x | gene2==x)[,c("gene1","gene2")]
})

View(GS)
lg <- lapply(1:nrow(GS), function(x)
       c(GS$gene1[x], GS$gene2[x])
       )


library(plyr)
sbjF <- "/media/respaldo4t/RNAseq/44129/44129__Fusions_discarded.tsv"
dt <- data.table::fread(file = sbjF, nrows=2)
cn <- colnames(dt)
ldf <- ldply(lg, function(x){
    dt <- data.table::fread(.MakeCommand(x[1], file = sbjF))
    if(nrow(dt)==0){
      return(dt)
    }else{
      if(!is.na(x[2])){
        dt<- subset(dt,V2==x[2])
      }
    }
    
    dt
})
id <- apply(ldf[,c("split_reads1", "split_reads2", "discordant_mates")], MARGIN = 1, FUN = function(x) any(x>0))
colnames(ldf) <- stringr::str_remove_all(cn,"#")
# head(ldf)
openxlsx::write.xlsx(ldf[id,],stringr::str_replace_all(sbjF,".tsv",".xlsx"),overwrite = T)
# openxlsx::write.xlsx(ldf[!id,],"44129_DiscardedNosplit.xlsx",overwrite = T)


lapply(ldf, dim)

command <- .MakeCommand("ATG7", file = "/media/respaldo4t/RNAseq/44129/44129__Fusions_discarded.tsv")

data.table::fread(.MakeCommand("ATG7", file = "/media/respaldo4t/RNAseq/44129/44129__Fusions_discarded.tsv"))
data.table::fread("grep -E -w ^ATG7 /media/respaldo4t/RNAseq/44129/44129__Fusions_discarded.tsv")[,1:4]





dt <- data.table::fread(file, nrows=2)
cn <- colnames(dt)
gene <- "CUX1"
command <- sprintf(
  "grep -P -w  '(?=^((?!RP11).)*$)%s' %s",gene,
  paste(file, collapse = " ")
)
command
dt <- data.table::fread(command)
colnames(dt) <- stringr::str_remove_all(cn, "#")
View(dt)


SAF <- rtracklayer::readGFF("/media/respaldo4t/Paquetes/Software/Annotation/GENCODE19.gtf")
SAF <- subset(SAF, type == "exon" & gene_type %in% c("protein_coding","pseudogene","processed_transcript","IG_C_gene","antisense",
                                                      "TR_V_gene","polymorphic_pseudogene","TR_C_gene"))
head(SAF)
