#' @importFrom Rbowtie2 bowtie2_build
#' @importFrom Biostrings masks DNAStringSet injectHardMask
#' @importFrom BSgenome available.genomes
#' @importFrom BSgenome installed.genomes
#' @importFrom BiocManager install
setClass(Class = ".ConfigClass",
         slots = list(
             configList = "list",
             validAttr = "list",
             validWriteAttr = "list",
             curOrgDb = "character",
             curTxDb = "character"
             )
         )


setMethod(f = "initialize",
          signature = ".ConfigClass",
          definition = function(.Object,...){
              .Object@configList<-list(threads=1,tmpdir=".",refdir=NULL,genome=NULL,knownGene=NULL,bsgenome=NULL,annoDb=NULL,bt2Idx=NULL,DHS=NULL,blacklist=NULL,SNP=NULL)
              .Object@validAttr<-list(threads="numeric",tmpdir="character",refdir="character",genome="character",knownGene="TxDb",bsgenome="BSgenome",annoDb="OrgDb",bt2Idx="character",DHS="character",blacklist="character",SNP="character")
              .Object@validWriteAttr<-list(threads="numeric",tmpdir="character",refdir="character",genome="character")
              .Object
          })


setGeneric(name = "isValidAttr",
           def = function(.Object,item,...){
               standardGeneric("isValidAttr")
           })
setMethod(f = "isValidAttr",
          signature = ".ConfigClass",
          definition = function(.Object,item,...){
              if(is.null(.Object@validAttr[[item]])){
                  stop(paste(item,"is not an attribute"))
              }
          })

setGeneric(name = "isWriteValidAttr",
           def = function(.Object,item,...){
               standardGeneric("isWriteValidAttr")
           })
setMethod(f = "isWriteValidAttr",
          signature = ".ConfigClass",
          definition = function(.Object,item,...){
              if(is.null(.Object@validWriteAttr[[item]])){
                  stop(paste(item,"is not a writeable attribute"))
              }
          })

setGeneric(name = "isValidVal",
           def = function(.Object,item,val,...){
               standardGeneric("isValidVal")
           })
setMethod(f = "isValidVal",
          signature = ".ConfigClass",
          definition = function(.Object,item,val,...){
              isWriteValidAttr(.Object,item)
              if(is.null(val)){
                  stop("val can not be NULL")
              }
              if(.Object@validAttr[[item]]!=class(val)){
                  stop(paste(item,"is requied to be",.Object@validAttr[[item]],",\"",val,"\" is ",class(val)))
              }
              if(item=="refdir"||item=="tmpdir"){
                  checkPathExist0(.Object,val)
                  val<-normalizePath(val)
              }
              if(item=="genome"){
                  msgBoxBegin()
                  message("Reference data configuraion start ...")
                  message("Configure bsgenome ...")
                  bsgenomename<- BSgenome::available.genomes()[grepl(paste0(val,"$"),available.genomes())]
                  if(length(bsgenomename)==0){
                      stop("bsgenome not")
                  }
                  
                  bsgenomeinstall <- BSgenome::installed.genomes()[grepl(paste0(val,"$"),installed.genomes())]
                  
                  if(length(bsgenomeinstall)==0){
                      BiocManager::install(bsgenomename)
                  }
                  .Object@configList[["bsgenome"]]<-getBSgenome(val)
                  if(is.null(.Object@configList[["refdir"]])){
                      stop("'refdir' should be configured before 'genome'")
                  }
                  fileprefix<-file.path(.Object@configList[["refdir"]],val)
                  ##annoDb:orgdb
                  message("Configure annoDb ...")
                  tryCatch({
                      .Object@configList[["annoDb"]]<-GetOrgDb(.Object,val)
                  },
                  error = function(e){
                      message(as.character(e))
                      stop(paste("To install the package, type:",
                                "library(\"BiocManager\")",
                                sprintf("BiocManager::install(\"org.%s.eg.db\")",unlist(strsplit(as.character(e),"org\\.|\\.eg\\.db"))[2]),
                                sep="\n"))
                  })
                  message("Configure knownGene ...")
                  tryCatch({
                      .Object@configList[["knownGene"]]<-GetTxDb(.Object,val)
                  },
                  error = function(e){
                      message(as.character(e))
                      stop(paste("To install the package, type:",
                                 "library(\"BiocManager\")",
                                 sprintf("BiocManager::install(\"TxDb.%sGene\")",unlist(strsplit(as.character(e),"TxDb\\.|Gene"))[2]),
                                 sep="\n"))
                  })
                  
                  #genome fasta
                  message("Generate genome fasta file ...")
                  fastaFilePath<-paste0(fileprefix,".fa")
                  fastaFilePathlock<-paste0(fileprefix,".fa.lock")
                  if(file.exists(fastaFilePathlock)){
                      unlink(fastaFilePath)
                      unlink(fastaFilePathlock)
                  }
                  if(!file.exists(fastaFilePath)){
                      file.create(fastaFilePathlock)
                      BSgenomeSeqToFasta(.Object,.Object@configList[["bsgenome"]],fastaFilePath)
                      unlink(fastaFilePathlock)
                  }

                  #bowtie2 index
                  message("Configure bowtie2 index ...")
                  fileprefixlock<-paste0(fileprefix,".fa.bt2.lock")
                  if(file.exists(fileprefixlock)){
                      unlink(paste0(fileprefix,".1.bt2"))
                      unlink(paste0(fileprefix,".2.bt2"))
                      unlink(paste0(fileprefix,".3.bt2"))
                      unlink(paste0(fileprefix,".4.bt2"))
                      unlink(paste0(fileprefix,".rev.1.bt2"))
                      unlink(paste0(fileprefix,".rev.2.bt2"))
                      unlink(fileprefixlock)
                  }
                  if(!(file.exists(paste0(fileprefix,".1.bt2"))&&file.exists(paste0(fileprefix,".2.bt2"))&&
                       file.exists(paste0(fileprefix,".3.bt2"))&&file.exists(paste0(fileprefix,".4.bt2"))&&
                       file.exists(paste0(fileprefix,".rev.1.bt2"))&&file.exists(paste0(fileprefix,".rev.1.bt2")))){
                      message("Build bowtie2 index ...")
                      message("It may take more than one hour with single thread. Please wait.")
                      message("Using multi-threads will be faster.")
                      file.create(fileprefixlock)
                      message(paste("--threads",as.character(.obtainConfigure("threads"))))
                      bowtie2_build(fastaFilePath,fileprefix,"-q","--threads",as.character(.obtainConfigure("threads")),overwrite=TRUE)
                      unlink(fileprefixlock)
                  }
                  .Object@configList[["bt2Idx"]]<-fileprefix
                  #kownGene
                  # message("Configure knownGene ...")
                  # knownGeneFilePath<-paste0(fileprefix,".knownGene.sqlite")
                  # knownGeneFilePathlock<-paste0(fileprefix,".knownGene.sqlite.lock")
                  # if(file.exists(knownGeneFilePathlock)){
                  #     unlink(knownGeneFilePath)
                  #     unlink(knownGeneFilePathlock)
                  # }
                  # if(!file.exists(knownGeneFilePath)){
                  #     file.create(knownGeneFilePathlock)
                  #     if(sum(val==c("hg19","hg38","mm10","mm9"))){
                  #           saveDb(makeTxDbFromUCSC(genome=val, tablename="knownGene"), file=knownGeneFilePath)
                  #     }else{
                  #         saveDb(makeTxDbFromUCSC(genome=val, tablename="refGene"), file=knownGeneFilePath) 
                  #     }
                  #     unlink(knownGeneFilePathlock)
                  # }
                  # .Object@configList[["knownGene"]]<-loadDb(knownGeneFilePath)
                  
                  blacklistFilePath<-paste0(fileprefix,".blacklist.bed")
                  DHSFilePath<-paste0(fileprefix,".DHS.bed")
                  if(sum(val==c("hg19","hg38","mm10","mm9"))){
                      message("Configure blacklist and DHS ...")
                      downloadFilePathlock<-paste0(fileprefix,".download.lock")
                      if(file.exists(downloadFilePathlock)){
                          unlink(blacklistFilePath)
                          unlink(DHSFilePath)
                          unlink(downloadFilePathlock)
                      }
                      if(!file.exists(blacklistFilePath)||!file.exists(blacklistFilePath)){
                          file.create(downloadFilePathlock)
                          #DHS
                          if(!file.exists(blacklistFilePath)){
                              download.file(url = sprintf("https://wzthu.github.io/esATAC/refdata/%s.DHS.bed.gz",val),
                                            destfile = paste0(DHSFilePath,".gz"),method = getOption("download.file.method"))
                              gunzip(paste0(DHSFilePath,".gz"),destname=DHSFilePath,overwrite=TRUE)
                          }
                          #blacklist
                          if(!file.exists(blacklistFilePath)){
                              download.file(url = sprintf("https://wzthu.github.io/esATAC/refdata/%s.blacklist.bed.gz",val),
                                            destfile = paste0(blacklistFilePath,".gz"),method = getOption("download.file.method"))
                              gunzip(paste0(blacklistFilePath,".gz"),destname=blacklistFilePath,overwrite=TRUE)
                          }
                          unlink(downloadFilePathlock)
                      }
                      #.Object@configList[["DHS"]]<-DHSFilePath
                      #.Object@configList[["blacklist"]]<-blacklistFilePath
                  }
                  if(file.exists(blacklistFilePath)){
                      .Object@configList[["blacklist"]]<-blacklistFilePath
                  }
                  if(file.exists(DHSFilePath)){
                      .Object@configList[["DHS"]]<-DHSFilePath
                  }
                  snpFilePath<-paste0(fileprefix,".snp.txt")
                  if(val=="hg19"){
                      message("Configure SNP ...")
                      downloadFilePathlock<-paste0(fileprefix,".download1.lock")
                      if(file.exists(downloadFilePathlock)){
                          unlink(downloadFilePathlock)
                      }
                      if(!file.exists(snpFilePath)){
                          file.create(downloadFilePathlock)
                          download.file(url = sprintf("https://wzthu.github.io/esATAC/refdata/%s.snp.txt.gz",val),
                                        destfile = paste0(snpFilePath,".gz"),method = getOption("download.file.method"))
                          gunzip(paste0(snpFilePath,".gz"),destname=snpFilePath,overwrite=TRUE)
                          unlink(downloadFilePathlock)
                      }
                      #.Object@configList[["SNP"]]<-snpFilePath
                  }
                  if(file.exists(snpFilePath)){
                      .Object@configList[["SNP"]]<-snpFilePath
                  }
                  

                  message("Reference data configuraion done")
                  msgBoxDone()
                  
              }
              
              .Object
          })


setGeneric(name = "checkPathExist0",
           def = function(.Object,filePath,...){
               standardGeneric("checkPathExist0")
           })
setMethod(f = "checkPathExist0",
          signature = ".ConfigClass",
          definition = function(.Object,filePath,...){
              if(!is.null(filePath)){
                  if(!dir.exists(dirname(filePath))){
                      stop(paste("error, path does not exist:",filePath))
                  }
              }
          })

setGeneric(name = "BSgenomeSeqToFasta",
           def = function(.Object,bsgenome, outFile,...){
               standardGeneric("BSgenomeSeqToFasta")
           })
setMethod(f = "BSgenomeSeqToFasta",
          signature = ".ConfigClass",
          definition = function(.Object,bsgenome, outFile,...){
              if(!is(bsgenome, "BSgenome")){stop("The variable 'bsgenome' is not a BSgenome")}
              append <- FALSE
              for(chrT in seqnames(bsgenome)){
                  if(is.null(masks(bsgenome[[chrT]])))
                      chrSeq <- DNAStringSet(bsgenome[[chrT]])
                  else
                      chrSeq <- DNAStringSet(injectHardMask(bsgenome[[chrT]], letter="N"))
                  names(chrSeq) <- chrT
                  writeXStringSet(chrSeq, filepath=outFile, format="fasta", append=append)
                  append <- TRUE
              }
              return(outFile)
          })

setGeneric(name = "GetTxDb",
           def = function(.Object,genome,...){
               standardGeneric("GetTxDb")
           })
setMethod(f = "GetTxDb",
          signature = ".ConfigClass",
          definition = function(.Object,genome,...){
              if(genome == "hg19"){
                  .Object@curTxDb <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
              }else if(genome == "hg38"){
                  .Object@curTxDb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
              }else if(genome == "mm9"){
                  .Object@curTxDb <- "TxDb.Mmusculus.UCSC.mm9.knownGene"
              }else if(genome == "mm10"){
                  .Object@curTxDb <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
              }else if(genome == "danRer10"){
                  .Object@curTxDb <- "TxDb.Drerio.UCSC.danRer10.refGene"
              }else if(genome == "galGal5"){
                  .Object@curTxDb <- "TxDb.Ggallus.UCSC.galGal5.refGene"
              }else if(genome == "galGal4"){
                  .Object@curTxDb <- "TxDb.Ggallus.UCSC.galGal4.refGene"
              }else if(genome == "rheMac3"){
                  .Object@curTxDb <- "TxDb.Mmulatta.UCSC.rheMac3.refGene"
              }else if(genome == "rheMac8"){
                  .Object@curTxDb <- "TxDb.Mmulatta.UCSC.rheMac8.refGene"
              }else if(genome == "rn5"){
                  .Object@curTxDb <- "TxDb.Rnorvegicus.UCSC.rn5.refGene"
              }else if(genome == "rn6"){
                  .Object@curTxDb <- "TxDb.Rnorvegicus.UCSC.rn6.refGene"
              }else if(genome == "sacCer3"){
                  .Object@curTxDb <- "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"
              }else if(genome == "sacCer2"){
                  .Object@curTxDb <- "TxDb.Scerevisiae.UCSC.sacCer2.sgdGene"
              }else if(genome == "susScr3"){
                  .Object@curTxDb <- "TxDb.Sscrofa.UCSC.susScr3.refGene"
              }else {
                  warning(paste0("TxDb Annotation package does not support for ",genome))
                  return(NULL)
              }
              tryCatch({
                  library(.Object@curTxDb,character.only = TRUE)
              },
              error = function(e){
                  BiocManager::install(.Object@curTxDb)
                  library(.Object@curTxDb,character.only = TRUE)
              })
              return(get0(.Object@curTxDb))
          })

setGeneric(name = "GetOrgDb",
           def = function(.Object,genome,...){
               standardGeneric("GetOrgDb")
           })
setMethod(f = "GetOrgDb",
          signature = ".ConfigClass",
          definition = function(.Object,genome,...){
              if(genome == "hg19"||genome == "hg38"){
                  .Object@curOrgDb <- "org.Hs.eg.db"
              }else if(genome == "mm10" || genome == "mm9"){
                  .Object@curOrgDb <- "org.Mm.eg.db"
              }else if(genome == "danRer10"){
                  .Object@curOrgDb <- "org.Dr.eg.db"
              }else if(genome == "galGal5" || genome == "galGal4"){
                  .Object@curOrgDb <- "org.Gg.eg.db"
              }else if(genome == "rheMac3" || genome == "rheMac8"){
                  .Object@curOrgDb <- "org.Mmu.eg.db"
              }else if(genome == "panTro4" ){
                  .Object@curOrgDb <- "org.Pt.eg.db"
              }else if(genome == "rn6" || genome == "rn5"){
                  .Object@curOrgDb <- "org.Rn.eg.db"
              }else if(genome == "sacCer3" || genome == "sacCer2"){
                  .Object@curOrgDb <- "org.Sc.sgd.db"
              }else if(genome == "susScr3"){
                  .Object@curOrgDb <- "org.Ss.eg.db"
              }else {
                  warning(paste0("OrgDb Annotation package does not support for ",genome))
                  return(NULL)
              }
              tryCatch({
                  library(.Object@curOrgDb,character.only = TRUE)
              },
              error = function(e){
                  BiocManager::install(.Object@curOrgDb)
                  library(.Object@curOrgDb,character.only = TRUE)
              })
              #return(get0(.Object@curOrgDb))
              return(.Object@curOrgDb)
          })


#' @name configureValue
#' @aliases getConfigure
#' @aliases setConfigure
#' @aliases getAllConfigure
#' @title Global parameters configure
#' @description
#' These functions are used to configure and
#' query global parameters. The items include
#' "threads","tmpdir","refdir","genome","knownGene",
#' "bsgenome","bt2Idx","DHS" and "blacklist".
#' "threads","tmpdir","refdir","genome" are setable
#' and getable. While the others are readable only.
#' You should consider to configure these parameters
#' before starting the workflow.
#' 
#' For get \code{getConfigure} and \code{setConfigure}:
#' \describe{
#'   \item{"refdir"}{\code{Character} scalar, the path for reference data being installed to and storage.} 
#'   \item{"genome"}{\code{Character} scalar, the genome(like hg19, mm10, etc.) reference data in "refdir" to be used in the pipeline.} 
#'   \item{"tmpdir"}{\code{Character} scalar, the temporary file storage path} 
#'   \item{"threads"}{\code{Integer} scalar, the max threads allowed to be created} 
#' }
#' 
#' For get \code{getConfigure} only:
#' \describe{
#'   \item{"knownGene"}{\code{TxDb} scalar, known gene TxDb object} 
#'   \item{"bsgenome"}{\code{BSGenome} scalar, BSGenome object}
#'   \item{"bt2Idx"}{\code{Character} scalar, bowtie2 index path prefix} 
#'   \item{"DHS"}{\code{Character} scalar, DHS BED file path}
#'   \item{"blacklist"}{\code{Character} scalar, blacklist BED file path}  
#' }
#' @param item \code{Character} scalar.
#' The items that are setable or getable including
#' "threads","tmpdir","refdir","genome","knownGene",
#' "bsgenome","bt2Idx","DHS" and "blacklist".
#' @param val \code{Character} or \code{Integer} scalar.
#' The items value that are setable including
#' "threads","tmpdir","refdir","genome"
#' @return \code{Character} scalar for getting
#' "tmpdir","refdir", bt2Idx","DHS" and "blacklist",
#' all of them are file or directory path.
#' Getting "genome" will return the genome tag like "hg19" "mm10" etc.
#' \code{Integer} scalar for getting "threads",
#' the max threads can be created by the process.
#' \code{TxDb} object for getting "knownGene".
#' \code{BSGenome} object for getting "bsgenome"
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{atacBedUtils}}

#' @examples
#' getAllConfigure()
#'
#' getConfigure("threads")
#'
#' options(atacConf=setConfigure("tmpdir",tempdir()))
#'

#' @importFrom Rcpp  evalCpp
#' @importFrom igraph  graph
#' @importFrom igraph vertex.attributes
#' @importFrom igraph vertex.attributes<-
#' @importFrom igraph are.connected
#' @importFrom Rcpp sourceCpp
#' @importFrom rJava .jpackage
#' @importFrom rJava .jnew
#' @importFrom rJava .jcall
#' @importFrom rtracklayer import
#' @importFrom rtracklayer export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_ribbon
#' @importFrom DiagrammeR render_graph
#' @importFrom DiagrammeR to_igraph
#' @importFrom DiagrammeR select_nodes
#' @importFrom DiagrammeR trav_in
#' @importFrom DiagrammeR trav_out
#' @importFrom DiagrammeR set_node_attrs_ws
#' @importFrom DiagrammeR clear_selection
#' @importFrom DiagrammeR create_node_df
#' @importFrom DiagrammeR create_edge_df
#' @importFrom DiagrammeR create_graph
#' @importFrom DiagrammeR add_global_graph_attrs
#' @importFrom DiagrammeR export_graph
#' @importFrom magrittr %>%
#' @importFrom digest digest
#' @importFrom BSgenome getBSgenome
#' @importFrom Biostrings writeXStringSet
#' @importFrom GenomeInfoDb seqnames
#' @importFrom AnnotationDbi saveDb
#' @importFrom AnnotationDbi loadDb
#' @importFrom GenomicFeatures makeTxDbFromUCSC
#' @importFrom R.utils isGzipped
#' @importFrom R.utils gunzip
#' @importFrom R.utils isBzipped
#' @importFrom R.utils bunzip2
#' @importFrom GenomicRanges coverage
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocGenerics subset
#' @importFrom rmarkdown render
#' @importFrom knitr knit
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics plot
#' @importFrom stats binom.test
#' @importFrom stats fft
#' @importFrom stats runif
#' @importFrom utils download.file
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @import IRanges
#' @import S4Vectors
#' @import tools
#' @import Rsamtools
#' @import methods
# @importFrom markdown markdownToHTML
#' @useDynLib esATAC





#' @rdname configureValue
#' @export
getAllConfigure<-function(){
    .configObj <- getOption("atacConf")
    if(is.null(.configObj)){
        .configObj <- new(".ConfigClass")
    }
    .configObj@configList

}
#' @rdname configureValue
#' @export
getConfigure <- function(item = c("threads","tmpdir","refdir","genome","knownGene","bsgenome","annoDb","bt2Idx","DHS","blacklist","SNP")){
    .configObj <- getOption("atacConf")
    if(is.null(.configObj)){
        .configObj <- new(".ConfigClass")
    }
    isValidAttr(.configObj,item);
    if(item=="tmpdir"||item=="refdir"){
        if(is.null(.configObj@configList[[item]])){
            return(NULL)
        }else{
            return(normalizePath(.configObj@configList[[item]]))
        }
    }
    return(.configObj@configList[[item]]);
}


#' @rdname configureValue
#' @export
setConfigure<- function(item = c("threads","tmpdir","refdir","genome"),val){
    .configObj <- getOption("atacConf")
    if(is.null(.configObj)){
        .configObj <- new(".ConfigClass")
    }
    if(!is.null(val)){
        .configObj<-isValidVal(.configObj,item,val);
        .configObj@configList[[item]]<-val
        .configObj
    }else{
        warning("val is NULL, nothing to do")
    }
    .configObj
}

# @rdname configureValue
# @export
#.configObj<-new(".ConfigClass")


setAllConfigure<-function(threads=NULL,tmpdir=NULL,refdir=NULL,genome=NULL){
    setConfigure("threads",threads)
    setConfigure("tmpdir",tmpdir)
    setConfigure("refdir",refdir)
    setConfigure("genome",genome)
}

.obtainConfigure<-function(item = c("threads","tmpdir","refdir","genome","knownGene","bsgenome","annoDb","bt2Idx","DHS","blacklist","SNP")){
    .configObj <- getOption("atacConf")
    if(is.null(.configObj)){
        .configObj <- new(".ConfigClass")
    }
    isValidAttr(.configObj,item);
    if(item=="tmpdir"||item=="refdir"){
        val<-normalizePath(.configObj@configList[[item]])
    }else{
        val<-.configObj@configList[[item]]
    }
    if(is.null(val)){
        stop(paste(item,sprintf("has not been configured yet! Please call 'options(atacConf=setConfigure(\"%s\",\"set your value here\"))' to configure first",item)))
    }else{
        return(val)
    }
}

msgBoxBegin<-function(){
    message(">>>>>>========================================")
}

msgBoxDone<-function(){
    message("========================================<<<<<<")
    message(" ")
}






