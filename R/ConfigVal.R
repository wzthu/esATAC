.ConfigClass<-R6Class(classname = ".ConfigClass",
    public = list(
        initialize = function(){
            private$validAttr=list(threads="numeric",tmpdir="character",refdir="character",genome="character",knownGene="TxDb",bsgenome="BSgenome",bt2Idx="character",DHS="character",blacklist="character")
            private$validWriteAttr=list(threads="numeric",tmpdir="character",refdir="character",genome="character")
        },
        getAllConfigure = function(){
            print(private$configList)
        },
        getConfigure = function(item = c("threads","tmpdir","refdir","genome","knownGene","bsgenome","bt2Idx","DHS","blacklist")){
            private$isValidAttr(item);
            if(item=="tmpdir"||item=="refdir"){
                return(normalizePath(private$configList[[item]]))
            }
            return(private$configList[[item]]);
        },
        setConfigure = function(item = c("threads","tmpdir","refdir","genome"),val){
            private$isValidVal(item,val);
            private$configList[[item]]<-val;
        }
    ),
    private = list(
        configList=list(threads=detectCores(),tmpdir=".",refdir=NULL,genome=NULL,knownGene=NULL,bsgenome=NULL,bt2Idx=NULL,DHS=NULL,blacklist=NULL),
        validAttr=NULL,
        validWriteAttr=NULL,
        isValidAttr=function(item){
            if(is.null(private$validAttr[[item]])){
                stop(paste(item,"is not an attribute"))
            }
        },
        isWriteValidAttr=function(item){
            if(is.null(private$validWriteAttr[[item]])){
                stop(paste(item,"is not a writeable attribute"))
            }
        },
        isValidVal=function(item,val){
            private$isWriteValidAttr(item)
            if(is.null(val)){
                stop("val can not be NULL")
            }
            if(private$validAttr[[item]]!=class(val)){
                stop(paste(item,"is requied to be",private$validAttr[[item]],",\"",val,"\" is ",class(val)))
            }
            if(item=="refdir"||item=="tmpdir"){
                private$checkPathExist(val)
                val<-normalizePath(val)
            }
            if(item=="genome"){
                private$configList[["bsgenome"]]<-getBSgenome(val)
                if(is.null(private$configList[["refdir"]])){
                    stop("'refdir' should be configured before 'genome'")
                }
                fileprefix<-file.path(private$configList[["refdir"]],val)
                #genome fasta
                fastaFilePath<-paste0(fileprefix,".fa")
                fastaFilePathlock<-paste0(fileprefix,".fa.lock")
                if(file.exists(fastaFilePathlock)){
                    unlink(fastaFilePath)
                    unlink(fastaFilePathlock)
                }
                if(!file.exists(fastaFilePath)){
                    file.create(fastaFilePathlock)
                    private$BSgenomeSeqToFasta(private$configList[["bsgenome"]],fastaFilePath)
                    unlink(fastaFilePathlock)
                }

                #bowtie2 index
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
                if(!(file.exists(paste0(fileprefix,".1.bt2"))&&file.exists(paste0(fileprefix,".2.bt2"))||
                     file.exists(paste0(fileprefix,".3.bt2"))&&file.exists(paste0(fileprefix,".4.bt2"))||
                     file.exists(paste0(fileprefix,".rev.1.bt2"))&&file.exists(paste0(fileprefix,".rev.1.bt2")))){
                   file.create(fileprefixlock)
                    bowtie2_build(fastaFilePath,fileprefix,"--threads",as.character(.obtainConfigure("threads")),overwrite=TRUE)
                    unlink(fileprefixlock)
                }
                private$configList[["bt2Idx"]]<-fileprefix
                #kownGene
                knownGeneFilePath<-paste0(fileprefix,".knownGene.sqlite")
                knownGeneFilePathlock<-paste0(fileprefix,".knownGene.sqlite.lock")
                if(file.exists(knownGeneFilePathlock)){
                    unlink(knownGeneFilePath)
                    unlink(knownGeneFilePathlock)
                }
                if(!file.exists(knownGeneFilePath)){
                    file.create(knownGeneFilePathlock)
                    saveDb(makeTxDbFromUCSC(genome=val, tablename="knownGene"), file=knownGeneFilePath)
                    unlink(knownGeneFilePathlock)
                }
                private$configList[["knownGene"]]<-loadDb(knownGeneFilePath)
                if(sum(val==c("hg19","hg38","mm10","mm9"))){
                    blacklistFilePath<-paste0(fileprefix,".blacklist.bed")
                    DHSFilePath<-paste0(fileprefix,".DHS.bed")
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
                            download.file(url = sprintf("http://bioinfo.au.tsinghua.edu.cn/member/zwei/refdata/%s.DHS.bed",val),
                                          destfile = DHSFilePath,method = getOption("download.file.method"))
                        }
                        #blacklist
                        if(!file.exists(blacklistFilePath)){
                            download.file(url = sprintf("http://bioinfo.au.tsinghua.edu.cn/member/zwei/refdata/%s.blacklist.bed",val),
                                          destfile = blacklistFilePath,method = getOption("download.file.method"))
                        }
                        unlink(downloadFilePathlock)
                    }
                    private$configList[["DHS"]]<-DHSFilePath
                    private$configList[["blacklist"]]<-blacklistFilePath
                }


            }
        },
        checkFileExist = function(filePath){
          if(!is.null(filePath)){
            if(!file.exists(filePath)){
              stop(paste("error, file does not exist:",filePath))
            }
          }
        },
        checkPathExist = function(filePath){
          if(!is.null(filePath)){
            if(!dir.exists(dirname(filePath))){
              stop(paste("error, path does not exist:",filePath))
            }
          }
        },
        BSgenomeSeqToFasta = function(bsgenome, outFile)
        {
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
        },
        GetOrgDb = function(genome){
            if(genome == hg19||genome == hg38){
                return("org.Hs.eg.db")
            }
        }

    )
)

.configObj<-.ConfigClass$new()


#' @name configureValue
#' @aliases getConfigure
#' @aliases setConfigure
#' @title generate BigWig file from BED file
#' @description 
#' This function is used to generate BigWig file 
#' from BED reads file.
#' The BigWig file can be shown reads coverage on genome browser.
#' @param atacProc \code{\link{ATACProc}} object scalar. 
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}, 
#' \code{\link{atacBedUtils}}.
#' @param item \code{Character} scalar. 
#' Bed file input path. 
#' @param val \code{Character} scalar. 
#' BigWig file output path.
#' @return \code{Character} scalar or TxDb object.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacSamToBed}} 
#' \code{\link{atacBedUtils}}


getAllConfigure<-function(){
    .configObj$getAllConfigure();
}
#' @rdname configureValue
#' @export 
getConfigure <- function(item = c("threads","tmpdir","refdir","genome","knownGene","bsgenome","bt2Idx","DHS","blacklist")){
    return(.configObj$getConfigure(item));
}


#' @rdname configureValue
#' @export 
setConfigure<- function(item = c("threads","tmpdir","refdir","genome"),val){
    if(is.null(val)){
        return()
    }
    .configObj$setConfigure(item,val);
}

setAllConfigure<-function(threads=NULL,tmpdir=NULL,refdir=NULL,genome=NULL){
    setConfigure("threads",threads)
    setConfigure("tmpdir",tmpdir)
    setConfigure("refdir",refdir)
    setConfigure("genome",genome)
}

.obtainConfigure<-function(item = c("threads","tmpdir","refdir","genome","knownGene","bsgenome","bt2Idx","DHS","blacklist")){
    val<-.configObj$getConfigure(item);
    if(is.null(val)){
        stop(paste(item,"has not been configured yet! Please call 'setConfigure' to configure first"))
    }else{
        return(val)
    }
}









