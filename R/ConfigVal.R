.ConfigClass<-R6Class(classname = ".ConfigClass",
    public = list(
        initialize = function(){
            private$validAttr=list(threads="numeric",tmpdir="character",datadir="character",genome="character",knownGene="TxDb",bsgenome="BSgenome",bt2Idx="character")
            private$validWriteAttr=list(threads="numeric",tmpdir="character",datadir="character",genome="character")
        },
        getAllConfigure = function(){
            print(private$configList)
        },
        getConfigure = function(item = c("threads","tmpdir","datadir","genome","knownGene","bsgenome","bt2Idx")){
            private$isValidAttr(item);
            if(item=="tmpdir"||item=="datadir"){
                return(normalizePath(private$configList[[item]]))
            }
            return(private$configList[[item]]);
        },
        setConfigure = function(item = c("threads","tmpdir","datadir","genome"),val){
            private$isValidVal(item,val);
            private$configList[[item]]<-val;
        }
    ),
    private = list(
        configList=list(threads=detectCores(),tmpdir=".",datadir=NULL,genome=NULL,knownGene=NULL,bsgenome=NULL,bt2Idx=NULL),
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
            if(item=="datadir"||item=="tmpdir"){
                private$checkPathExist(val)
                val<-normalizePath(val)
            }
            if(item=="genome"){
                private$configList[["bsgenome"]]<-getBSgenome(val)
                if(is.null(private$configList[["datadir"]])){
                    stop("'datadir' should be configured before 'genome'")
                }
                fileprefix<-file.path(private$configList[["datadir"]],val)
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
                    .bowtie2_build_call(fastaFilePath,fileprefix,c("--threads",as.character(.obtainConfigure("threads"))))
                    #buildBowtie2Index(fastaFilePath)################################赶紧实现
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
        }

    )
)

.configObj<-.ConfigClass$new()

getAllConfigure<-function(){
    .configObj$getAllConfigure();
}

getConfigure <- function(item = c("threads","tmpdir","datadir","genome","kownGene","bsgenome","bt2Idx")){
    return(.configObj$getConfigure(item));
}



setConfigure<- function(item = c("threads","tmpdir","datadir","genome"),val){
    if(is.null(val)){
        return()
    }
    .configObj$setConfigure(item,val);
}

setAllConfigure<-function(threads=NULL,tmpdir=NULL,datadir=NULL,genome=NULL){
    setConfigure("threads",threads)
    setConfigure("tmpdir",tmpdir)
    setConfigure("datadir",datadir)
    setConfigure("genome",genome)
}

.obtainConfigure<-function(item = c("threads","tmpdir","datadir","genome","kownGene","bsgenome","bt2Idx")){
    val<-.configObj$getConfigure(item);
    if(is.null(val)){
        stop(paste(item,"has not been configured yet! Please call 'setConfigure' to configure first"))
    }else{
        return(val)
    }
}









