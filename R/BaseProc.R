ATACProc <- R6Class(
  classname = "ATACProc",
  public = list(
    fileType = list(fq="fq",fastq="fq",fa="fa",fasta="fa"),
    initialize = function(procName,editable,atacProcs){
      private$timeStamp <- format(Sys.time(),"[%H%M%S]")
      private$graphMng<-GraphMng$new()
      private$procName<-procName
      private$editable<-editable
      if(private$editable){
        private$finish <- TRUE;
      }
      argSize <- length(atacProcs)
      if(argSize>=1&&!is.null(atacProcs[[1]])){
        if(!atacProcs[[1]]$isReady()){
          stop(paste(atacProcs[[1]]$getProcName(),"is not ready"))
        }
        if(!private$graphMng$checkRelation1(atacProcs[[1]]$getProcName(),procName)){
          stop(paste(atacProcs[[1]]$getProcName(),"is not valid input"))
        }
        private$singleEnd<-atacProcs[[1]]$isSingleEnd()
      }else if(argSize>=2&&!is.null(atacProcs[[2]])){
          if(!atacProcs[[2]]$isReady()){
              stop(paste(atacProcs[[2]]$getProcName(),"is not ready"))
          }
          if(!private$graphMng$checkRelation2(atacProcs[[2]]$getProcName(),procName)){
              stop(paste(atacProcs[[2]]$getProcName(),"is not valid input"))
          }
      }
    },
    process = function(){
        if(private$editable){
            stop("The \"processing\" method can not be call in editable result object");
        }
        if(private$checkMD5Cache()){
            message(paste0("The process:`",private$procName,"` was finished. Nothing to do."))
            message("If you need to redo, please call 'clearProcCache(YourObject)' or 'YourObject$clearCache()'")
            private$finish<-TRUE
            return(FALSE)
        }else{
            private$writeLog(as.character(Sys.time()))
            if(private$singleEnd){
                private$writeLog(paste0("start processing(single end data): ", private$procName))
            }else{
                private$writeLog(paste0("start processing(paired end data): ", private$procName))
            }
            #private$paramlistbk<-private$paramlist
            private$processing()
            private$setFinish()
            return(TRUE)
        }
    },
    

    getNextProcList = function(){
      return(private$graphMng$getNextList())
    },
    getProcName = function(){
      return(private$procName)
    },
    printMap = function(preProc=FALSE,nextProc=TRUE,curProc=TRUE,display=TRUE){
      private$graphMng$printMap(procName=private$procName,preProc=preProc,nextProc=nextProc,curProc=curProc,display=display)
    },
    getParam = function(item){
      return(private$paramlist[[item]])
    },
    getParamItems = function(){
        return(names(private$paramlist))
    },
    setResultParam = function(item,val){
      if(!private$editable){
        stop("This object can not be edited");
      }
      private$paramlist[[item]]<-val
    },
    isReady = function(){
      if(private$editable){
        return(TRUE);
      }else if(private$finish){
        return(TRUE);
      }else{
        return(FALSE);
      }

    },
    clearCache = function(){
        if(!unlink(private$getParamMD5Path())){
            message("Chache has been cleared")
        }else{
            message("Chache does not exist. Nothing has been done.")
        }
        rslist<-grep("(o|O)utput",names(private$paramlist))
        for(i in 1:length(rslist)){
            unlink(private$paramlist[[rslist[i]]])
        }
        private$finish<-FALSE
    },
    isSingleEnd = function(){
        return(private$singleEnd)
    },
    getReportVal = function(item){
        if(private$finish){
            if(sum(item == private$getReportItemsImp())>0){
                return(private$getReportValImp(item))
            }else{
                stop(paste0(item," is not an item of report value.")) 
            }
        }else{
            stop("Unfinished process is not available for report value.")
        }
        
    },
    getReportItems = function(){
        if(private$finish){
            return(private$getReportItemsImp())
        }else{
            stop("Unfinished process is not available for report items.")
        }
    }
  ),
  private = list(
    paramlist = list(),
    paramlistbk = list(),
    reportVal = list(),
    procName = "",
    completObj = TRUE,
    editable = FALSE,
    finish = FALSE,
    graphMng = NULL,
    singleEnd = FALSE,
    logRecord= NULL,
    timeStamp="",
    timeStampPattern="(\\[\\d\\d\\d\\d\\d\\d\\])*",
    getAutoPath = function(originPath,regexProcName,suffix){
        if(!is.null(originPath)){
            prefix<-private$getBasenamePrefix(originPath,regexProcName)
            return(file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),suffix)))
        }else{
            return(NULL)
        }
    },
    paramValidation = function(){
        if(!private$checkMD5Cache()){
            private$checkAllPath()
        }
        if(!private$editable){
            private$checkRequireParam();
        }
    },
    checkRequireParam = function(){
      stop("checkRequireParam function has not been implemented")
    },
    checkAllPath = function(){
      stop("checkAllPath function has not been implemented")
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
    checkFileCreatable = function(filePath){
        if(!is.null(filePath)){
            if(file.exists(filePath)){
                private$writeLog(paste("file exist:",filePath,". It may be overwrited in processing"));
            }else if(!file.create(filePath)){
                stop(paste("cannot create file '",filePath,"', No such file or directory or permission denied"));
            }else{
                unlink(filePath)
            }
        }
    },
    checkParam = function(paramlist,paramPattern){
        rs<-grepl(paramPattern, paramlist)
        if(sum(rs)>0){
            banp=paste(paramlist[rs], collapse = "'/'")
            stop(sprintf("Parameter(s) '%s' are not acceptable in paramList. it should be set as fix parameter.",banp))
        }

    },
    getSuffix = function(filePath){
        filename<-basename(filePath)
        lst=strsplit(filename,"\\.")[[1]]
        if(length(lst)==1){
            return(NULL)
        }else{
            return(lst[length(lst)])
        }
    },
    getSuflessFileName = function(filePath){
        sfx=private$getSuffix(filePath)
        if(is.null(sfx)){
            return(filePath)
        }else {
            return(strsplit(filePath,paste0(".",sfx)))
        }
    },
    getParamMD5Path = function(){
        paramstr=c(private$procName)
        for(n in sort(names(private$paramlist))){
            paramstr<-c(paramstr,n)
            paramstr<-c(paramstr,private$paramlist[[n]])
        }
        rslist<-grep("((o|O)utput)|((i|I)nput)",names(private$paramlist))
        flag = FALSE
        for(i in 1:length(rslist)){
            filelist<-private$paramlist[[rslist[i]]]
            for(j in 1:length(filelist)){
                if(!is.character(filelist[j])){
                    flag = TRUE
                    break
                }
                fileinfo<-file.info(filelist[j])
                if(is.na(fileinfo$isdir)){
                    paramstr<-c(paramstr,runif(1))
                    flag = TRUE
                    break
                }
                if(!fileinfo$isdir){
                    paramstr<-c(paramstr,fileinfo$size)
                }
            }
            if(flag){
                break
            }
        }
        md5code<-substr(digest(object = paramstr,algo = "md5"),1,8)
        curtmpdir<-.obtainConfigure("tmpdir")
        md5filepath<-file.path(curtmpdir,paste(private$procName,md5code,"log",sep = "."))
        return(md5filepath)
    },
    setFinish = function(){
        private$finish<-TRUE
        private$writeLog(as.character(Sys.time()))
        private$writeLog("processing finished")
        logFilePath<-private$getParamMD5Path()
        write.table(private$logRecord,logFilePath,quote = FALSE,row.names = FALSE,col.names = FALSE)
    },
    checkMD5Cache = function(){
        if(file.exists(private$getParamMD5Path())){
            return(TRUE)
        }else{
            return(FALSE)
        }
    },
    getBasenamePrefix = function(filepath,words){
        return(basename(gsub(paste0("[.]",words,".*"),"",filepath)))
    },
    getPathPrefix = function(filepath,words){
        return(gsub(paste0("[.]",words,".*"),"",filepath))
    },
    writeLog = function(msg,isWarnning=FALSE,appendLog=TRUE){
        if(isWarnning){
            warning(msg)
            msg<-paste0("Warning:",msg)
        }else{
            message(msg)
        }
        if(appendLog){
            private$logRecord<-c(private$logRecord,msg)
        }else{
            private$logRecord<-msg
        }

    },
    processing = function(){
        stop("processing function has not been implemented")
    },
    getReportValImp = function(item){
        return(private$reportVal[[item]])
    },
    getReportItemsImp = function(){
        return(names(private$reportVal))
    }
    



  )

)
#' @name ATACProc
#' @aliases printMap
#' @aliases printNextList
#' @aliases printPrevList
#' @aliases clearProcCache
#' @aliases process
#' @aliases getNextProcList
#' @aliases getProcName
#' @aliases getParam
#' @aliases getParamItems
#' @aliases isReady
#' @aliases isSingleEnd
#' @aliases getReportVal
#' @aliases getReportItems
#' @title Methods for ATACProc objects
#' @description 
#' You can call ATACProc objects operation methods below to 
#' obtain information in objects. Or, you call also call
#' "atacProc$method(parameters)" to do the same things.
#' @details 
#' ATACProc is a R6ClassGenerator for generating ATACProc R6 objects. 
#' All ATACProc objects generated by its subclasses R6ClassGenerator.
#' You can only use the ATACProc objects returned by any functions 
#' rather than use ATACProc R6ClassGenerator to generate yourself.
#' @return The value return by object methods.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacSamToBed}} 
#' \code{\link{atacBedUtils}}

#' @param atacProc ATACProc object return by functions or Character scalar.
#' @param preProc \code{Logitcal} scalar. 
#' show the available upstream processes if TRUE
#' @param nextProc \code{Logitcal} scalar.
#' show the available downstream processes if TRUE
#' @param curProc \code{Logitcal} scalar.
#' show the current process of parameter \code{atacProc} if TRUE
#' @param display \code{Logitcal} scalar.
#' Save to pdf file if FALSE.
#' @param item \code{Characters} scalar
#' The parameters name
#' @return the function and result of functions 
#' @rdname ATACProc
#' @export 
printMap <-function(atacProc=NULL,preProc=FALSE,nextProc=TRUE,curProc=TRUE,display=TRUE){
    if(is.null(atacProc)){
        GraphMng$new()$printMap(display=display)
    }else if(is.character(atacProc)){
        GraphMng$new()$printMap(atacProc,preProc,nextProc,curProc,display=display)
    }else{
        atacProc$printMap(preProc,nextProc,curProc,display=display)
        invisible(atacProc)
    }
}
#' @rdname ATACProc
#' @return \item{printNextList}{Print the map to show available downstream process}
#' @export 
printNextList<-function(atacProc){
    if(is.character(atacProc)){
        GraphMng$new()$getNextProcs(atacProc)
    }else{
        GraphMng$new()$getNextProcs(atacProc$getProcName())
        invisible(atacProc)
    }
}
#' @rdname ATACProc
#' @return \item{printPrevList}{Print the map to show available upstream process}
#' @export 
printPrevList<-function(atacProc){
    if(is.character(atacProc)){
        GraphMng$new()$getPrevProcs(atacProc)
    }else{
        GraphMng$new()$getPrevProcs(atacProc$getProcName())
        invisible(atacProc)
    }
}
#' @rdname ATACProc
#' @return \item{clearProcCache}{Clear cache of atacProc object}
#' @export 
clearProcCache <- function(atacProc){
    invisible(atacProc$clearCache())
}
#' @rdname ATACProc
#' @return \item{process}{Call this function to redo processing }
#' @export 
process <- function(atacProc){
    invisible(atacProc$process())
}

#' @rdname ATACProc
#' @return \item{getNextProcList}{Get list of available downstream process}
#' @export 
getNextProcList <- function(atacProc){
    return(atacProc$getNextProcList())
}


#' @rdname ATACProc
#' @return \item{getProcName}{get atacProc Characher name}
#' @export 
getProcName = function(atacProc){
    return(atacProc$getProcName())
}
#' @rdname ATACProc
#' @return \item{getParam}{Get parameter value setted by process function}
#' @export 
getParam = function(atacProc,item){
    return(atacProc$getParam(item))
}
#' @rdname ATACProc
#' @return \item{getParamItems}{Get parameter name list}
#' @export 
getParamItems = function(atacProc){
    return(atacProc$getParamItems())
}
#' @rdname ATACProc
#' @return \item{isReady}{Is the process ready}
#' @export 
isReady = function(atacProc){
    return(atacProc$isReady())
    
}
#' @rdname ATACProc
#' @return \item{isSingleEnd}{Single end data if TRUE else FALSE}
#' @export 
isSingleEnd = function(atacProc){
    return(atacProc$isSingleEnd())
}
#' @rdname ATACProc
#' @return \item{getReportVal}{Get report value of item}
#' @export 
getReportVal = function(atacProc,item){
    return(atacProc$getReportVal(item))
    
}
#' @rdname ATACProc
#' @return \item{getReportItems}{Get all items that can be reported}
#' @export 
getReportItems = function(atacProc){
    return(atacProc$getReportItems())
}