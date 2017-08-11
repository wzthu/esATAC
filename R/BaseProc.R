BaseProc <- R6Class(
  classname = "BaseProc",
  public = list(
    fileType = list(fq="fq",fastq="fq",fa="fa",fasta="fa"),
    initialize = function(procName,editable,atacProcs){
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
    processing = function(){
        if(private$editable){
            stop("The \"processing\" method can not be call in editable result object");
        }
        if(private$checkMD5Cache()){
            message(paste0("This process:`",private$procName,"` was finished. Nothing to be done."))
            message("If you need to redo, please call 'YourObject$clearCache()'")
            private$finish<-TRUE
            return(FALSE)
        }else{
            return(TRUE)
        }
    },

    getNextProcList = function(){
      return(private$graphMng$getNextList())
    },
    getProcName = function(){
      return(private$procName)
    },
    printMap = function(preProc=FALSE,nextProc=TRUE,curProc=TRUE){
      private$graphMng$printMap(procName=private$procName,preProc=preProc,nextProc=nextProc,curProc=curProc)
    },
    finalize = function(){
      #rm(private$graphMng)
      #rm(private$paramlist)
      #rm(private$procName)
    },
    getParam = function(item){
      return(private$paramlist[[item]])
    },
    setResultParam = function(){
      if(!private$editable){
        stop("This object can not be edited");
      }
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
            unlink(private$paramlist[rslist[i]])
        }
    },
    isSingleEnd = function(){
        return(private$singleEnd)
    }
  ),
  private = list(
    paramlist = list(),
    procName = "",
    completObj = TRUE,
    editable = FALSE,
    finish = FALSE,
    graphMng = NULL,
    singleEnd = FALSE,
    checkRequireParam = function(){
      stop("processing function has not been implemented")
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
                warning(paste("file exist:",filePath,", we will overwrite the file"));
            }else if(!file.create(filePath)){
                stop(paste("cannot create file '",filePath,"', No such file or directory"));
            }
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
        md5code<-digest(object = paramstr,algo = "md5")
        curtmpdir<-.obtainConfigure("tmpdir")
        md5filepath<-file.path(curtmpdir,paste(private$procName,md5code,"log",sep = "."))
        return(md5filepath)
    },
    setFinish = function(){
        private$finish<-TRUE
        file.create(private$getParamMD5Path())
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
    }



  )

)


