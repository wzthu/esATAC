.ConfigClass<-R6Class(classname = ".ConfigClass",
  public = list(
    getAllConfigure = function(){
      print(private$configList)
    },
    getConfigure = function(item = c("threads","tmpdir","datadir","genome","knownGene")){
      private$isValidAttr(item);
      return(private$configList[[item]]);
    },
    setConfigure = function(item = c("threads","tmpdir","datadir","genome","knownGene"),val){
      private$isValidVal(item,val);
      private$configList[[item]]<-val;
    }

  ),
  private = list(
    configList=list(threads=1,tmpdir=NULL),
    validAttr=list(threads="numeric",tmpdir="character",datadir="character",genome="character",knownGene="TxDb"),
    isValidAttr=function(item){
      if(is.null(private$validAttr[[item]])){
        stop(paste(item,"is not a attribute"))
      }
    },
    isValidVal=function(item,val){
      private$isValidAttr(item)
      if(!is.null(val)&&private$validAttr[[item]]!=class(val)){
        stop(paste(item,"is requied to be",private$validAttr[[item]],",\"",val,"\" is ",mode(val)))
      }
      if(item=="thread"&&is.null(val)){
          stop("thread can not be NULL");
      }
      if(item=="datadir"){
        private$checkPathExist(val)
      }
      if(item=="genome"){
        validgenome <- c("hg19","hg38","mm9","mm10")
        if(sum(validgenome==val)<1){
          stop(paste(val,"is invalid genome type, only",validgenome,"are supported"));
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
    }

  )
)

.configObj<-.ConfigClass$new()

getAllConfigure<-function(){
  .configObj$getAllConfigure();
}

getConfigure <- function(item = c("threads","tmpdir","datadir","genome","knownGene")){
  return(.configObj$getConfigure(item));
}

setConfigure<- function(item = c("threads","tmpdir","datadir","genome","knownGene"),val){
  .configObj$setConfigure(item,val);
}

.obtainConfigure<-function(item = c("threads","tmpdir","datadir","genome","knownGene")){
    val<-.configObj$getConfigure(item);
    if(is.null(val)){
        stop(paste(item,"has not been configured yet! Please call 'setConfigure' to configure first"))
    }else{
        return(val)
    }
}



