Renamer <-R6Class(
  classname = "Renamer",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, fastqOutput1=NULL, fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL,editable=FALSE){
      super$initialize("Renamer",editable,list(arg1=atacProc))
      print("RenamerInitCall")
      if(!is.null(atacProc)){
        private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1");
        private$paramlist[["fastqInput2"]] <- atacProc$getParam("fastqOutput2");
      }
      
      if(!is.null(fastqInput1)){
        private$paramlist[["fastqInput1"]] <- fastqInput1;
      }
      if(!is.null(fastqInput2)){
        private$paramlist[["fastqInput2"]] <- fastqInput2;
      }
      
      if(is.null(fastqOutput1)){
        private$paramlist[["fastqOutput1"]] <- paste(private$paramlist[["fastqInput1"]],".renamed",sep="");
      }else{
        private$paramlist[["fastqOutput1"]] <- fastqOutput1;
      }
      if(is.null(fastqOutput2)){
        if(!is.null(private$paramlist[["fastqInput2"]])){
          private$paramlist[["fastqOutput2"]] <- paste(private$paramlist[["fastqInput2"]],".renamed",sep="");
        }
      }else{
        private$paramlist[["fastqOutput2"]] <- fastqOutput2;
      }
      
      private$checkFileExist(private$paramlist[["fastqInput1"]]);
      private$checkFileExist(private$paramlist[["fastqInput2"]]);
      
      private$checkPathExist(private$paramlist[["fastqOutput1"]]);
      private$checkPathExist(private$paramlist[["fastqOutput2"]]);
      
      private$checkRequireParam()
      print("finishRenamerInitCall")
    },
    processing = function(){
      super$processing()
      .renamer_call(private$paramlist[["fastqInput1"]],private$paramlist[["fastqOutput1"]])
      if(!is.null(private$paramlist[["fastqInput2"]])){
        .renamer_call(private$paramlist[["fastqInput2"]],private$paramlist[["fastqOutput2"]])
      }
      private$finish<-TRUE
    },
    setResultParam = function(fastqOutput1, fastqOutput2=NULL){
      super$setResultParam();
      private$paramlist[["fastqOutput1"]] <- fastqOutput1
      private$paramlist[["fastqOutput2"]] <- fastqOutput2
    }
  ),
  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["fastqInput1"]])){
        stop("fastqInput1 is required.")
      }
      
    }
  )

)
