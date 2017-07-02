RemoveAdapter <-R6Class(
  classname = "RemoveAdapter",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc,adapter1=NULL,adapter2=NULL,fastqOutput1=NULL,reportPrefix=NULL,
                          fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL,editable=FALSE){
      super$initialize("RemoveAdapter",editable,list(arg1=atacProc))
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
        private$paramlist[["fastqOutput1"]] <- paste(private$paramlist[["fastqInput1"]],".clipped",sep="");
      }else{
        private$paramlist[["fastqOutput1"]] <- fastqOutput1;
      }
      if(is.null(fastqOutput2)){
        private$paramlist[["fastqOutput2"]] <- paste(private$paramlist[["fastqInput2"]],".clipped",sep="");
      }else{
        private$paramlist[["fastqOutput2"]] <- fastqOutput2;
      }
      if(is.null(reportPrefix)){
        private$paramlist[["reportPrefix"]] <- paste(private$paramlist[["fastqInput1"]],".report",sep="");
      }else{
        private$paramlist[["reportPrefix"]] <- reportPrefix;
      }
      
      private$paramlist[["adapter1"]] <- adapter1
      private$paramlist[["adapter2"]] <- adapter2
      
      
      private$checkFileExist(private$paramlist[["fastqInput1"]]);
      private$checkFileExist(private$paramlist[["fastqInput2"]]);
      private$checkPathExist(private$paramlist[["fastqOutput1"]]);
      private$checkPathExist(private$paramlist[["fastqOutput2"]]);
      private$checkPathExist(private$paramlist[["reportPrefix"]]);
      private$checkRequireParam();
    },
    processing = function(){
      super$processing()
      if(is.null(private$paramlist[["adapter1"]])){
        adapters<-.identify_adapters_call(private$paramlist[["fastqInput1"]],
                                         private$paramlist[["fastqInput2"]],threads=getConfigure("threads"))
        private$paramlist[["adapter1"]] <- adapters$adapter1
        private$paramlist[["adapter2"]] <- adapters$adapter2
      }
      .remove_adapters_call(inputFile1=private$paramlist[["fastqInput1"]],adapter1=private$paramlist[["adapter1"]],
                           outputFile1 = private$paramlist[["fastqOutput1"]],basename = private$paramlist[["reportPrefix"]],
                           inputFile2=private$paramlist[["fastqInput2"]],adapter2=private$paramlist[["adapter2"]],
                           outputFile2 = private$paramlist[["fastqOutput2"]],threads=getConfigure("threads"))
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
      if(is.null(private$paramlist[["fastqInput2"]])&&is.null(adapter1)){
        stop("Parameter \"adapter1\" is requied for single end sequencing.")
      }
      
    }
  )
  

)

