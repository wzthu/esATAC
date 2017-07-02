Rsortbam <- R6::R6Class(
  classname = "Rsortbam",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, inputbam = NULL, outputbam = NULL, editable=FALSE){
      super$initialize("Rsortbam",editable,list(arg1=atacProc))
      print("RsortbamInitCall")
      private$checkRequireParam()
      if((!is.null(atacProc)) && (class(atacProc)[1] == "SamToBam")){ # atacproc from SamToBam
        private$paramlist[["bamInput"]] <- atacProc$getParam("bamOutput")
      }else if((!is.null(atacProc)) && (class(atacProc)[1] != "SamToBam")){ # atacproc not from SamToBam, error!
        stop("Input class must be got from 'SamToBam'!")
      }else if(is.null(atacProc)){ # input
        private$paramlist[["bamInput"]] <- inputbam
      }
      if(is.null(outputbam)){
        private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["bamInput"]], ".sorted", collapse = "")
      }else{
        private$paramlist[["bamOutput"]] <- outputbam
      }
      # file check
      private$checkFileExist(private$paramlist[["bamInput"]])
      private$checkPathExist(private$paramlist[["bamOutput"]])
      print("finishRsortbamInitCall")

    },

    processing = function(){
      super$processing()
      Rsamtools::sortBam(file = private$paramlist[["bamInput"]], destination = private$paramlist[["bamOutput"]])
      private$finish <- TRUE
    },

    setResultParam = function(bamFilePath){
      super$setResultParam();
      private$paramlist[["bamOutput"]] <- bamFilePath
    }

  ),

  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
    }
  )


)
