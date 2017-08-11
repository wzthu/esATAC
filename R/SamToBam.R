SamToBam <-R6::R6Class(
  classname = "SamToBam",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, samfile = NULL, bamfile = NULL, editable=FALSE){
      super$initialize("SamToBam",editable,list(arg1=atacProc))
      print("SamToBamInitCall")
      private$checkRequireParam()
#      if((!is.null(atacProc)) && (class(atacProc)[1] == "Mapping")){ # atacproc from mapping
      print(atacProc$getParam("samOutput"))
      if((!is.null(atacProc)) ){
        private$paramlist[["samInput"]] <- atacProc$getParam("samOutput")
#      }else if((!is.null(atacProc)) && (class(atacProc)[1] != "Mapping")){ # atacproc not from mapping, error!
#        stop("Input class must be got from 'Mapping'!")
      }else if(is.null(atacProc)){ # input
        private$paramlist[["samInput"]] <- samfile
      }
      if(is.null(bamfile)){
        private$paramlist[["bamOutput"]] <- private$paramlist[["samInput"]]
      }else{
        private$paramlist[["bamOutput"]] <- bamfile
      }
      # file check
      private$checkFileExist(private$paramlist[["samInput"]])
      private$checkPathExist(private$paramlist[["bamOutput"]])

      print("finishSamToBamInitCall")
    },


    processing = function(){
      if(!super$processing()){
        return()
      }
      Rsamtools::asBam(file = private$paramlist[["samInput"]], destination = private$paramlist[["bamOutput"]])
      private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["samInput"]], ".bam", collapse = "")
      private$paramlist[["baiOutput"]] <- paste0(private$paramlist[["samInput"]], ".bam.bai", collapse = "")
      private$setFinish()
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
