SamToBed <- R6::R6Class(
  classname = "SamToBed",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, samfile = NULL, bedfile = NULL, editable=FALSE){
      super$initialize("SamToBed",editable,list(arg1=atacProc))
      print("SamToBedInitCall")

      if(!is.null(atacProc)){
        private$paramlist[["SamInput"]] <- atacProc$getParam("samOutput");
        if(is.null(bedfile)){
          private$paramlist[["BedOutput"]] <- paste0(private$paramlist[["SamInput"]], ".bed",collapse = "")
        }else{
          private$paramlist[["BedOutput"]] <- bedfile
        }
        # file check
        private$checkFileExist(private$paramlist[["SamInput"]]);
        private$checkPathExist(private$paramlist[["BedOutput"]]);
      }else{
        private$paramlist[["SamInput"]] <- samfile;
        if(is.null(bedfile)){
          private$paramlist[["BedOutput"]] <- paste0(private$paramlist[["SamInput"]], ".bed",collapse = "")
        }else{
          private$paramlist[["BedOutput"]] <- bedfile
        }
        # file check
        private$checkFileExist(private$paramlist[["SamInput"]]);
        private$checkPathExist(private$paramlist[["BedOutput"]]);
      }
      private$checkRequireParam()
      print("finishMappingInitCall")
    },

    processing = function(){
      super$processing()
      .sam2bed_call(samfile = private$paramlist[["SamInput"]], bedfile = private$paramlist[["BedOutput"]])
      private$finish <- TRUE
    },

    setResultParam = function(bedFilePath){
      super$setResultParam();
      private$paramlist[["BedOutput"]] <- bedFilePath
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
