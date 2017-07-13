CutSitePre <- R6::R6Class(
  classname = "CutSitePre",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, InputFile = NULL, OutputFile = NULL, prefix = NULL, editable = FALSE){
      super$initialize("CutSitePre",editable,list(arg1=atacProc))
      print("CutSitePreInitCall")

      # necessary and unchanged parameters, this should be tested
      private$paramlist[["Input"]] <- InputFile
      if(is.null(prefix)){
        private$paramlist[["Output"]] <- paste0(OutputFile, "/", "output_", collapse = "")
      }else{
        private$paramlist[["Output"]] <- paste0(OutputFile, "/", prefix, "_", collapse = "")
      }
      private$checkRequireParam()
      private$checkFileExist(private$paramlist[["Input"]])
      private$checkPathExist(private$paramlist[["Output"]])
      print("finishCutSitePreInitCall")
    }, # initialization end

    processing = function(){
      super$processing()
      .CutSite_call(InputFile = private$paramlist[["Input"]], OutputFile = private$paramlist[["Output"]])
      private$finish <- TRUE
    }, # processing end

    setResultParam = function(){
      print("This function is not using!")
    }  # setResultParam end

  ), # public end


  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["Input"]])){
        stop("Parameter InputFile is required!")
      }
      if(is.null(private$paramlist[["Output"]])){
        stop("Parameter OutputFile is required!")
      }
    }
  ) # private end

) # class end
