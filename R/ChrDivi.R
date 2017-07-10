ChrDivi <- R6::R6Class(
  classname = "ChrDivi",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, ReadsIfile = NULL, ReadsOpath = NULL, editable = FALSE){
      super$initialize("ChrDivi",editable,list(arg1=atacProc))
      print("ChrDiviInitCall")
      # add "/" in the output path
      if(substr(ReadsOpath, nchar(ReadsOpath), nchar(ReadsOpath)) != "/"){
        ReadsOpath <- paste0(ReadsOpath, "/", collapse = "")
      }
      # necessary and unchanged parameters, this should be tested
      private$paramlist[["ReadsIfile"]] <- ReadsIfile
      private$paramlist[["ReadsOpath"]] <- ReadsOpath
      # parameter check
      private$checkRequireParam()
      private$checkFileExist(private$paramlist[["ReadsIfile"]])
      if(!dir.exists(private$paramlist[["ReadsOpath"]])){
        stop(paste("error, path does not exist:",private$paramlist[["ReadsOpath"]]))
      }
      print("finishChrDiviInitCall")
    },

    processing = function(){
      super$processing()
      .chr_separate_call(ReadsIfile = private$paramlist[["ReadsIfile"]], ReadsOpath = private$paramlist[["ReadsOpath"]])
      private$finish <- TRUE
    },

    setResultParam = function(){
      print("This function is not using!")
    }

  ), #public end

  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["ReadsIfile"]])){
        stop("Parameter ReadsIfile is required!")
      }
      if(is.null(private$paramlist[["ReadsOpath"]])){
        stop("Parameter ReadsOpath is required!")
      }
    }
  ) # private end

)
