RPeakAnno <- R6::R6Class(
  classname = "RPeakAnno",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, Input = NULL, Output = NULL, editable = FALSE){
      super$initialize("RPeakAnno",editable,list(arg1=atacProc))
      print("RPeakAnnoInitCall")

      # necessary and unchanged parameters, this should be tested
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }

      private$paramlist[["Input"]] <- Input
      private$checkRequireParam()

      if(is.null(Output)){
        private$paramlist[["Output"]] <- paste0(dirname(Input), "/output", collapse = "")
      }else{
        private$paramlist[["Output"]] <- Output
      }

      private$checkFileExist(private$paramlist[["Input"]])
      private$checkPathExist(private$paramlist[["Output"]])
      print("finishRPeakAnnoInitCall")

    }, # initialization end

    processing = function(){
      if(!super$processing()){
        return()
      }
      peak <- ChIPseeker::readPeakFile(private$paramlist[["Input"]])
      peakAn <- ChIPseeker::annotatePeak(peak)
      write.table(x = peakAn, file = private$paramlist[["Output"]], append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
      private$setFinish()
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
        stop("Parameter Input is required!")
      }
    }
  ) # private end

) # R6 class end
