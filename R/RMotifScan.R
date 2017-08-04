RMotifScan <- R6::R6Class(
  classname = "RMotifScan",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, Seq_file = NULL, Motif = NULL,
                          output_file = NULL, editable = FALSE){
      super$initialize("RMotifScan",editable,list(arg1=atacProc))
      print("RMotifScanInitCall")

      # necessary and unchanged parameters, this should be tested
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }
      private$paramlist[["Seq_file"]] <- Seq_file
      private$paramlist[["Motif"]] <- Motif
      private$checkRequireParam()

      if(is.null(output_file)){
        private$paramlist[["output_file"]] <- paste0(dirname(Seq_file), "/output", collapse = "")
      }else{
        private$paramlist[["output_file"]] <- output_file
      }

      private$checkFileExist(private$paramlist[["Seq_file"]])
      private$checkPathExist(private$paramlist[["output_file"]])
      print("finishRMotifScanInitCall")
    }, # initialization end


    processing = function(){
      super$processing()
      ref <- Biostrings::readDNAStringSet(private$paramlist[["Seq_file"]])
      pwm <- TFBSTools::toPWM(private$paramlist[["Motif"]])
      sitesetList <- TFBSTools::searchSeq(pwm, ref)
      tmp <- as(sitesetList, "data.frame")
      write.table(tmp, private$paramlist[["output_file"]], row.names = FALSE,  quote = FALSE)
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
      if(is.null(private$paramlist[["Seq_file"]])){
        stop("Parameter Seq_file is required!")
      }
      if(is.null(private$paramlist[["Motif"]])){
        stop("Parameter Motif is required!")
      }
    }
  ) # private end


) # R6 class end
