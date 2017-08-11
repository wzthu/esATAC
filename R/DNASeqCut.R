DNASeqCut <- R6::R6Class(
  classname = "DNASeqCut",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, ref_path = NULL, save_path = NULL,
                          bed_path = NULL, save_format = NULL,
                          editable = FALSE){
      super$initialize("DNASeqCut",editable,list(arg1=atacProc))
      print("DNASeqCutInitCall")
      # necessary and unchanged parameters, this should be tested
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }
      private$paramlist[["ref_path"]] <- ref_path
      private$paramlist[["save_path"]] <- save_path
      private$paramlist[["bed_path"]] <- bed_path
      private$paramlist[["save_format"]] <- save_format
      # check parameter
      private$checkRequireParam()
      private$checkFileExist(private$paramlist[["ref_path"]])
      private$checkFileExist(private$paramlist[["bed_path"]])
      private$checkPathExist(private$paramlist[["save_path"]])
      print("finishDNASeqCutInitCall")

    }, # # initialization end

    processing = function(){
      if(!super$processing()){
        return()
      }
      Sequence_Cut(private$paramlist[["ref_path"]], private$paramlist[["save_path"]],
                   private$paramlist[["bed_path"]], private$paramlist[["save_format"]])
      private$setFinish()
    }, # processing end

    setResultParam = function(){
      print("This function is not using!")
    }  # setResultParam end

  ), # public end


  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return()
      }
      if(is.null(private$paramlist[["save_format"]])){
        stop("Parameter save_format is required!")
      }
    }
  ) # private end

) # class end
