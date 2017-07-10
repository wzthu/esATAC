CountCutSiteR <- R6::R6Class(
  classname = "CountCutSiteR",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, ForwReads_dir = NULL, RevReads_dir = NULL, Motif_dir = NULL, ForwMatrix_dir = NULL,
                          RevMatrix_dir = NULL, Motif_length = NULL, Strand_length = NULL, editable = FALSE){
      super$initialize("CountCutSiteR",editable,list(arg1=atacProc))
      print("CountCutSiteRInitCall")
      # necessary and unchanged parameters, this should be tested
      private$paramlist[["FR_dir"]] <- ForwReads_dir
      private$paramlist[["RR_dir"]] <- RevReads_dir
      private$paramlist[["M_dir"]] <- Motif_dir
      private$paramlist[["FM_dir"]] <- ForwMatrix_dir
      private$paramlist[["RM_dir"]] <- RevMatrix_dir
      private$paramlist[["M_len"]] <- Motif_length
      private$paramlist[["S_len"]] <- Strand_length
      private$checkRequireParam()
      if(!dir.exists(private$paramlist[["FR_dir"]])){
        stop(paste("error, path does not exist:",private$paramlist[["FR_dir"]]))
      }
      if(!dir.exists(private$paramlist[["RR_dir"]])){
        stop(paste("error, path does not exist:",private$paramlist[["RR_dir"]]))
      }
      if(!dir.exists(private$paramlist[["M_dir"]])){
        stop(paste("error, path does not exist:",private$paramlist[["M_dir"]]))
      }
      if(!dir.exists(private$paramlist[["FM_dir"]])){
        stop(paste("error, path does not exist:",private$paramlist[["FM_dir"]]))
      }
      if(!dir.exists(private$paramlist[["RM_dir"]])){
        stop(paste("error, path does not exist:",private$paramlist[["RM_dir"]]))
      }
      print("finishMappingInitCall")

    }, # end of initialize

    processing = function(){



    }, # end of processing


    setResultParam = function(bedFilePath){



    }  # end of setResultParam

  ), # end of public

  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return()
      }

      if(is.null(private$paramlist[["FR_dir"]])){
        stop("Forward strand reads file path is required!")
      }

      if(is.null(private$paramlist[["RR_dir"]])){
        stop("Reverse strand reads file path is required!")
      }

      if(is.null(private$paramlist[["M_dir"]])){
        stop("Motif position file path is required!")
      }

      if(is.null(private$paramlist[["FM_dir"]])){
        stop("Forward strand cut site matrix file path is required!")
      }

      if(is.null(private$paramlist[["RM_dir"]])){
        stop("Reverse strand cut site matrix file path is required!")
      }

      if(is.null(private$paramlist[["M_len"]])){
        stop("Motif length is required!")
      }

      if(is.null(private$paramlist[["S_len"]])){
        stop("Strand length is required(how many bp do you want to count
           upstream/downstream of the motif center)!")
      }

    }

  ) # end of private



) # end of R6 class
