BamToBed <- R6::R6Class(
  classname = "BamToBed",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, bamfile = NULL, bedfile = NULL, editable=FALSE){
      super$initialize("BamToBed",editable,list(arg1=atacProc))
      print("BamToBedInitCall")
      private$checkRequireParam()

      if((!is.null(atacProc)) && (class(atacProc)[1] == "SamToBam")){ # atacproc from mapping
        private$paramlist[["bamInput"]] <- atacProc$getParam("bamOutput")
      }else if((!is.null(atacProc)) && (class(atacProc)[1] != "SamToBam")){ # atacproc not from mapping, error!
        stop("Input class must be got from 'SamToBam'!")
      }else if(is.null(atacProc)){ # input
        private$paramlist[["bamInput"]] <- bamfile
      }

      if(is.null(bedfile)){
        private$paramlist[["bedOutput"]] <- paste0(private$paramlist[["bamInput"]], ".bed", collapse = "")
      }else{
        private$paramlist[["bedOutput"]] <- bedfile
      }

      private$checkFileExist(private$paramlist[["bamInput"]])
      private$checkPathExist(private$paramlist[["bedOutput"]])

      print("finishBamToBedInitCall")
    }, # initialize end


    processing = function(){
      super$processing()
      rtracklayer::export(con = private$paramlist[["bedOutput"]],
                          object = rtracklayer::import(con = private$paramlist[["bamInput"]], format = "bam", paired = TRUE, use.names = TRUE),
                          format = "bed")
      private$finish <- TRUE
    }, # processing end




    setResultParam = function(bedFilePath){
      super$setResultParam()
      private$paramlist[["bamOutput"]] <- bedFilePath
    }

  ), # public end



  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return()
      }
    }

  ) # private end










) # R6 class end
