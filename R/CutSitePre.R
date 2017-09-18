CutSitePre <- R6::R6Class(
  classname = "CutSitePre",
  inherit = ATACProc,

  public = list(
    initialize = function(atacProc, bedInput = NULL, csOutput.dir = NULL, prefix = NULL, editable = FALSE){
      super$initialize("CutSitePre",editable,list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){ # get parameter from class BamToBed or SamToBed
        private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput")
        regexProcName <- sprintf("(bed|%s)", atacProc$getProcName())
      }else{
        private$paramlist[["bedInput"]] <- bedInput
        regexProcName <- "(bed)"
      }
      # unnecessary parameters
      if(is.null(csOutput.dir)){
        dir_name <- private$getBasenamePrefix(private$paramlist[["bedInput"]], regexProcName)
        private$paramlist[["csOutput.dir"]] <- file.path(.obtainConfigure("tmpdir"),
                                                         dir_name)
        dir.create(private$paramlist[["csOutput.dir"]])
      }else{
        private$paramlist[["csOutput.dir"]] <- csOutput.dir
      }
      if(is.null(prefix)){
        private$paramlist[["csOutput"]] <- paste0(private$paramlist[["csOutput.dir"]], "/", "Cutsite", collapse = "")
      }else{
        private$paramlist[["csOutput"]] <- paste0(private$paramlist[["csOutput.dir"]], "/", prefix, collapse = "")
      }
      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["bedInput"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["csOutput"]]))
      Cut_Site_Number <- .CutSite_call(InputFile = private$paramlist[["bedInput"]], OutputFile = private$paramlist[["csOutput"]])
      print(sprintf("Total Cut Site Number is: %d", Cut_Site_Number))
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["bedInput"]])){
        stop("Parameter bedInput is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["bedInput"]])
      private$checkPathExist(private$paramlist[["csOutput"]])
    } # checkAllPath end

  ) # private end

) # class end


#' extract cut site information
#' @param atacProc Do not use this parameter, we will add nore functions in the future!
#' @param bedInput Input file path, the No.1-3 column is chromatin name, start cut site, end cut site.
#' The merged bed file is recommented.
#' @param csOutput.dir The output path, an empty folder would be great.
#' @param prefix Output file name prefix, e.g. prefix_chr*.bed, default "Cutsite".
#' @export
atacCutSitePre <- function(atacProc = NULL, bedInput = NULL, csOutput.dir = NULL, prefix = NULL){
  tmp <- CutSitePre$new(atacProc, bedInput, csOutput.dir, prefix)
  tmp$process()
  invisible(tmp)
}

