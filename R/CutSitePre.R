CutSitePre <- R6::R6Class(
  classname = "CutSitePre",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, bedInput = NULL, csOutput = NULL, prefix = NULL, editable = FALSE){
      super$initialize("CutSitePre",editable,list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){ # get parameter from class BamToBed or SamToBed
        private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
      }else{
        private$paramlist[["bedInput"]] <- bedInput
      }
      # unnecessary parameters
      if(is.null(prefix)){
        private$paramlist[["csOutput"]] <- paste0(csOutput, "/", "output", collapse = "")
      }else{
        private$paramlist[["csOutput"]] <- paste0(csOutput, "/", prefix, collapse = "")
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
      .CutSite_call(InputFile = private$paramlist[["bedInput"]], OutputFile = private$paramlist[["csOutput"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["bedInput"]])){
        stop("Parameter bedInput is required!")
      }
      if(is.null(private$paramlist[["csOutput"]])){
        stop("Parameter csOutput is required! An empty file is recommented!")
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
#' @param csOutput The output path, an empty folder would be great.
#' @param prefix Output file name prefix, e.g. prefix_chr*.bed, default "output".
#' @export
atacCutSitePre <- function(atacProc = NULL, bedInput = NULL, csOutput = NULL, prefix = NULL){
  tmp <- CutSitePre$new(atacProc, bedInput, csOutput, prefix)
  tmp$process()
  return(tmp)
}

