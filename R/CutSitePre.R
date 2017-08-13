CutSitePre <- R6::R6Class(
  classname = "CutSitePre",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, InputFile = NULL, OutputFile = NULL, prefix = NULL, editable = FALSE){
      super$initialize("CutSitePre",editable,list(arg1=atacProc))

      # necessary parameters
      private$paramlist[["Input"]] <- InputFile
      # unnecessary parameters
      if(is.null(prefix)){
        private$paramlist[["Output"]] <- paste0(OutputFile, "/", "output_", collapse = "")
      }else{
        private$paramlist[["Output"]] <- paste0(OutputFile, "/", prefix, "_", collapse = "")
      }

      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["Input"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["Output"]]))
      .CutSite_call(InputFile = private$paramlist[["Input"]], OutputFile = private$paramlist[["Output"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["Input"]])){
        stop("Parameter InputFile is required!")
      }
      if(is.null(private$paramlist[["Output"]])){
        stop("Parameter OutputFile is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["Input"]])
      private$checkPathExist(private$paramlist[["Output"]])
    } # checkAllPath end

  ) # private end

) # class end


#' extract cut site information
#' @param atacProc Do not use this parameter, we will add nore functions in the future!
#' @param InputFile Input file path, the No.1-3 column is chromatin name, start cut site, end cut site.
#' @param OutputFile The output path, an empty folder would be great.
#' @param prefix Output file name prefix, e.g. prefix_chr*.bed, default "output".
#' @export
atacCutSitePre <- function(atacProc = NULL, InputFile = NULL, OutputFile = NULL, prefix = NULL){
  tmp <- CutSitePre$new(atacProc, InputFile, OutputFile, prefix)
  tmp$process()
  return(tmp)
}

