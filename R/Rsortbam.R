Rsortbam <- R6::R6Class(
  classname = "Rsortbam",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, inputbam = NULL, outputbam = NULL, editable=FALSE){
      super$initialize("Rsortbam",editable,list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){ # atacproc from SamToBam
        private$paramlist[["bamInput"]] <- atacProc$getParam("bamOutput")
      }else if(is.null(atacProc)){ # input
        private$paramlist[["bamInput"]] <- inputbam
      }
      # unnecessary parameters
      if(is.null(outputbam)){
        private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["bamInput"]], ".sorted", collapse = "")
      }else{
        private$paramlist[["bamOutput"]] <- outputbam
      }

      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end

  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["bamInput"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["bamOutput"]]))
      Rsamtools::sortBam(file = private$paramlist[["bamInput"]], destination = private$paramlist[["bamOutput"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["bamInput"]])){
        stop("Parameter inputbam is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["bamInput"]])
      private$checkPathExist(private$paramlist[["bamOutput"]])
    } # checkAllPath end

  ) # private end

) # class end

#' sortbam using Rbowtie::sortBam
#' the output bam file do not have bam header, can not use QCreport function
#' @export
atacBamSort <- function(atacProc = NULL, inputbam = NULL, outputbam = NULL){
  tmp <- Rsortbam$new(atacProc, inputbam, outputbam)
  tmp$process()
  return(tmp)
}
