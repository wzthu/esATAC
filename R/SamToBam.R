SamToBam <-R6::R6Class(
  classname = "SamToBam",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, samfile = NULL, bamfile = NULL, editable=FALSE){
      super$initialize("SamToBam",editable,list(arg1=atacProc))

      # necessary parameters
      if((!is.null(atacProc)) ){
        private$paramlist[["samInput"]] <- atacProc$getParam("samOutput")
      }else if(is.null(atacProc)){ # input
        private$paramlist[["samInput"]] <- samfile
      }
      # unnecessary parameters
      if(is.null(bamfile)){
        private$paramlist[["bamOutput_tmp"]] <- private$paramlist[["samInput"]]
        private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["samInput"]], ".bam", collapse = "")
        private$paramlist[["baiOutput"]] <- paste0(private$paramlist[["samInput"]], ".bam.bai", collapse = "")
      }else{
        private$paramlist[["bamOutput_tmp"]] <- bamfile
        private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["bamOutput_tmp"]], ".bam", collapse = "")
        private$paramlist[["baiOutput"]] <- paste0(private$paramlist[["bamOutput_tmp"]], ".bam.bai", collapse = "")
      }


      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end

  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["samInput"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["bamOutput"]]))
      Rsamtools::asBam(file = private$paramlist[["samInput"]], destination = private$paramlist[["bamOutput_tmp"]])

    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["samInput"]])){
        stop("Parameter samInput is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["samInput"]])
      private$checkPathExist(private$paramlist[["bamOutput"]])
    } # checkAllPath end
  ) # private end

) # class end

#' sam2bam using Rbowtie::asBam
#' @export
atacSam2Bam <- function(atacProc = NULL, samfile = NULL, bamfile = NULL){
  tmp <- SamToBam$new(atacProc, samfile, bamfile)
  tmp$process()
  return(tmp)
}
