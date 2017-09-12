Rsortbam <- R6::R6Class(
  classname = "Rsortbam",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, bamInput = NULL, bamOutput = NULL, editable=FALSE){
      super$initialize("Rsortbam", editable, list(arg1 = atacProc))

      # necessary parameters
      if(!is.null(atacProc)){ # atacproc from SamToBam
        private$paramlist[["bamInput"]] <- atacProc$getParam("bamOutput")
        regexProcName <- sprintf("(bam|%s)", atacProc$getProcName())
      }else if(is.null(atacProc)){ # input
        private$paramlist[["bamInput"]] <- inputbam
        regexProcName <- "(bam)"
      }
      # unnecessary parameters
      if(is.null(bamOutput)){
        prefix <- private$getBasenamePrefix(private$paramlist[["bamInput"]], regexProcName)
        private$paramlist[["bamOutput_tmp"]] <- file.path(.obtainConfigure("tmpdir"),
                                                          paste0(prefix, ".", self$getProcName()))
        private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["bamOutput_tmp"]], ".bam", collapse = "")
      }else{
        name_split <- unlist(base::strsplit(x = bamOutput, split = ".", fixed = TRUE))
        suffix <- tail(name_split, 1)
        if(suffix == "bam"){
          private$paramlist[["bamOutput_tmp"]] <- paste0(name_split[-length(name_split)], collapse = ".")
          private$paramlist[["bamOutput"]] <- bamOutput
        }else{
          private$paramlist[["bamOutput_tmp"]] <- bamOutput
          private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["bamOutput_tmp"]], ".bam", collapse = "")
        }
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
      Rsamtools::sortBam(file = private$paramlist[["bamInput"]], destination = private$paramlist[["bamOutput_tmp"]])
      Rsamtools::indexBam(files = private$paramlist[["bamOutput"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["bamInput"]])){
        stop("Parameter bamInput or atacProc is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["bamInput"]])
      private$checkPathExist(private$paramlist[["bamOutput"]])
    } # checkAllPath end

  ) # private end

) # class end

#' @name BamSort
#' @title Sort bam file.
#' @description This function sort a bam file and create a new .bai index.
#' @param atacProc Result from SamToBam.
#' @param samfile Path to bam file.
#' @param bamfile Where you want save the new bam file.
#' @export
atacBamSort <- function(atacProc = NULL, bamInput = NULL, bamOutput = NULL){
  tmp <- Rsortbam$new(atacProc, bamInput, bamOutput)
  tmp$process()
  invisible(tmp)
}
