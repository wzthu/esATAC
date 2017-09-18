SamToBam <-R6::R6Class(
  classname = "SamToBam",
  inherit = ATACProc,
  public = list(
    initialize = function(atacProc, samfile = NULL, bamfile = NULL, editable=FALSE){
      super$initialize("SamToBam", editable, list(arg1 = atacProc))

      # necessary parameters
      if((!is.null(atacProc)) ){
        private$paramlist[["samInput"]] <- atacProc$getParam("samOutput")
        regexProcName <- sprintf("(sam|%s)", atacProc$getProcName())
      }else if(is.null(atacProc)){ # input
        private$paramlist[["samInput"]] <- samfile
        regexProcName <- "(sam)"
      }
      # unnecessary parameters
      if(is.null(bamfile)){
        prefix <- private$getBasenamePrefix(private$paramlist[["samInput"]], regexProcName)
        private$paramlist[["bamOutput_tmp"]] <- file.path(.obtainConfigure("tmpdir"),
                                                          paste0(prefix, ".", self$getProcName()))
        private$paramlist[["bamOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                      paste0(prefix, ".", self$getProcName(), ".bam"))
        private$paramlist[["baiOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                      paste0(prefix, ".", self$getProcName(), ".bam.bai"))
      }else{
        name_split <- unlist(base::strsplit(x = bamfile, split = ".", fixed = TRUE))
        suffix <- tail(name_split, 1)
        if(suffix == "bam"){
          private$paramlist[["bamOutput_tmp"]] <- paste0(name_split[-length(name_split)], collapse = ".")
          private$paramlist[["bamOutput"]] <- bamfile
          private$paramlist[["baiOutput"]] <- paste0(private$paramlist[["bamOutput"]], ".bai", collapse = "")
        }else{
          private$paramlist[["bamOutput_tmp"]] <- bamfile
          private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["bamOutput_tmp"]], ".bam", collapse = "")
          private$paramlist[["baiOutput"]] <- paste0(private$paramlist[["bamOutput"]], ".bai", collapse = "")
        }

      }

      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end

  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["samInput"]]))
      private$writeLog(sprintf("Bam destination:%s", private$paramlist[["bamOutput"]]))
      private$writeLog(sprintf("Bai destination:%s", private$paramlist[["baiOutput"]]))
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

#' @name atacSam2Bam
#' @title Convert sam format to bam format
#' @description This function convert sam file(from mapping) to bam file
#' (with a .bai index).
#' @param atacProc Result from Mapping.
#' @param samfile Path to sam file.
#' @param bamfile Where you want save the bam file.
#' @export
atacSam2Bam <- function(atacProc = NULL, samfile = NULL, bamfile = NULL){
  tmp <- SamToBam$new(atacProc, samfile, bamfile)
  tmp$process()
  invisible(tmp)
}
