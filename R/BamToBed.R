BamToBed <- R6::R6Class(
  classname = "BamToBed",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, bamInput = NULL, bedOutput = NULL, editable=FALSE){
      super$initialize("BamToBed", editable, list(arg1 = atacProc))

      # necessary parameters
      if(!is.null(atacProc)){ # atacproc from mapping
        private$paramlist[["bamInput"]] <- atacProc$getParam("bamOutput")
        regexProcName <- sprintf("(bam|%s)", atacProc$getProcName())
      }else if(is.null(atacProc)){ # input
        private$paramlist[["bamInput"]] <- bamInput
        regexProcName <- "(bam)"
      }
      # unnecessary parameters
      if(is.null(bedOutput)){
        prefix <- private$getBasenamePrefix(private$paramlist[["bamInput"]], regexProcName)
        private$paramlist[["bedOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                      paste0(prefix, ".", self$getProcName(), ".bed"))
      }else{
        name_split <- unlist(base::strsplit(x = bedOutput, split = ".", fixed = TRUE))
        suffix <- tail(name_split, 1)
        if(suffix == "bed"){
          private$paramlist[["bedOutput"]] <- bedOutput
        }else{
          private$paramlist[["bedOutput"]] <- paste0(bedOutput, ".bed", collapse = "")
        }
      }

      # parameter check
      private$paramValidation()

    } # initialize end

  ), # public end



  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["bamInput"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["bedOutput"]]))
      rtracklayer::export(con = private$paramlist[["bedOutput"]],
                          object = rtracklayer::import(con = private$paramlist[["bamInput"]], format = "bam", paired = TRUE, use.names = TRUE),
                          format = "bed")
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["bamInput"]])){
        stop("Parameter bamInput is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["bamInput"]])
      private$checkFileCreatable(private$paramlist[["bedOutput"]])
    }

  ) # private end


) # R6 class end


#' @name atacBam2Bed
#' @title Convert bam format to bed format.
#' @description This function convert a bam file to a bed file.
#' @param atacProc Result from atacSam2Bam or atacBamSort.
#' @param bamInput Path to bam file.
#' @param bedOutput Where you want save the bed file.
#' @export
atacBam2Bed <- function(atacProc = NULL, bamInput = NULL, bedOutput = NULL){
  tmp <- BamToBed$new(atacProc, bamInput, bedOutput)
  tmp$process()
  return(tmp)
}
