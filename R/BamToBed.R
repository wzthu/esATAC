BamToBed <- R6::R6Class(
  classname = "BamToBed",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, bamInput = NULL, bedOutput = NULL, editable=FALSE){
      super$initialize("BamToBed", editable, list(arg1 = atacProc))

      # necessary parameters
      # CHECK atacProc!!!!!!!!!!!!!!!!!!!!!!
      if(!is.null(atacProc)){ # atacproc from mapping
        private$paramlist[["bamInput"]] <- atacProc$getParam("bamOutput")
      }else if(is.null(atacProc)){ # input
        private$paramlist[["bamInput"]] <- bamInput
      }
      # unnecessary parameters
      if(is.null(bedOutput)){
        private$paramlist[["bedOutput"]] <- paste0(private$paramlist[["bamInput"]], ".bed", collapse = "")
      }else{
        private$paramlist[["bedOutput"]] <- bedOutput
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


#' bam2bed using rtracklayer package
#' @export
atacBam2Bed <- function(atacProc = NULL, bamInput = NULL, bedOutput = NULL){
  tmp <- BamToBed$new(atacProc, bamInput, bedOutput)
  tmp$process()
  return(tmp)
}
