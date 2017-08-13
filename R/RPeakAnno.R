RPeakAnno <- R6::R6Class(
  classname = "RPeakAnno",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, Input = NULL, Output = NULL, editable = FALSE){
      super$initialize("RPeakAnno",editable,list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }
      private$paramlist[["Input"]] <- Input
      # unnecessary parameters
      if(is.null(Output)){
        private$paramlist[["Output"]] <- paste0(dirname(Input), "/output", collapse = "")
      }else{
        private$paramlist[["Output"]] <- Output
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
      peak <- ChIPseeker::readPeakFile(private$paramlist[["Input"]])
      peakAn <- ChIPseeker::annotatePeak(peak)
      write.table(x = peakAn, file = private$paramlist[["Output"]], append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["Input"]])){
        stop("Parameter Input is required!")
      }
    },

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["Input"]])
      private$checkPathExist(private$paramlist[["Output"]])
    } # checkAllPath end

  ) # private end

) # R6 class end

#' Using ChIPseeker to annotate peak file.
#'
#' Just for test, other parameters will be add.
#' @param Input peak file.
#' @param Output Output file, default Input_path/output.
PeakAnno <- function(atacProc = NULL, Input = NULL,
                     Output = NULL){
  tmp <- RPeakAnno$new(atacProc, Input, Output)
  tmp$process()
  return(tmp)
}
