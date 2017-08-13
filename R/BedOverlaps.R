BedOverlaps <- R6::R6Class(
  classname = "BedOverlaps",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, BedInput1 = NULL, BedInput2 = NULL, Output = NULL,
                          n.col = NULL, editable=FALSE){
      super$initialize("BedOverlaps",editable,list(arg1=atacProc))

      if(!is.null(atacProc)){
        # Not using now
      }else{ # atacProc is NULL
        # necessary parameters
        private$paramlist[["BedInput1"]] <- BedInput1
        private$paramlist[["BedInput2"]] <- BedInput2
        private$paramlist[["col_num"]] <- n.col
        # unnecessary parameters
        if(is.null(Output)){
          tmp_path <- dirname(private$paramlist[["BedInput1"]])
          private$paramlist[["Output"]] <- paste0(tmp_path, "/Output.bed", collapse = "")
        }else{
          private$paramlist[["Output"]] <- Output
        }
      }

      # parameter check
      private$paramValidation()

    } # initialization end

  ), # public end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["BedInput1"]]))
      private$writeLog(sprintf("source:%s", private$paramlist[["BedInput2"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["Output"]]))

      BedIntersect(Input1 = private$paramlist[["BedInput1"]], Input2 = private$paramlist[["BedInput2"]],
                   bed_output = private$paramlist[["Output"]], n.col = private$paramlist[["col_num"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["BedInput1"]])){
        stop("Parameter BedInput1 is requied!");
      }
      if(is.null(private$paramlist[["BedInput2"]])){
        stop("Parameter BedInput2 is requied!");
      }
      if(is.null(private$paramlist[["col_num"]])){
        stop("Parameter n.col is requied!");
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["BedInput1"]])
      private$checkFileExist(private$paramlist[["BedInput2"]])
      private$checkFileCreatable(private$paramlist[["Output"]])
    }

  ) # private end

) # R6 class end


#' Using GenomicRanges function to find overlap of 2 bed files.
#'
#' Find the overlap peak of the Input1 file compared with Input2 file.
#' The two input file must have the same column, because the function will processing 3 and 6 column file with
#' difference methods according to the strand column.
#' The function using 0-based coordinate as bed file.
#' @param atacProc Do not use this parameter, we will add nore functions in the future!
#' @param BedInput1 The peak bed file(3-6 columns) that you want to find difference compared with the reference bed file.
#' @param BedInput2 The reference peak bed file.
#' @param Output The output bed file.
#' @param n.col How many columns of the input bed file, only support 3 or 6.
#' @export
atacBedOverlaps <- function(atacProc = NULL, BedInput1 = NULL, BedInput2 = NULL, Output = NULL,
                            n.col = NULL){
  tmp <- BedOverlaps$new(atacProc, BedInput1, BedInput2, Output, n.col)
  tmp$process()
  return(tmp)
}
