BedOverlaps <- R6::R6Class(
  classname = "BedOverlaps",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, bedInput1 = NULL, bedInput2 = NULL,
                          bedOutput = NULL, editable=FALSE){
      super$initialize("BedOverlaps", editable, list(arg1 = atacProc))

      if(!is.null(atacProc)){
        # Not using now
      }else{ # atacProc is NULL
        # necessary parameters
        private$paramlist[["bedInput1"]] <- bedInput1
        private$paramlist[["bedInput2"]] <- bedInput2
        # unnecessary parameters
        if(is.null(bedOutput)){
          tmp_path <- dirname(private$paramlist[["bedInput1"]])
          private$paramlist[["bedOutput"]] <- paste0(tmp_path, "/Output.bed", collapse = "")
        }else{
          private$paramlist[["bedOutput"]] <- bedOutput
        }
      }

      # parameter check
      private$paramValidation()

    } # initialization end

  ), # public end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["bedInput1"]]))
      private$writeLog(sprintf("source:%s", private$paramlist[["bedInput2"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["bedOutput"]]))

      bed1 <- rtracklayer::import(con = private$paramlist[["bedInput1"]], format = "bed")
      bed2 <- rtracklayer::import(con = private$paramlist[["bedInput2"]], format = "bed")
      output_data <- subsetByOverlaps(bed1, bed2, ignore.strand = TRUE)
      write.table(output_data, file = private$paramlist[["bedOutput"]],
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["bedInput1"]])){
        stop("Parameter bedInput1 is requied!");
      }
      if(is.null(private$paramlist[["bedInput2"]])){
        stop("Parameter bedInput2 is requied!");
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["bedInput1"]])
      private$checkFileExist(private$paramlist[["bedInput2"]])
      private$checkFileCreatable(private$paramlist[["bedOutput"]])
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
#' @param BedInput1 The bed file that you want to find difference compared with the reference bed file.
#' @param BedInput2 The reference bed file.
#' @param Output The output bed file.
#' @export
atacBedOverlaps <- function(atacProc = NULL, bedInput1 = NULL,
                            bedInput2 = NULL, bedOutput = NULL){
  tmp <- BedOverlaps$new(atacProc, bedInput1, bedInput2, bedOutput)
  tmp$process()
  return(tmp)
}
