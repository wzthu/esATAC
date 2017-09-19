RPeakComp <- R6::R6Class(
  classname = "RPeakComp",
  inherit = ATACProc,
  public = list(
    initialize = function(atacProc, bedInput1 = NULL, bedInput2 = NULL,
                          bedOutput = NULL, operation = NULL,
                          editable = FALSE){
      super$initialize("RPeakComp", editable, list(arg1 = atacProc))

      if(!is.null(atacProc)){
        # Not using now
      }else{ # atacProc is NULL
        # necessary parameters
        private$paramlist[["bedInput1"]] <- bedInput1
        private$paramlist[["bedInput2"]] <- bedInput2
        private$paramlist[["operation"]] <- operation
        # unnecessary parameters
        if(private$paramlist[["operation"]] == "overlap"){
          if(is.null(bedOutput)){
            tmp_path <- dirname(private$paramlist[["bedInput1"]])
            private$paramlist[["bedOutput"]] <- paste0(tmp_path, "/Output_overlap",
                                                       collapse = "")
          }else{
            private$paramlist[["bedOutput"]] <- bedOutput
          }
        }else if(private$paramlist[["operation"]] == "diff"){
          private$paramlist[["bedOutput"]] <- list()
          if(is.null(bedOutput)){
            private$paramlist[["bedOutput"]][[1]] <- paste(private$paramlist[["bedInput1"]],
                                                           "_diff", sep = "")
            private$paramlist[["bedOutput"]][[2]] <- paste(private$paramlist[["bedInput2"]],
                                                           "_diff", sep = "")
          }else{
            private$paramlist[["bedOutput"]][[1]] <- bedOutput[[1]]
            private$paramlist[["bedOutput"]][[2]] <- bedOutput[[2]]
          }
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
      if(private$paramlist[["operation"]] == "overlap"){
        output_data <- GenomicRanges::intersect(bed1, bed2, ignore.strand = TRUE)
        rtracklayer::export(object = output_data,
                            con = private$paramlist[["bedOutput"]],
                            format = "bed")
      }else{
        overlaps <- GenomicRanges::findOverlaps(query = bed1, subject = bed2,
                                                ignore.strand = TRUE)
        bed1_index <- seq(length(bed1))
        bed2_index <- seq(length(bed2))
        bed1_output <- bed1[setdiff(bed1_index, S4Vectors::queryHits(overlaps))]
        bed2_output <- bed2[setdiff(bed2_index, S4Vectors::subjectHits(overlaps))]
        rtracklayer::export(object = bed1_output,
                            con = private$paramlist[["bedOutput"]][[1]],
                            format = "bed")
        rtracklayer::export(object = bed2_output,
                            con = private$paramlist[["bedOutput"]][[2]],
                            format = "bed")
      }

    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["bedInput1"]])){
        stop("Parameter bedInput1 is requied!")
      }
      if(is.null(private$paramlist[["bedInput2"]])){
        stop("Parameter bedInput2 is requied!")
      }
      if(is.null(private$paramlist[["operation"]])){
        stop("Parameter operation is requied!")
        stopifnot(private$paramlist[["operation"]] %in% c("overlap", "diff"))
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["bedInput1"]])
      private$checkFileExist(private$paramlist[["bedInput2"]])
    }

  ) # private end

) # R6 class end

#' Peak comparation
#'
#' Find the overlap or differential pewks of the imput files.
#' The function use 0-based coordinate in bed file.
#' @param atacProc Do not use this parameter, we will add nore functions in the future!
#' @param bedInput1 The bed file that you want to find difference compared with the reference bed file.
#' @param bedInput2 The reference bed file.
#' @param bedOutput The output bed file.
#' If parameter operation is "diff", bedOutput should be a list containing 2 file path.
#' Default bedInput1_diff and bedInput2_diff.
#' If parameter operation is "overlap", bedOutput should be a vector.
#' Default Output_overlap.
#' @param operation "overlap" or "diff". If choose "overlap", the output will be a bed file contain the exact overlap
#' of the two peak set. If choose "diff", the output will be 2 bed files contain peaks without any overlap of
#' the two input bed files.
#' @export
PeakComp <- function(atacProc = NULL, bedInput1 = NULL, bedInput2 = NULL,
                     bedOutput = NULL, operation = NULL){
  tmp <- RPeakComp$new(atacProc, bedInput1, bedInput2, bedOutput, operation)
  tmp$process()
  invisible(tmp)
}
