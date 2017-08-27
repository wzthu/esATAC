BedOverlaps <- R6::R6Class(
  classname = "BedOverlaps",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, bedInput1 = NULL, bedInput2 = NULL, bedOutput = NULL,
                          n.col = NULL, editable=FALSE){
      super$initialize("BedOverlaps",editable,list(arg1=atacProc))

      if(!is.null(atacProc)){
        # Not using now
      }else{ # atacProc is NULL
        # necessary parameters
        private$paramlist[["bedInput1"]] <- bedInput1
        private$paramlist[["bedInput2"]] <- bedInput2
        private$paramlist[["col_num"]] <- n.col
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

      BedIntersect(Input1 = private$paramlist[["bedInput1"]], Input2 = private$paramlist[["bedInput2"]],
                   bed_output = private$paramlist[["bedOutput"]], n.col = private$paramlist[["col_num"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["bedInput1"]])){
        stop("Parameter bedInput1 is requied!");
      }
      if(is.null(private$paramlist[["bedInput2"]])){
        stop("Parameter bedInput2 is requied!");
      }
      if(is.null(private$paramlist[["col_num"]])){
        stop("Parameter n.col is requied!");
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
#' @param Input1 The peak bed file(3-6 columns) that you want to find difference compared with the reference bed file.
#' @param Input2 The reference peak bed file.
#' @param bed_output The output bed file.
#' @param n.col How many columns of the input bed file.
#' @export
BedIntersect <- function(Input1, Input2, bed_output, n.col){
  bed_1 <- read.table(Input1, header = FALSE)
  bed_2 <- read.table(Input2, header = FALSE)
  if(n.col == 3){
    colnames(bed_1) <- c('chr', 'start', 'end')
    colnames(bed_2) <- c('chr', 'start', 'end')
    data1 <- with(bed_1, GRanges(chr, IRanges(start, end)))
    data2 <- with(bed_2, GRanges(chr, IRanges(start, end)))
    output_data <- subsetByOverlaps(data1, data2, ignore.strand = TRUE)
    write.table(output_data, file = bed_output,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }else if(n.col == 6){
    colnames(bed_1) <- c('chr', 'start', 'end', 'id', 'score', 'strand')
    colnames(bed_2) <- c('chr', 'start', 'end', 'id', 'score', 'strand')
    data1 <- with(bed_1, GRanges(chr, IRanges(start, end), strand, score, id = id))
    data2 <- with(bed_2, GRanges(chr, IRanges(start, end), strand, score, id = id))
    output_data <- subsetByOverlaps(data1, data2, ignore.strand = FALSE)
    write.table(output_data, file = bed_output,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}



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
atacBedOverlaps <- function(atacProc = NULL, bedInput1 = NULL, bedInput2 = NULL, bedOutput = NULL,
                            n.col = NULL){
  tmp <- BedOverlaps$new(atacProc, bedInput1, bedInput2, bedOutput, n.col)
  tmp$process()
  return(tmp)
}
