RSNPs <- R6::R6Class(
  classname = "RSNPs",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, snp.info = NULL, peak.info = NULL,
                          annoOutput = NULL, editable = FALSE){
      super$initialize("RSNPs",editable,list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){
        private$paramlist[["peak.info"]] <- atacProc$getParam("bedOutput")
      }else{
        private$paramlist[["peak.info"]] <- peak.info
      }
      private$paramlist[["snp.info"]] <- snp.info
      # unnecessary parameters
      if(is.null(annoOutput)){
        private$paramlist[["annoOutput"]] <- paste0(dirname(private$paramlist[["peak.info"]]),
                                                    "/SNPAnno", collapse = "")
      }else{
        private$paramlist[["annoOutput"]] <- annoOutput
      }
      # parameter check
      private$paramValidation()
    } # initialization end

  ), # piblic end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("SNP source:%s", private$paramlist[["snp.info"]]))
      private$writeLog(sprintf("Peak source:%s", private$paramlist[["peak.info"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["annoOutput"]]))
      SNP_info <- read.table(file = private$paramlist[["snp.info"]],
                             header = FALSE)
      peak_info <- read.table(file = private$paramlist[["peak.info"]],
                              header = FALSE)
      snp_gr <- with(SNP_info, GRanges(V1, IRanges(V2 - 1, V2)))
      peak_gr <- with(peak_info, GRanges(V1, IRanges(V2, V3)))
      overlaps <- GenomicRanges::findOverlaps(query = peak_gr, subject = snp_gr)
      output <- cbind(peak_info[S4Vectors::queryHits(overlaps), ],
                      SNP_info[S4Vectors::subjectHits(overlaps), ])
      write.table(x = output, file = private$paramlist[["annoOutput"]],
                  quote = FALSE, sep = "\t", row.names = FALSE,
                  col.names = FALSE)
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["peak.info"]])){
        stop("Parameter peak.info is required!")
      }
      if(is.null(private$paramlist[["snp.info"]])){
        stop("Parameter snp.info is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["snp.info"]])
      private$checkPathExist(private$paramlist[["annoOutput"]])
    } # checkAllPath end

  ) # private end

) # class end


#' Find snps in the given peak file
#'
#' This function do not consider strand.
#' @param atacProc Result from function "PeakCallingFseq".
#' @param snp.info Path to your snps file.the first 2 column must be snp position
#' (chromatin and  position).
#' @param peak.info Path to your peak file.
#' @annoOutput Where to save annotation information, default: peak.info_dir/SNPAnno.
SNPAnno <- function(atacProc = NULL, snp.info = NULL, peak.info = NULL,
                   annoOutput = NULL){
  tmp <- RSNPs$new(atacProc, snp.info, peak.info, annoOutput)
  tmp$process()
  return(tmp)
}
