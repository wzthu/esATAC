CutSiteCountR <- R6::R6Class(
  classname = "CutSiteCountR",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProcCutSite, csInput = NULL,
                          motifInput = NULL, matrixOutput = NULL, motifLength = NULL,
                          strandLength = NULL, FootPrint = NULL, editable = FALSE){
      super$initialize("CutSiteCountR", editable, list(arg1 = atacProcCutSite))

      # necessary parameters
      if(!is.null(atacProcCutSite)){
        private$paramlist[["csInput"]] <- atacProcCutSite$getParam("csOutput");
      }else{
        private$paramlist[["csInput"]] <- csInput
      }
      private$paramlist[["motifInput"]] <- motifInput
      private$paramlist[["matrixOutput"]] <- matrixOutput
      private$paramlist[["motifLength"]] <- motifLength
      private$paramlist[["strandLength"]] <- strandLength
      private$paramlist[["FootPrint"]] <- FootPrint
      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end

  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("Cut site:%s", private$paramlist[["csInput"]]))
      private$writeLog(sprintf("Motif:%s", private$paramlist[["motifInput"]]))
      private$writeLog(sprintf("Matrix destination:%s", private$paramlist[["matrixOutput"]]))
      tmp_dir <- paste(tempdir(), "/", Sys.getpid(), sep="")
      # using tmp dir to save temp data
      dir.create(tmp_dir, FALSE, TRUE, "0700")
      .chr_separate_call(ReadsIfile = private$paramlist[["motifInput"]],
                         ReadsOpath = tmp_dir,
                         Name = "/Motif")
      motif_tmp <- paste(tmp_dir, "/Motif", sep = "")

      chr <- as.list(c(1:22, "X", "Y"))
      for(i in seq(1:24)){
        echo_str <- paste("Now, processing chr", chr[[i]], "......", sep = "")
        print(echo_str)
        CutSiteInput <- paste0(private$paramlist[["csInput"]], "_chr", chr[[i]], ".cs", collapse = "")
        MotifInput <- paste0(motif_tmp, "_chr", chr[[i]], ".bed", collapse = "")
        MatrixOutput <- paste0(private$paramlist[["matrixOutput"]], "_chr", chr[[i]], ".matrix", collapse = "")
        .CutSiteCount(readsfile = CutSiteInput, motiffile = MotifInput, matrixfile = MatrixOutput,
                      motif_len = private$paramlist[["motifLength"]], strand_len = private$paramlist[["strandLength"]])
        if(i == 1){
          data <- try(read.table(MatrixOutput), silent = TRUE)
          if(inherits(data, "try-error")){
            data <- data.frame()
          }
        }else{
          temp <- try(read.table(MatrixOutput), silent = TRUE)
          if(inherits(temp, "try-error")){
            temp <- data.frame()
          }
          data <- rbind(data, temp)
        }
      }
      OverallMatrix <- paste0(private$paramlist[["matrixOutput"]], "data", ".matrix", collapse = "")
      write.table(data, file = OverallMatrix, row.names = FALSE, col.names = FALSE, quote = FALSE)
      if(private$paramlist[["FootPrint"]]){
        fp <- apply(data, 2, sum)
        plot(fp, type = "l", col = "blue", lwd = 2, xlab = "Relative Distance From Motif (bp)", ylab = "Cut Site Count", xaxt = "n", yaxt = "n")
        axis(1, at = seq(1, private$paramlist[["strandLength"]], len = 3),
             labels = -(private$paramlist[["strandLength"]] + 1 - seq(1, private$paramlist[["strandLength"]] + 1, len = 3)),
             padj = -1.0, tck = -0.01)
        axis(1, at = private$paramlist[["strandLength"]] + private$paramlist[["motifLength"]] + seq(1, private$paramlist[["strandLength"]], len = 3),
             labels = seq(0, private$paramlist[["strandLength"]], len = 3),
             padj = -1.0, tck = -0.01)
        axis(2, padj = 1.0,tck = -0.02)
        abline(v = c(private$paramlist[["strandLength"]], private$paramlist[["strandLength"]] + private$paramlist[["motifLength"]] + 1),
               lty = 2)
      }
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["csInput"]])){
        stop("Parameter csInput is required!")
      }
      if(is.null(private$paramlist[["motifInput"]])){
        stop("Parameter motifInput is required!")
      }
      if(is.null(private$paramlist[["matrixOutput"]])){
        stop("Parameter matrixOutput is required!")
      }
      if(is.null(private$paramlist[["motifLength"]])){
        stop("Parameter motifLength is required!")
      }
      if(is.null(private$paramlist[["strandLength"]])){
        stop("Parameter strandLength is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkPathExist(private$paramlist[["csInput"]])
      private$checkPathExist(private$paramlist[["motifInput"]])
      private$checkPathExist(private$paramlist[["matrixOutput"]])
    } # checkAllPath end

  ) # private end

) # class end


#' Counting cut site around motif.
#' @param atacProcCutSite Result from function atacCutSitePre.
#' @param atacProcMotif Result from function MotifScan.
#' @param csInput Your cut site information file(from atacCutSitePre function) path with prefix.
#' e.g. "/your_cut_site_information_path/prefix"
#' @param motifInput Your motif information file(from MotifScan function, set position = TRUE).
#' The first 4 columns of the motif file must be "chr start_site end_site strand".
#' @param matrixOutput The output path with a prefix, an empty folder would be great.
#' e.g. "/where_you_want_to_save_output/prefix"
#' @param motifLength Motif length.
#' @param strandLength How many bp(base pair) do you want to count up/downstream of the motif.
#' @param FootPrint TRUE or FALSE, plot footprint or not.
#' @export
atacCutSiteCount <- function(atacProcCutSite = NULL, csInput = NULL,
                             motifInput = NULL, matrixOutput = NULL, motifLength = NULL,
                             strandLength = NULL, FootPrint = FALSE){
  tmp <- CutSiteCountR$new(atacProcCutSite, csInput, motifInput,
                           matrixOutput, motifLength, strandLength, FootPrint)
  tmp$process()
  return(tmp)
}
