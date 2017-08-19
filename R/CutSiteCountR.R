CutSiteCountR <- R6::R6Class(
  classname = "CutSiteCountR",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, CutSiteFile = NULL, MotifFile = NULL, MatrixPath = NULL,
                          motif_length = NULL, strand_length = NULL, FootPrint = FALSE, editable = FALSE){
      super$initialize("CutSiteCountR",editable,list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }
      private$paramlist[["CutSiteFile"]] <- CutSiteFile
      private$paramlist[["MotifFile"]] <- MotifFile
      private$paramlist[["MatrixPath"]] <- MatrixPath
      private$paramlist[["motif_length"]] <- motif_length
      private$paramlist[["strand_length"]] <- strand_length
      private$paramlist[["FootPrint"]] <- FootPrint
      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end

  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("Cut site:%s", private$paramlist[["CutSiteFile"]]))
      private$writeLog(sprintf("Motif:%s", private$paramlist[["MotifFile"]]))
      private$writeLog(sprintf("Matrix destination:%s", private$paramlist[["MatrixPath"]]))
      tmp_dir <- paste(tempdir(), "/", Sys.getpid(), sep="")
      # using tmp dir to save temp data
      dir.create(tmp_dir, FALSE, TRUE, "0700")
      .chr_separate_call(ReadsIfile = private$paramlist[["MotifFile"]],
                         ReadsOpath = tmp_dir,
                         Name = "/Motif")
      motif_tmp <- paste(tmp_dir, "/Motif", sep = "")

      chr <- as.list(c(1:22, "X", "Y"))
      for(i in seq(1:24)){
        echo_str <- paste("Now, processing chr", chr[[i]], "......", sep = "")
        print(echo_str)
        CutSiteInput <- paste0(private$paramlist[["CutSiteFile"]], "_chr", chr[[i]], ".bed", collapse = "")
        MotifInput <- paste0(motif_tmp, "_chr", chr[[i]], ".bed", collapse = "")
        MatrixOutput <- paste0(private$paramlist[["MatrixPath"]], "_chr", chr[[i]], ".matrix", collapse = "")
        .CutSiteCount(readsfile = CutSiteInput, motiffile = MotifInput, matrixfile = MatrixOutput,
                      motif_len = private$paramlist[["motif_length"]], strand_len = private$paramlist[["strand_length"]])
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
      OverallMatrix <- paste0(private$paramlist[["MatrixPath"]], "data", ".matrix", collapse = "")
      write.table(data, file = OverallMatrix, row.names = FALSE, col.names = FALSE, quote = FALSE)
      if(private$paramlist[["FootPrint"]]){
        fp <- apply(data, 2, sum)
        plot(fp, type = "l", col = "blue", lwd = 2, xlab = "Relative Distance From Motif (bp)", ylab = "Cut Site Count", xaxt = "n", yaxt = "n")
        axis(1, at = seq(1, private$paramlist[["strand_length"]], len = 3),
             labels = -(private$paramlist[["strand_length"]] + 1 - seq(1, private$paramlist[["strand_length"]] + 1, len = 3)),
             padj = -1.0, tck = -0.01)
        axis(1, at = private$paramlist[["strand_length"]] + private$paramlist[["motif_length"]] + seq(1, private$paramlist[["strand_length"]], len = 3),
             labels = seq(0, private$paramlist[["strand_length"]], len = 3),
             padj = -1.0, tck = -0.01)
        axis(2, padj = 1.0,tck = -0.02)
        abline(v = c(private$paramlist[["strand_length"]], private$paramlist[["strand_length"]] + private$paramlist[["motif_length"]] + 1),
               lty = 2)
      }
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["CutSiteFile"]])){
        stop("Parameter CutSiteFile is required!")
      }
      if(is.null(private$paramlist[["MotifFile"]])){
        stop("Parameter MotifFile is required!")
      }
      if(is.null(private$paramlist[["MatrixPath"]])){
        stop("Parameter MatrixPath is required!")
      }
      if(is.null(private$paramlist[["motif_length"]])){
        stop("Parameter motif_length is required!")
      }
      if(is.null(private$paramlist[["strand_length"]])){
        stop("Parameter strand_length is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkPathExist(private$paramlist[["CutSiteFile"]])
      private$checkPathExist(private$paramlist[["MotifFile"]])
      private$checkPathExist(private$paramlist[["MatrixPath"]])
    } # checkAllPath end

  ) # private end

) # class end


#' Counting cut site around motif.
#' @param atacProc Do not use this parameter, we will add nore functions in the future!
#' @param CutSiteFile Your cut site infoemation file(from atacCutSitePre function) path with prefix.
#' e.g. "/your_cut_site_information_path/prefix"
#' @param MotifFile Your cut site infoemation file(from MotifScan function, set position = TRUE).
#' The first 4 columns of the motif file must be "chr start_site end_site strand".
#' @param MatrixPath The output path with a prefix, an empty folder would be great.
#' e.g. "/where_you_want_to_save_output/prefix"
#' @param motif_length Motif length.
#' @param strand_length How many bp(base pair) do you want to count up/downstream of the motif.
#' @export
atacCutSiteCount <- function(atacProc = NULL, CutSiteFile = NULL, MotifFile = NULL, MatrixPath = NULL,
                             motif_length = NULL, strand_length = NULL, FootPrint = FALSE){
  tmp <- CutSiteCountR$new(atacProc, CutSiteFile, MotifFile, MatrixPath,
                           motif_length, strand_length, FootPrint)
  tmp$process()
  return(tmp)
}
