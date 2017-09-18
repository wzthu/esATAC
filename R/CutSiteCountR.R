CutSiteCountR <- R6::R6Class(
  classname = "CutSiteCountR",
  inherit = ATACProc,
  public = list(
    initialize = function(atacProcCutSite, atacProcMotifScan, csInput = NULL,
                          motif_info = NULL, chr = c(1:22, "X", "Y"), matrixOutput = NULL,
                          strandLength = NULL, FootPrint = TRUE, editable = FALSE){
      super$initialize("CutSiteCountR", editable, list(arg1 = atacProcCutSite, arg2 = atacProcMotifScan))

      # necessary parameters
      if(!is.null(atacProcCutSite)){
        private$paramlist[["csInput"]] <- atacProcCutSite$getParam("csOutput");
      }else{
        private$paramlist[["csInput"]] <- csInput
      }
      if(!is.null(atacProcMotifScan)){
        private$paramlist[["motif_info"]] <- readRDS(atacProcMotifScan$getParam("rdsOutput"))
      }else{
        private$paramlist[["motif_info"]] <- read.table(motif_info)
      }
      private$paramlist[["chr"]] <- as.list(chr)
      private$paramlist[["strandLength"]] <- strandLength
      private$paramlist[["FootPrint"]] <- FootPrint

      if(is.null(matrixOutput)){
        private$paramlist[["matrixOutput"]] <- paste(.obtainConfigure("tmpdir"),
                                                     "/Footprint", sep = "")
      }else{
        private$paramlist[["matrixOutput"]] <- matrixOutput
      }
      dir.create(private$paramlist[["matrixOutput"]])

      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end

  private = list(
    processing = function(){
      private$writeLog(paste0("Now, start processing!"))
      motif_num <- nrow(private$paramlist[["motif_info"]])
      for(i in seq(motif_num)){
        motif_name <- private$paramlist[["motif_info"]][i,1]
        motif_file <- private$paramlist[["motif_info"]][i,2]
        motif_length <- private$paramlist[["motif_info"]][i,3]
        matrixsave.dir <- file.path(private$paramlist[["matrixOutput"]], motif_name)
        dir.create(matrixsave.dir)
        footprint.path <- file.path(.obtainConfigure("tmpdir"), paste(motif_name, ".pdf", sep = ""))
        # start!
        private$writeLog(sprintf("Start Processing %s", motif_name))
        private$writeLog(sprintf("Matrix Destination:%s", matrixsave.dir))
        private$writeLog(sprintf("Footprint PDF Destination:%s", matrixsave.dir))
        tmp_dir <- paste(tempdir(), "/", Sys.getpid(), sep="")
        # using tmp dir to save temp data
        dir.create(tmp_dir, FALSE, TRUE, "0700")
        .chr_separate_call(ReadsIfile = motif_file,
                           ReadsOpath = tmp_dir,
                           Name = "/Motif")
        motif_tmp <- paste(tmp_dir, "/Motif", sep = "")

        chr <- private$paramlist[["chr"]]
        chr_len <- length(chr)
        for(i in seq(1:chr_len)){
          echo_str <- paste("Now, processing chr", chr[[i]], "......", sep = "")
          print(echo_str)
          CutSiteInput <- paste0(private$paramlist[["csInput"]], "_chr", chr[[i]], ".cs", collapse = "")
          MotifInput <- paste0(motif_tmp, "_chr", chr[[i]], ".bed", collapse = "")
          MatrixOutput <- paste0(matrixsave.dir, "/", motif_name , "_chr", chr[[i]], ".matrix", collapse = "")
          .CutSiteCount(readsfile = CutSiteInput, motiffile = MotifInput, matrixfile = MatrixOutput,
                        motif_len = motif_length, strand_len = private$paramlist[["strandLength"]])
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
        if(private$paramlist[["FootPrint"]]){
          fp <- apply(data, 2, sum)
          pdf(file = footprint.path)
          plot(fp, type = "l", col = "blue", lwd = 2, xlab = "Relative Distance From Motif (bp)", ylab = "Cut Site Count", xaxt = "n", yaxt = "n")
          axis(1, at = seq(1, private$paramlist[["strandLength"]], len = 3),
               labels = -(private$paramlist[["strandLength"]] + 1 - seq(1, private$paramlist[["strandLength"]] + 1, len = 3)),
               padj = -1.0, tck = -0.01)
          axis(1, at = private$paramlist[["strandLength"]] + motif_length + seq(1, private$paramlist[["strandLength"]], len = 3),
               labels = seq(0, private$paramlist[["strandLength"]], len = 3),
               padj = -1.0, tck = -0.01)
          axis(2, padj = 1.0,tck = -0.02)
          abline(v = c(private$paramlist[["strandLength"]], private$paramlist[["strandLength"]] + motif_length + 1),
                 lty = 2)
          dev.off()
        }

      }


    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["csInput"]])){
        stop("Parameter csInput is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkPathExist(private$paramlist[["csInput"]])
    } # checkAllPath end

  ) # private end

) # class end


#' Counting cut site around motif.
#' @param atacProcCutSite Result from function atacCutSitePre.
#' @param atacProcMotifScan Result from function MotifScan.
#' @param csInput Your cut site information file(from atacCutSitePre function) path with prefix.
#' e.g. "/your_cut_site_information_path/prefix"
#' @param motif_info A file as the same as rds as the output as MotifScan.
#' @param chr Which chromatin the program will processing.Default:c(1:22, "X", "Y").
#' @param matrixOutput The output directory, an empty folder would be great.
#' Default:tmpdir/Footprint
#' @param strandLength How many bp(base pair) do you want to count up/downstream of the motif.
#' default:100.
#' @param FootPrint TRUE or FALSE, plot footprint or not.
#' @export
atacCutSiteCount <- function(atacProcCutSite = NULL, atacProcMotifScan = NULL, csInput = NULL,
                             motif_info = NULL, chr = c(1:22, "X", "Y"), matrixOutput = NULL,
                             strandLength = 100, FootPrint = TRUE){
  tmp <- CutSiteCountR$new(atacProcCutSite, atacProcMotifScan, csInput,
                           motif_info, chr, matrixOutput, strandLength, FootPrint)
  tmp$process()
  invisible(tmp)
}
