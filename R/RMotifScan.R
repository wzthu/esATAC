RMotifScan <- R6::R6Class(
  classname = "RMotifScan",
  inherit = ATACProc,

  public = list(
    initialize = function(atacProc, peak = NULL, genome = NULL,
                          motifPWM = NULL, min.score = NULL,
                          scanOutput = NULL, n.cores = NULL, editable = FALSE){
      super$initialize("RMotifScan", editable, list(arg1 = atacProc))

      # necessary parameters
      if(!is.null(atacProc)){
        private$paramlist[["peak"]] <- atacProc$getParam("bedOutput");
      }else{
        private$paramlist[["peak"]] <- peak
      }
      if(!is.null(genome)){
        private$paramlist[["genome"]] <- genome
      }else{
        private$paramlist[["genome"]] <- .obtainConfigure("bsgenome")
        print(private$paramlist[["genome"]])
      }
      private$paramlist[["motifPWM"]] <- motifPWM
      private$paramlist[["motifPWM.len"]] <- lapply(X = private$paramlist[["motifPWM"]], FUN = ncol)
      private$paramlist[["min.score"]] <- min.score

      # unnecessary parameters
      if(is.null(scanOutput)){
        private$paramlist[["scanOutput"]] <- dirname(private$paramlist[["peak"]])
      }else{
        private$paramlist[["scanOutput"]] <- scanOutput
      }
      private$paramlist[["rdsOutput"]] <- file.path(private$paramlist[["scanOutput"]], "/RMotifScan.rds")

      if(is.null(n.cores)){
        private$paramlist[["n.cores"]] <- .obtainConfigure("threads")
        print(private$paramlist[["n.cores"]])
      }else{
        private$paramlist[["n.cores"]] <- n.cores
      }

      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end

  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("peak file:%s", private$paramlist[["peak"]]))
      private$writeLog(sprintf("Output destination:%s", private$paramlist[["scanOutput"]]))
      # running
      cl <- makeCluster(private$paramlist[["n.cores"]])
      sitesetList <- parLapply(cl = cl,
                               X = private$paramlist[["motifPWM"]],
                               fun = Biostrings::matchPWM,
                               subject = private$paramlist[["genome"]],
                               min.score = private$paramlist[["min.score"]],
                               with.score = TRUE)
      stopCluster(cl)
      n_motif <- length(sitesetList)
      peak <- rtracklayer::import(private$paramlist[["peak"]])
      save_info <- data.frame()
      for(i in seq(n_motif)){
        motif_name <- names(sitesetList[i])
        output_data <- IRanges::subsetByOverlaps(x = sitesetList[[i]],
                                                 ranges = peak,
                                                 ignore.strand = TRUE)
        output_data <- sort(x = output_data, ignore.strand = TRUE)
        output_data <- as.data.frame(output_data)
        output_data <- within(output_data, rm(width))
        output_path <- file.path(private$paramlist[["scanOutput"]], motif_name)
        motif_len <- private$paramlist[["motifPWM.len"]][[motif_name]]
        save_info[i, 1] <- motif_name
        save_info[i, 2] <- R.utils::getAbsolutePath(output_path)
        save_info[i, 3] <- motif_len
        write.table(x = output_data, file = output_path, row.names = FALSE,
                    col.names = FALSE, quote = FALSE)
      }
      saveRDS(object = save_info, file = private$paramlist[["rdsOutput"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["peak"]])){
        stop("Parameter peak is required!")
      }
      if(is.null(private$paramlist[["motifPWM"]])){
        stop("Parameter motifPWM is required!")
        if(!is.list(private$paramlist[["motifPWM"]])){
          stop("Parameter motifPWM must be a list!")
        }
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["peak"]])
      private$checkPathExist(private$paramlist[["scanOutput"]])
    } # checkAllPath end

  ) # private end

) # R6 class end

#' @name atacMotifScan
#' @aliases atacMotifScan
#' @aliases motifscan
#' @title Search Motif Position in Given Regions
#' @description
#' Search motif position in given genome regions according PWM matrix.
#' @param atacProc \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}.
#' @param peak \code{Character} scalar.
#' Input region path. UCSC bed file is recommented. Other file should be able
#' to import as \link[GenomicRanges]{GRanges} objects through
#' \link[rtracklayer]{import}.
#' @param genome A DNAString object.
#' @param motifPWM \code{list} scalar. Default: from \code{\link{setConfigure}}.
#' Every element in the \code{list} contains a motif PWM matrix.
#' e.g. pwm <- list("CTCF" = CTCF_PWMmatrix)
#' @param min.score The minimum score for counting a match. Can be given as a
#' character string containing a percentage (e.g. "85%") of the highest
#' possible score or as a single number.
#' @param scanOutput \code{Character} scalar.
#' the output file directory. This function will use the index in motifPWM as
#' the file name to save the motif position information in separate files.
#' @param n.cores How many core to run this function.
#' Default: from \code{\link{setConfigure}}.
#' @details This function scan motif position in a given genome regions.
#' @return An invisible \code{\link{ATACProc}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @seealso
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacCutSiteCount}}
#' \link[Biostrings]{matchPWM}
#' \link[IRanges]{subsetByOverlaps}
#'

#' @rdname atacMotifScan
#' @export
atacMotifScan <- function(atacProc, peak = NULL, genome = NULL,
                      motifPWM = NULL, min.score = "85%", scanOutput = NULL,
                      n.cores = NULL){
  tmp <- RMotifScan$new(atacProc, peak, genome, motifPWM, min.score,
                        scanOutput, n.cores)
  tmp$process()
  invisible(tmp)
}

#' @rdname atacMotifScan
#' @export
motifscan <- function(peak, genome = NULL,
                      motifPWM = NULL, min.score = "85%", scanOutput = NULL,
                      n.cores = NULL){
  tmp <- RMotifScan$new(atacProc = NULL, peak, genome, motifPWM, min.score,
                        scanOutput, n.cores)
  tmp$process()
  invisible(tmp)
}
