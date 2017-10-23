RMotifScan <- R6::R6Class(
    classname = "RMotifScan",
    inherit = ATACProc,

    public = list(
        initialize = function(atacProc, peak = NULL, genome = NULL,
                              motifPWM = NULL, min.score = NULL,
                              scanO.dir = NULL, n.cores = NULL,
                              prefix = NULL, editable = FALSE){
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
            if(is.null(prefix)){
                private$paramlist[["prefix"]] <- "motifscan"
            }else{
                private$paramlist[["prefix"]] <- prefix
            }
            # unnecessary parameters
            if(is.null(scanO.dir)){
                private$paramlist[["scanO.dir"]] <- dirname(private$paramlist[["peak"]])
            }else{
                private$paramlist[["scanO.dir"]] <- scanO.dir
            }
            private$paramlist[["rdsOutput"]] <- paste(
                private$paramlist[["scanO.dir"]],
                "/", private$paramlist[["prefix"]], "_",
                "RMotifScan.rds",
                sep = ""
            )

            if(is.null(n.cores)){
                private$paramlist[["n.cores"]] <- .obtainConfigure("threads")
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
            private$writeLog(sprintf("Output destination:%s", private$paramlist[["scanO.dir"]]))
            # running
            peak <- rtracklayer::import(private$paramlist[["peak"]])
            save_info <- data.frame()


            # processing 2*n.core motifs in each turn
            k <- private$paramlist[["n.cores"]] * 2
            n_motif <- length(private$paramlist[["motifPWM"]])
            motif_in_group <- split(private$paramlist[["motifPWM"]],
                         rep(1:ceiling(n_motif/k), each=k)[1:n_motif])
            n_group <- length(motif_in_group)

            # write order(motif index) while writing save_info
            WriteMotifOrder <- 1
            cl <- makeCluster(private$paramlist[["n.cores"]])
            for(i in seq(n_group)){
                thisGroup.motif <- motif_in_group[[i]]
                thisGroup.motifname <- names(thisGroup.motif)
                thisGroup.motifnum <- length(thisGroup.motif)
                thisGroup.motifinfo <- paste("Now, processing the following motif: ",
                                             paste(thisGroup.motifname, collapse = ","),
                                             sep = "")
                print(thisGroup.motifinfo)
                sitesetList <- parLapply(cl = cl,
                                                  X = thisGroup.motif,
                                                  fun = Biostrings::matchPWM,
                                                  subject = private$paramlist[["genome"]],
                                                  min.score = private$paramlist[["min.score"]],
                                                  with.score = TRUE)

                for(i in seq(thisGroup.motifnum)){
                    motif_name <- names(sitesetList[i])
                    output_data <- IRanges::subsetByOverlaps(x = sitesetList[[i]],
                                                             ranges = peak,
                                                             ignore.strand = TRUE)
                    output_data <- sort(x = output_data, ignore.strand = TRUE)
                    output_data <- as.data.frame(output_data)
                    output_data <- within(output_data, rm(width))
                    output_path <- paste(private$paramlist[["scanO.dir"]],
                                         "/", private$paramlist[["prefix"]], "_",
                                         motif_name, sep = "")
                    motif_len <- private$paramlist[["motifPWM.len"]][[motif_name]]
                    save_info[WriteMotifOrder, 1] <- motif_name
                    save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output_path)
                    save_info[WriteMotifOrder, 3] <- motif_len
                    WriteMotifOrder <- WriteMotifOrder + 1
                    write.table(x = output_data, file = output_path, row.names = FALSE,
                                col.names = FALSE, quote = FALSE)
                }
            }
            stopCluster(cl)

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
            private$checkPathExist(private$paramlist[["scanO.dir"]])
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
#' @param scanO.dir \code{Character} scalar.
#' the output file directory. This function will use the index in motifPWM as
#' the file name to save the motif position information in separate files.
#' @param n.cores How many core to run this function.
#' Default: from \code{\link{setConfigure}}.
#' @param prefix prefix for Output file.
#' @details This function scan motif position in a given genome regions.
#' @return An invisible \code{\link{ATACProc}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' \dontrun{
#' # library(BSgenome.Hsapiens.UCSC.hg19)
#' # library(R.utils)
#' # p1bz <- system.file("extdata", "chr20_sample_peak.bed.bz2", package="ATACpipe")
#' # peak1_path <- as.vector(bunzip2(filename = p1bz,
#' # destname = file.path(getwd(), "chr20_sample_peak.bed"),
#' # ext="bz2", FUN = bzfile, overwrite=TRUE, remove = FALSE))
#' # pwm <- readRDS(system.file("extdata", "motifPWM.rds", package="ATACpipe"))
#' # motifscan(peak = peak1_path, genome = BSgenome.Hsapiens.UCSC.hg19,
#' # motifPWM = pwm, prefix = "test")
#' }
#'
#' @seealso
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacCutSiteCount}}
#' \link[Biostrings]{matchPWM}
#' \link[IRanges]{subsetByOverlaps}
#'

#' @rdname atacMotifScan
#' @export
atacMotifScan <- function(atacProc, peak = NULL, genome = NULL,
                          motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                          n.cores = NULL, prefix = NULL){
    tmp <- RMotifScan$new(atacProc, peak, genome, motifPWM, min.score,
                          scanO.dir, n.cores, prefix)
    tmp$process()
    invisible(tmp)
}

#' @rdname atacMotifScan
#' @export
motifscan <- function(peak, genome = NULL,
                      motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                      n.cores = NULL, prefix = NULL){
    tmp <- RMotifScan$new(atacProc = NULL, peak, genome, motifPWM, min.score,
                          scanO.dir, n.cores, prefix)
    tmp$process()
    invisible(tmp)
}
