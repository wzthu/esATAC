setClass(Class = "RMotifScan",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "RMotifScan",
    definition = function(.Object, atacProc, ..., peak = NULL, genome = NULL,
                          motifPWM = NULL, min.score = NULL,
                          scanO.dir = NULL, n.cores = NULL,
                          prefix = NULL, editable = FALSE){
        .Object <- init(.Object, "RMotifScan", editable, list(arg1 = atacProc))

        if(!is.null(atacProc)){
            .Object@paramlist[["peak"]] <- getParam(atacProc, "bedOutput");
        }else{
            .Object@paramlist[["peak"]] <- peak
        }

        if(!is.null(genome)){
            .Object@paramlist[["genome"]] <- genome
        }else{
            .Object@paramlist[["genome"]] <- .obtainConfigure("bsgenome")
        }

        .Object@paramlist[["motifPWM"]] <- motifPWM
        .Object@paramlist[["motifPWM.len"]] <- lapply(X = .Object@paramlist[["motifPWM"]], FUN = ncol)
        .Object@paramlist[["min.score"]] <- min.score

        if(is.null(prefix)){
            .Object@paramlist[["prefix"]] <- "MotifScan"
        }else{
            .Object@paramlist[["prefix"]] <- prefix
        }

        if(is.null(scanO.dir)){
            .Object@paramlist[["scanO.dir"]] <- dirname(.Object@paramlist[["peak"]])
        }else{
            .Object@paramlist[["scanO.dir"]] <- scanO.dir
        }

        .Object@paramlist[["rdsOutput"]] <- paste(
            .Object@paramlist[["scanO.dir"]],
            "/", .Object@paramlist[["prefix"]], "_",
            "RMotifScan.rds",
            sep = ""
        )

        if(is.null(n.cores)){
            .Object@paramlist[["n.cores"]] <- .obtainConfigure("threads")
        }else{
            .Object@paramlist[["n.cores"]] <- n.cores
        }

        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RMotifScan",
    definition = function(.Object,...){
        .Object <- writeLog(.Object, paste0("processing file:"))
        .Object <- writeLog(.Object, sprintf("peak file:%s", .Object@paramlist[["peak"]]))
        .Object <- writeLog(.Object, sprintf("Output destination:%s", .Object@paramlist[["scanO.dir"]]))

        # running
        cl <- makeCluster(.Object@paramlist[["n.cores"]])
        sitesetList <- parLapply(cl = cl,
                                 X = .Object@paramlist[["motifPWM"]],
                                 fun = Biostrings::matchPWM,
                                 subject = .Object@paramlist[["genome"]],
                                 min.score = .Object@paramlist[["min.score"]],
                                 with.score = TRUE)
        stopCluster(cl)
        n_motif <- length(sitesetList)
        peak <- rtracklayer::import(.Object@paramlist[["peak"]])
        save_info <- data.frame()
        for(i in seq(n_motif)){
            motif_name <- names(sitesetList[i])
            output_data <- IRanges::subsetByOverlaps(x = sitesetList[[i]],
                                                     ranges = peak,
                                                     ignore.strand = TRUE)
            output_data <- sort(x = output_data, ignore.strand = TRUE)
            output_data <- as.data.frame(output_data)
            output_data <- within(output_data, rm(width))
            output_path <- paste(.Object@paramlist[["scanO.dir"]],
                                 "/", .Object@paramlist[["prefix"]], "_",
                                 motif_name, sep = "")
            motif_len <- .Object@paramlist[["motifPWM.len"]][[motif_name]]
            save_info[i, 1] <- motif_name
            save_info[i, 2] <- R.utils::getAbsolutePath(output_path)
            save_info[i, 3] <- motif_len
            write.table(x = output_data, file = output_path, row.names = FALSE,
                        col.names = FALSE, quote = FALSE)
        }
        saveRDS(object = save_info, file = .Object@paramlist[["rdsOutput"]])
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "RMotifScan",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["peak"]])){
            stop("Parameter peak is required!")
        }
        if(is.null(.Object@paramlist[["motifPWM"]])){
            stop("Parameter motifPWM is required!")
            if(!is.list(.Object@paramlist[["motifPWM"]])){
                stop("Parameter motifPWM must be a list!")
            }
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "RMotifScan",
    definition = function(.Object,...){
        checkFileExist(.Object, .Object@paramlist[["peak"]]);
        checkPathExist(.Object, .Object@paramlist[["scanO.dir"]]);
    }
)



#' @title Search Motif Position in Given Regions
#' @description
#' Search motif position in given genome regions according PWM matrix.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
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
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(R.utils)
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="esATAC")
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN = bzfile, overwrite=TRUE, remove = FALSE))
#' pwm <- readRDS(system.file("extdata", "motifPWM.rds", package="esATAC"))
#' #motifscan(peak = peak1_path, genome = BSgenome.Hsapiens.UCSC.hg19,
#' #motifPWM = pwm, prefix = "test")
#'
#'
#' @seealso
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacCutSiteCount}}
#' \link[Biostrings]{matchPWM}
#' \link[IRanges]{subsetByOverlaps}
#' @importFrom rtracklayer import
#' @importFrom IRanges subsetByOverlaps
#' @importFrom R.utils getAbsolutePath

#' @name atacMotifScan
#' @export
#' @docType methods
#' @rdname atacMotifScan-methods
setGeneric("atacMotifScan",
           function(atacProc, peak = NULL, genome = NULL,
                    motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                    n.cores = NULL, prefix = NULL) standardGeneric("atacMotifScan"))


#' @rdname atacMotifScan-methods
#' @aliases atacMotifScan
setMethod(
    f = "atacMotifScan",
    signature = "ATACProc",
    function(atacProc,
             peak = NULL,
             genome = NULL,
             motifPWM = NULL,
             min.score = "85%",
             scanO.dir = NULL,
             n.cores = NULL,
             prefix = NULL){
        atacproc <- new(
            "RMotifScan",
            atacProc = atacProc,
            peak = peak,
            genome = genome,
            motifPWM = motifPWM,
            min.score = min.score,
            scanO.dir = scanO.dir,
            n.cores = n.cores,
            prefix = prefix)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname atacMotifScan-methods
#' @export
motifscan <- function(peak,
                      genome = NULL,
                      motifPWM = NULL,
                      min.score = "85%",
                      scanO.dir = NULL,
                      n.cores = NULL,
                      prefix = NULL){
    atacproc <- new(
        "RMotifScan",
        atacProc = NULL,
        peak = peak,
        genome = genome,
        motifPWM = motifPWM,
        min.score = min.score,
        scanO.dir = scanO.dir,
        n.cores = n.cores,
        prefix = prefix)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
