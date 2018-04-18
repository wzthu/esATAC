setClass(Class = "RMotifScan",
         contains = "ATACProc"
)

setMethod(
    f = "initialize",
    signature = "RMotifScan",
    definition = function(.Object, atacProc, ..., peak = NULL, genome = NULL,
                          motifs = NULL, p.cutoff = NULL, scanO.dir = NULL,
                          prefix = NULL, editable = FALSE){
        .Object <- init(.Object, "RMotifScan", editable, list(arg1 = atacProc))

        # necessary parameters
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
        .Object@paramlist[["motifs"]] <- motifs
        .Object@paramlist[["motifs.len"]] <- lapply(X = .Object@paramlist[["motifs"]], FUN = ncol)
        .Object@paramlist[["p.cutoff"]] <- p.cutoff
        if(is.null(prefix)){
            .Object@paramlist[["prefix"]] <- "MotifScan"
        }else{
            .Object@paramlist[["prefix"]] <- prefix
        }
        # unnecessary parameters
        if(is.null(scanO.dir)){
            .Object@paramlist[["scanO.dir"]] <- paste(tools::file_path_sans_ext(.Object@paramlist[["peak"]]),
                                                      "_", .Object@paramlist[["prefix"]],
                                                      "_MotifScanOutput", sep = "")
            dir.create(.Object@paramlist[["scanO.dir"]])
        }else{
            .Object@paramlist[["scanO.dir"]] <- scanO.dir
        }
        .Object@paramlist[["rdsOutput"]] <- paste(
            .Object@paramlist[["scanO.dir"]],
            "/", .Object@paramlist[["prefix"]], "_",
            "RMotifScan.rds",
            sep = ""
        )

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

        peaks <- rtracklayer::import(con = .Object@paramlist[["peak"]])
        SiteSetList <- motifmatchr::matchMotifs(pwms = .Object@paramlist[["motifs"]],
                                                subject = peaks,
                                                genome = .Object@paramlist[["genome"]],
                                                out = "positions",
                                                p.cutoff = .Object@paramlist[["p.cutoff"]])

        save_info <- data.frame()
        WriteMotifOrder <- 1

        for(i in seq(length(SiteSetList))){
            motif.name <- names(SiteSetList[i])
            motif.len <- .Object@paramlist[["motifs.len"]][[motif.name]]
            motif <- SiteSetList[[motif.name]]
            motif <- sort.GenomicRanges(x = motif, ignore.strand = TRUE)
            motif <- as.data.frame(motif)
            motif <- within(motif, rm(width))

            output.path <- paste(.Object@paramlist[["scanO.dir"]],
                                 "/", .Object@paramlist[["prefix"]], "_",
                                 motif.name, sep = "")

            save_info[WriteMotifOrder, 1] <- motif.name
            save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output.path)
            save_info[WriteMotifOrder, 3] <- motif.len
            WriteMotifOrder <- WriteMotifOrder + 1

            write.table(x = motif, file = output.path, row.names = FALSE,
                        col.names = FALSE, quote = FALSE, sep = "\t",
                        append = FALSE)
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
        if(is.null(.Object@paramlist[["motifs"]])){
            stop("Parameter motifs is required!")
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


#' @name RMotifScan
#' @title Search Motif Position in Given Regions
#' @description
#' Search motif position in genome according thr given motif and peak information.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}.
#' @param peak \code{Character} scalar.
#' Input region path. UCSC bed file is recommented. Other file should be able
#' to import as \link[GenomicRanges]{GRanges} objects through
#' \link[rtracklayer]{import}.
#' @param genome BSgenome object, Default: from \code{\link{setConfigure}}.
#' @param motifs either \link[TFBSTools]{PFMatrix}, \link[TFBSTools]{PFMatrixList},
#' \link[TFBSTools]{PWMatrix}, \link[TFBSTools]{PWMatrixList}.
#' @param p.cutoff p-value cutoff for returning motifs.
#' @param scanO.dir \code{Character} scalar.
#' the output file directory. This function will use the name in motifs as
#' the file name to save the motif position information in separate files.
#' @param prefix prefix for Output file.
#' @param ... Additional arguments, currently unused.
#' @details This function scan motif position in a given genome regions.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' \dontrun{
#' library(R.utils)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' peak.path <- system.file("extdata", "Example_peak1.bed.bz2", package="esATAC")
#' peak.path <- as.vector(bunzip2(filename = peak.path, destname = file.path(getwd(), "Example_peak1.bed"), ext="bz2", FUN=bzfile, overwrite=TRUE , remove = FALSE))
#'
#' motif <- readRDS(system.file("extdata", "MotifPFM.rds", package="esATAC"))
#'
#' motifscan(peak = peak.path, genome = BSgenome.Hsapiens.UCSC.hg19, motifs = motif)
#' }
#'
#'
#' @seealso
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacCutSiteCount}}
#'
#' @importFrom rtracklayer import
#' @importFrom IRanges subsetByOverlaps
#' @importFrom R.utils getAbsolutePath
#' @importFrom motifmatchr matchMotifs
#' @importFrom GenomicRanges sort.GenomicRanges
#'

setGeneric("atacMotifScan",
           function(atacProc, peak = NULL, genome = NULL,
                    motifs = NULL, p.cutoff = 1e-6, scanO.dir = NULL,
                    prefix = NULL, ...) standardGeneric("atacMotifScan"))



#' @rdname RMotifScan
#' @aliases atacMotifScan
#' @export
setMethod(
    f = "atacMotifScan",
    signature = "ATACProc",
    function(atacProc, peak = NULL, genome = NULL,
             motifs = NULL, p.cutoff = 1e-6, scanO.dir = NULL,
             prefix = NULL, ...){
        atacproc <- new(
            "RMotifScan",
            atacProc = atacProc,
            peak = peak,
            genome = genome,
            motifs = motifs,
            p.cutoff = p.cutoff,
            scanO.dir = scanO.dir,
            prefix = prefix)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname RMotifScan
#' @aliases motifscan
#' @export
motifscan <- function(peak = NULL, genome = NULL,
                      motifs = NULL, p.cutoff = 1e-6, scanO.dir = NULL,
                      prefix = NULL, ...){
    atacproc <- new(
        "RMotifScan",
        atacProc = NULL,
        peak = peak,
        genome = genome,
        motifs = motifs,
        p.cutoff = p.cutoff,
        scanO.dir = scanO.dir,
        prefix = prefix)
    atacproc <- process(atacproc)
    invisible(atacproc)
}


