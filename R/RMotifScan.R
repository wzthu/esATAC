setClass(Class = "RMotifScan",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "RMotifScan",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        peak <- allparam[["peak"]]
        genome <- allparam[["genome"]]
        motifs <- allparam[["motifs"]]
        p.cutoff <- allparam[["p.cutoff"]]
        scanO.dir <- allparam[["scanO.dir"]]
        prefix <- allparam[["prefix"]]

        
        atacProc <- NULL
        if(length(prevSteps) > 0){
            atacProc <- prevSteps[[1]]
        }
        
        
        # necessary parameters
        if(!is.null(atacProc)){
            input(.Object)[["peak"]] <- getParam(atacProc, "bedOutput");
        }else{
            input(.Object)[["peak"]] <- peak
        }
        if(!is.null(genome)){
            param(.Object)[["genome"]] <- genome
        }else{
            param(.Object)[["genome"]] <- getRefRc("bsgenome")
        }
        param(.Object)[["motifs"]] <- motifs
        param(.Object)[["motifs.len"]] <- lapply(X = param(.Object)[["motifs"]], FUN = ncol)
        param(.Object)[["p.cutoff"]] <- p.cutoff
        if(is.null(prefix)){
            param(.Object)[["prefix"]] <- "MotifScan"
        }else{
            param(.Object)[["prefix"]] <- prefix
        }
        # unnecessary parameters
        if(is.null(scanO.dir)){
            output(.Object)[["scanO.dir"]] <- getAutoPath(.Object, input(.Object)[["peak"]],"bed", "MotifScanOutput")
        }else{
            output(.Object)[["scanO.dir"]] <- scanO.dir
        }
        output(.Object)[["rdsOutput"]] <- getAutoPath(.Object, input(.Object)[["peak"]],"bed", "RMotifScan.rds")
        

        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RMotifScan",
    definition = function(.Object,...){
        dir.create(output(.Object)[["scanO.dir"]])

        peaks <- rtracklayer::import(con = input(.Object)[["peak"]])
        SiteSetList <- motifmatchr::matchMotifs(pwms = param(.Object)[["motifs"]],
                                                subject = peaks,
                                                genome = param(.Object)[["genome"]],
                                                out = "positions",
                                                p.cutoff = param(.Object)[["p.cutoff"]])

        save_info <- data.frame()
        WriteMotifOrder <- 1

        for(i in seq(length(SiteSetList))){
            motif.name <- names(SiteSetList[i])
            motif.len <- param(.Object)[["motifs.len"]][[motif.name]]
            motif <- SiteSetList[[motif.name]]
            motif <- sort.GenomicRanges(x = motif, ignore.strand = TRUE)
            motif <- as.data.frame(motif)
            motif <- within(motif, rm(width))

            output.path <- paste0(output(.Object)[["scanO.dir"]],
                                 "/", param(.Object)[["prefix"]], "_",
                                 motif.name)

            save_info[WriteMotifOrder, 1] <- motif.name
            save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output.path)
            save_info[WriteMotifOrder, 3] <- motif.len
            WriteMotifOrder <- WriteMotifOrder + 1

            write.table(x = motif, file = output.path, row.names = FALSE,
                        col.names = FALSE, quote = FALSE, sep = "\t",
                        append = FALSE)
        }
        saveRDS(object = save_info, file = output(.Object)[["rdsOutput"]])

        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "RMotifScan",
    definition = function(.Object, ...){
        .Object
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
#' to import as \code{\link{GRanges}} objects through
#' \code{\link{import}}.
#' @param genome BSgenome object, Default: from \code{\link{getRefRc}}.
#' @param motifs either\code{\link{PFMatrix}}, \code{\link{PFMatrixList}},
#' \code{\link{PWMatrix}}, \code{\link{PWMatrixList}}.
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
        allpara <- c(list(Class = "RMotifScan", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname RMotifScan
#' @aliases motifscan
#' @export
motifscan <- function(peak = NULL, genome = NULL,
                      motifs = NULL, p.cutoff = 1e-6, scanO.dir = NULL,
                      prefix = NULL, ...){
    allpara <- c(list(Class = "RMotifScan", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}


