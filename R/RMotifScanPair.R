mid_interval <- function(x, inerval = 100){
    x_mid <- as.integer((start(x) + end(x)) / 2)
    start(x) <- x_mid - inerval
    end(x) <- x_mid + inerval
    return(x)
}

setClass(Class = "RMotifScanPair",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "RMotifScanPair",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        peak1 <- allparam[["peak1"]]
        peak2 <- allparam[["peak2"]]
        background <- allparam[["background"]]
        genome <- allparam[["genome"]]
        motifs <- allparam[["motifs"]]
        p.cutoff <- allparam[["p.cutoff"]]
        prefix <- allparam[["prefix"]]

        atacProc <- NULL
        if(length(prevSteps)>0){
            atacProc <- prevSteps[[1]]
        }
        
        if(!is.null(atacProc)){
            input(.Object)[["peak1"]] <- getParam(atacProc, "bedOutput1")
            input(.Object)[["peak2"]] <- getParam(atacProc, "bedOutput2")
            input(.Object)[["background"]] <- getParam(atacProc, "bedOutput")
        }else{
            input(.Object)[["peak1"]] <- peak1
            input(.Object)[["peak2"]] <- peak2
            input(.Object)[["background"]] <- background
        }

        if(!is.null(genome)){
            param(.Object)[["genome"]] <- genome
        }else{
            param(.Object)[["genome"]] <- getRefRc("bsgenome")
        }

        param(.Object)[["motifs"]] <- motifs
        param(.Object)[["motifs.len"]] <- lapply(X = param(.Object)[["motifs"]], FUN = ncol)
        param(.Object)[["p.cutoff"]] <- p.cutoff

        output(.Object)[["rdsOutput.peak1"]] <- 
            getAutoPath(.Object, input(.Object)[["peak1"]], 
                        "peak1.bed", "peak1.RMotifScanPair.rds")
        output(.Object)[["rdsOutput.peak2"]] <- 
            getAutoPath(.Object, input(.Object)[["peak2"]], 
                        "peak1.bed", "peak2.RMotifScanPair.rds")
        output(.Object)[["rdsOutput.background"]] <- 
            getAutoPath(.Object, input(.Object)[["background"]], 
                        "overlap.bed", "background.RMotifScanPair.rds")
        
        
        if(is.null(prefix)){
            peak1.prefix <- getAutoPath(.Object, input(.Object)[["peak1"]], "peak1.bed", "peak1" )
            peak2.prefix <- getAutoPath(.Object, input(.Object)[["peak2"]], "peak2.bed", "peak2" )
            background.prefix <- getAutoPath(.Object, input(.Object)[["background"]], "overlap.bed", "background" )
            output(.Object)[["prefix"]] <- c(peak1.prefix, peak2.prefix, background.prefix)
        }else{
            output(.Object)[["prefix"]] <- prefix
        }

        
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RMotifScanPair",
    definition = function(.Object,...){
        
        dir.create(output(.Object)[["prefix"]][1])
        dir.create(output(.Object)[["prefix"]][2])
        dir.create(output(.Object)[["prefix"]][3])
       
        # original peak
        case.peak <- rtracklayer::import(input(.Object)[["peak1"]], format = "bed")
        ctrl.peak <- rtracklayer::import(input(.Object)[["peak2"]], format = "bed")
        backg.peak <- rtracklayer::import(input(.Object)[["background"]], format = "bed")

        # middle interval peak
        case_mid.peak <- mid_interval(x = case.peak)
        ctrl_mid.peak <- mid_interval(x = ctrl.peak)
        backg_mid.peak <- mid_interval(x = backg.peak)

        # number of peak
        case_peak.num <- length(case_mid.peak)
        ctrl_peak.num <- length(ctrl_mid.peak)
        backg_peak.num <- length(backg_mid.peak)

        # save info
        case_save_info <- data.frame()
        ctrl_save_info <- data.frame()
        backg_save_info <- data.frame()

        # start finding motif occurence
        # for case
        SiteSetList.case <- motifmatchr::matchMotifs(pwms = param(.Object)[["motifs"]],
                                                     subject = case.peak,
                                                     genome = param(.Object)[["genome"]],
                                                     out = "positions",
                                                     p.cutoff = param(.Object)[["p.cutoff"]])
        # for ctrl
        SiteSetList.ctrl <- motifmatchr::matchMotifs(pwms = param(.Object)[["motifs"]],
                                                     subject = ctrl.peak,
                                                     genome = param(.Object)[["genome"]],
                                                     out = "positions",
                                                     p.cutoff = param(.Object)[["p.cutoff"]])
        # for background
        SiteSetList.bg <- motifmatchr::matchMotifs(pwms = param(.Object)[["motifs"]],
                                                   subject = backg.peak,
                                                   genome = param(.Object)[["genome"]],
                                                   out = "positions",
                                                   p.cutoff = param(.Object)[["p.cutoff"]])

        # length check
        Motif.caselen <- length(SiteSetList.case)
        Motif.ctrllen <- length(SiteSetList.ctrl)
        Motif.bglen <- length(SiteSetList.bg)
        if(!(Motif.caselen == Motif.ctrllen & Motif.ctrllen == Motif.bglen & Motif.caselen == Motif.bglen)){
            stop("An unexpected error occured!")
        }

        # save info
        case_save_info <- data.frame()
        ctrl_save_info <- data.frame()
        backg_save_info <- data.frame()

        WriteMotifOrder <- 1

        # save files and compute motif enrichment
        for(i in seq(Motif.caselen)){
            motif.name <- names(SiteSetList.case[i])
            motif.len <- param(.Object)[["motifs.len"]][[motif.name]]

            # for case
            motif.case <- SiteSetList.case[[motif.name]]
            motif.case <- sort.GenomicRanges(x = motif.case, ignore.strand = TRUE)
            motif.case <- as.data.frame(motif.case)
            motif.case <- within(motif.case, rm(width))
            output.case <- file.path(output(.Object)[["prefix"]][1],
                                     paste0("peak1_",motif.name))
            case_save_info[WriteMotifOrder, 1] <- motif.name
            case_save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output.case)
            case_save_info[WriteMotifOrder, 3] <- motif.len
            write.table(x = motif.case, file = output.case, row.names = FALSE,
                        col.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)

            # for ctrl
            motif.ctrl <- SiteSetList.ctrl[[motif.name]]
            motif.ctrl <- sort.GenomicRanges(x = motif.ctrl, ignore.strand = TRUE)
            motif.ctrl <- as.data.frame(motif.ctrl)
            motif.ctrl <- within(motif.ctrl, rm(width))
            output.ctrl <- file.path(output(.Object)[["prefix"]][2],
                                     paste0("peak2_",motif.name))
            ctrl_save_info[WriteMotifOrder, 1] <- motif.name
            ctrl_save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output.ctrl)
            ctrl_save_info[WriteMotifOrder, 3] <- motif.len
            write.table(x = motif.ctrl, file = output.ctrl, row.names = FALSE,
                        col.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)

            # for background
            motif.bg <- SiteSetList.bg[[motif.name]]
            motif.bg <- sort.GenomicRanges(x = motif.bg, ignore.strand = TRUE)
            motif.bg <- as.data.frame(motif.bg)
            motif.bg <- within(motif.bg, rm(width))
            output.bg <- file.path(output(.Object)[["prefix"]][3],
                                   paste0("background_",motif.name))
            backg_save_info[WriteMotifOrder, 1] <- motif.name
            backg_save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output.bg)
            backg_save_info[WriteMotifOrder, 3] <- motif.len
            write.table(x = motif.bg, file = output.bg, row.names = FALSE,
                        col.names = FALSE, quote = FALSE, sep = "\t", append = FALSE)

            # processing motif enrichment
            case_overlap <- GenomicRanges::findOverlaps(query = case_mid.peak,
                                                        subject = SiteSetList.case[[motif.name]],
                                                        ignore.strand = TRUE)
            ctrl_overlap <- GenomicRanges::findOverlaps(query = ctrl_mid.peak,
                                                        subject = SiteSetList.ctrl[[motif.name]],
                                                        ignore.strand = TRUE)
            backg_overlap <- GenomicRanges::findOverlaps(query = backg_mid.peak,
                                                         subject = SiteSetList.bg[[motif.name]],
                                                         ignore.strand = TRUE)
            case.occur <- length(unique(S4Vectors::queryHits(case_overlap)))
            ctrl.occur <- length(unique(S4Vectors::queryHits(ctrl_overlap)))
            backg.occur <- length(unique(S4Vectors::queryHits(backg_overlap)))
            case.btest <- binom.test(x = case.occur, n = case_peak.num, p = backg.occur / backg_peak.num)
            ctrl.btest <- binom.test(x = ctrl.occur, n = ctrl_peak.num, p = backg.occur / backg_peak.num)
            case_save_info[WriteMotifOrder ,4] <- case.btest$p.value
            ctrl_save_info[WriteMotifOrder, 4] <- ctrl.btest$p.value

            WriteMotifOrder <- WriteMotifOrder + 1
        }

        saveRDS(object = case_save_info, file = output(.Object)[["rdsOutput.peak1"]])
        saveRDS(object = ctrl_save_info, file = output(.Object)[["rdsOutput.peak2"]])
        saveRDS(object = backg_save_info, file = output(.Object)[["rdsOutput.background"]])

        
        
        .Object
    }
)

setMethod(
    f = "genReport",
    signature = "RMotifScanPair",
    definition = function(.Object, ...){
        report(.Object)$rdsOutput.peak1 <- readRDS(output(.Object)[["rdsOutput.peak1"]])
        report(.Object)$rdsOutput.peak2 <- readRDS(output(.Object)[["rdsOutput.peak2"]])
        report(.Object)$rdsOutput.background <- readRDS(output(.Object)[["rdsOutput.background"]])
        .Object
    }
)


#' @name RMotifScanPair
#' @title Search Motif Position in Given Regions
#' @description
#' Search motif position in genome according thr given motif and peak information.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacpeakComp}}.
#' @param peak1 peak file path.
#' @param peak2 peak file path.
#' @param background background peak file path.
#' @param genome BSgenome object, Default: from \code{\link{getRefRc}}.
#' @param motifs either\code{\link{PFMatrix}}, \code{\link{PFMatrixList}},
#' \code{\link{PWMatrix}}, \code{\link{PWMatrixList}}.
#' @param p.cutoff p-value cutoff for returning motifs.
#' @param scanO.dir \code{Character} scalar.
#' the output file directory. This function will use the name in motifs as
#' the file name to save the motif position information in separate files.
#' @param prefix prefix for Output file. Order: peak1, peak2, backgroud.
#' @param ... Additional arguments, currently unused.
#' @details This function scan motif position in a given genome regions.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#'
#' @examples
#'
#' \dontrun{
#' library(R.utils)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="esATAC")
#' p2bz <- system.file("extdata", "Example_peak2.bed.bz2", package="esATAC")
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE , remove = FALSE))
#' peak2_path <- as.vector(bunzip2(filename = p2bz,
#' destname = file.path(getwd(), "Example_peak2.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' peakcom.output <- peakcomp(bedInput1 = peak1_path, bedInput2 = peak2_path,
#' olap.rate = 0.1)
#'
#' motif <- readRDS(system.file("extdata", "MotifPFM.rds", package="esATAC"))
#' output <- atacMotifScanPair(atacProc = peakcom.output,
#' genome = BSgenome.Hsapiens.UCSC.hg19, motifs = motif)
#'}
#'
#' @seealso
#' \code{\link{atacpeakComp}}
#'
#' @importFrom motifmatchr matchMotifs
#' @importFrom GenomicRanges sort.GenomicRanges
#' @importFrom rtracklayer import
#' @importFrom IRanges subsetByOverlaps

setGeneric("atacMotifScanPair",
           function(atacProc, peak1 = NULL, peak2 = NULL,
                    background = NULL, genome = NULL, motifs = NULL,
                    p.cutoff = 0.0001, scanO.dir = NULL, prefix = NULL, ...) standardGeneric("atacMotifScanPair"))

#' @rdname RMotifScanPair
#' @aliases atacMotifScanPair
#' @export
setMethod(
    f = "atacMotifScanPair",
    signature = "ATACProc",
    definition = function(atacProc, peak1 = NULL, peak2 = NULL,
                          background = NULL, genome = NULL, motifs = NULL,
                          p.cutoff = 0.0001, scanO.dir = NULL, prefix = NULL, ...){
        allpara <- c(list(Class = "RMotifScanPair", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname RMotifScanPair
#' @aliases motifscanpair
#' @export
motifscanpair <- function(peak1 = NULL, peak2 = NULL,
                          background = NULL, genome = NULL, motifs = NULL,
                          p.cutoff = 0.0001, scanO.dir = NULL, prefix = NULL, ...){
    allpara <- c(list(Class = "RMotifScanPair", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}








