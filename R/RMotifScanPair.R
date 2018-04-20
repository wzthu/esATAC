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
    f = "initialize",
    signature = "RMotifScanPair",
    definition = function(.Object, atacProc, ..., peak1 = NULL, peak2 = NULL,
                          background = NULL, genome = NULL, motifs = NULL,
                          p.cutoff = 0.0001, scanO.dir = NULL, prefix = NULL,
                          editable = FALSE){
        .Object <- init(.Object, "RMotifScanPair", editable, list(arg1 = atacProc))

        if(!is.null(atacProc)){
            .Object@paramlist[["peak1"]] <- getParam(atacProc, "bedOutput")[1]
            .Object@paramlist[["peak2"]] <- getParam(atacProc, "bedOutput")[2]
            .Object@paramlist[["background"]] <- getParam(atacProc, "bedOutput")[3]
        }else{
            .Object@paramlist[["peak1"]] <- peak1
            .Object@paramlist[["peak2"]] <- peak2
            .Object@paramlist[["background"]] <- background
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
            peak1.prefix <- tools::file_path_sans_ext(base::basename(.Object@paramlist[["peak1"]]))
            peak2.prefix <- tools::file_path_sans_ext(base::basename(.Object@paramlist[["peak2"]]))
            background.prefix <- tools::file_path_sans_ext(base::basename(.Object@paramlist[["background"]]))
            .Object@paramlist[["prefix"]] <- c(peak1.prefix, peak2.prefix, background.prefix)
        }else{
            .Object@paramlist[["prefix"]] <- prefix
        }

        if(is.null(scanO.dir)){
            pre_prefix <- paste(.Object@paramlist[["prefix"]], sep = "", collapse = "_")
            .Object@paramlist[["scanO.dir"]] <- paste(dirname(.Object@paramlist[["peak1"]]),
                                                      "/",
                                                      pre_prefix,
                                                      "_MotifScanOutput",
                                                      sep = "")
            dir.create(.Object@paramlist[["scanO.dir"]])
        }else{
            .Object@paramlist[["scanO.dir"]] <- scanO.dir
        }

        .Object@paramlist[["rdsOutput.peak1"]] <- paste(
            .Object@paramlist[["scanO.dir"]], "/",
            .Object@paramlist[["prefix"]][1],
            "_", "RMotifScanPair.rds", sep = ""
        )
        .Object@paramlist[["rdsOutput.peak2"]] <- paste(
            .Object@paramlist[["scanO.dir"]], "/",
            .Object@paramlist[["prefix"]][2],
            "_", "RMotifScanPair.rds", sep = ""
        )
        .Object@paramlist[["rdsOutput.background"]] <- paste(
            .Object@paramlist[["scanO.dir"]], "/",
            .Object@paramlist[["prefix"]][3],
            "_", "RMotifScanPair.rds", sep = ""
        )
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RMotifScanPair",
    definition = function(.Object,...){
        .Object <- writeLog(.Object, paste0("processing file:"))
        .Object <- writeLog(.Object, sprintf("peak1 file:%s", .Object@paramlist[["peak1"]]))
        .Object <- writeLog(.Object, sprintf("peak2 file:%s", .Object@paramlist[["peak2"]]))
        .Object <- writeLog(.Object, sprintf("background file:%s", .Object@paramlist[["background"]]))
        .Object <- writeLog(.Object, sprintf("Output destination:%s", .Object@paramlist[["scanO.dir"]]))

        # original peak
        case.peak <- rtracklayer::import(.Object@paramlist[["peak1"]], format = "bed")
        ctrl.peak <- rtracklayer::import(.Object@paramlist[["peak2"]], format = "bed")
        backg.peak <- rtracklayer::import(.Object@paramlist[["background"]], format = "bed")

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
        SiteSetList.case <- motifmatchr::matchMotifs(pwms = .Object@paramlist[["motifs"]],
                                                     subject = case.peak,
                                                     genome = .Object@paramlist[["genome"]],
                                                     out = "positions",
                                                     p.cutoff = .Object@paramlist[["p.cutoff"]])
        # for ctrl
        SiteSetList.ctrl <- motifmatchr::matchMotifs(pwms = .Object@paramlist[["motifs"]],
                                                     subject = ctrl.peak,
                                                     genome = .Object@paramlist[["genome"]],
                                                     out = "positions",
                                                     p.cutoff = .Object@paramlist[["p.cutoff"]])
        # for background
        SiteSetList.bg <- motifmatchr::matchMotifs(pwms = .Object@paramlist[["motifs"]],
                                                   subject = backg.peak,
                                                   genome = .Object@paramlist[["genome"]],
                                                   out = "positions",
                                                   p.cutoff = .Object@paramlist[["p.cutoff"]])

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
            motif.len <- .Object@paramlist[["motifs.len"]][[motif.name]]

            # for case
            motif.case <- SiteSetList.case[[motif.name]]
            motif.case <- sort.GenomicRanges(x = motif.case, ignore.strand = TRUE)
            motif.case <- as.data.frame(motif.case)
            motif.case <- within(motif.case, rm(width))
            output.case <- paste(.Object@paramlist[["scanO.dir"]],
                                 "/", .Object@paramlist[["prefix"]][1], "_",
                                 motif.name, sep = "")
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
            output.ctrl <- paste(.Object@paramlist[["scanO.dir"]],
                                 "/", .Object@paramlist[["prefix"]][2], "_",
                                 motif.name, sep = "")
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
            output.bg <- paste(.Object@paramlist[["scanO.dir"]],
                                 "/", .Object@paramlist[["prefix"]][3], "_",
                                 motif.name, sep = "")
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

        saveRDS(object = case_save_info, file = .Object@paramlist[["rdsOutput.peak1"]])
        saveRDS(object = ctrl_save_info, file = .Object@paramlist[["rdsOutput.peak2"]])
        saveRDS(object = backg_save_info, file = .Object@paramlist[["rdsOutput.background"]])

        .Object
    }
)



setMethod(
    f = "checkRequireParam",
    signature = "RMotifScanPair",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["peak1"]])){
            stop("Parameter atacProc or peak1 is required!")
        }
        if(is.null(.Object@paramlist[["peak2"]])){
            stop("Parameter atacProc or peak2 is required!")
        }
        if(is.null(.Object@paramlist[["background"]])){
            stop("Parameter atacProc or background is required!")
        }
        if(is.null(.Object@paramlist[["motifs"]])){
            stop("Parameter motifs is required!")
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "RMotifScanPair",
    definition = function(.Object,...){
        checkFileExist(.Object, .Object@paramlist[["peak1"]]);
        checkFileExist(.Object, .Object@paramlist[["peak2"]]);
        checkFileExist(.Object, .Object@paramlist[["background"]]);
        checkPathExist(.Object, .Object@paramlist[["scanO.dir"]]);
    }
)

setMethod(
    f = "getReportValImp",
    signature = "RMotifScanPair",
    definition = function(.Object, item){
        if(item == "rdsOutput.peak1"){
            fp <- readRDS(.Object@paramlist[["rdsOutput.peak1"]])
        }else if(item == "rdsOutput.peak2"){
            fp <- readRDS(.Object@paramlist[["rdsOutput.peak2"]])
        }else if(item == "rdsOutput.background"){
            fp <- readRDS(.Object@paramlist[["rdsOutput.background"]])
        }
        return(fp)
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "RMotifScanPair",
    definition = function(.Object){
        return(c("rdsOutput.peak1", "rdsOutput.peak2", "rdsOutput.background"))
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
#' @param genome BSgenome object, Default: from \code{\link{setConfigure}}.
#' @param motifs either \link[TFBSTools]{PFMatrix}, \link[TFBSTools]{PFMatrixList},
#' \link[TFBSTools]{PWMatrix}, \link[TFBSTools]{PWMatrixList}.
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
# \dontrun{
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
#}
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
        atacproc <- new(
            "RMotifScanPair",
            atacProc = atacProc,
            peak1 = peak1,
            peak2 = peak2,
            background = background,
            genome = genome,
            motifs = motifs,
            p.cutoff = p.cutoff,
            scanO.dir = scanO.dir,
            prefix = prefix)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname RMotifScanPair
#' @aliases motifscanpair
#' @export
motifscanpair <- function(peak1 = NULL, peak2 = NULL,
                          background = NULL, genome = NULL, motifs = NULL,
                          p.cutoff = 0.0001, scanO.dir = NULL, prefix = NULL, ...){
    atacproc <- new(
        "RMotifScanPair",
        atacProc = NULL,
        peak1 = peak1,
        peak2 = peak2,
        background = background,
        genome = genome,
        motifs = motifs,
        p.cutoff = p.cutoff,
        scanO.dir = scanO.dir,
        prefix = prefix)
    atacproc <- process(atacproc)
    invisible(atacproc)
}








