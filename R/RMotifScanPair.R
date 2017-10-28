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
                          background = NULL, genome = NULL,
                          motifPWM = NULL, min.score = "85%",
                          scanO.dir = NULL, n.cores = NULL,
                          prefix = NULL, editable = FALSE){
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

        .Object@paramlist[["motifPWM"]] <- motifPWM
        .Object@paramlist[["motifPWM.len"]] <- lapply(X = .Object@paramlist[["motifPWM"]], FUN = ncol)
        .Object@paramlist[["min.score"]] <- min.score

        if(is.null(prefix)){
            peak1.prefix <- tools::file_path_sans_ext(base::basename(.Object@paramlist[["peak1"]]))
            peak2.prefix <- tools::file_path_sans_ext(base::basename(.Object@paramlist[["peak2"]]))
            background.prefix <- tools::file_path_sans_ext(base::basename(.Object@paramlist[["background"]]))
            .Object@paramlist[["prefix"]] <- c(peak1.prefix, peak2.prefix, background.prefix)
        }else{
            .Object@paramlist[["prefix"]] <- prefix
        }

        if(is.null(n.cores)){
            .Object@paramlist[["n.cores"]] <- .obtainConfigure("threads")
        }else{
            .Object@paramlist[["n.cores"]] <- n.cores
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

        # running
        sitesetList <- list()
        n_motif <- length(.Object@paramlist[["motifPWM"]])
        k <- .Object@paramlist[["n.cores"]] * 2
        motif_in_group <- split(.Object@paramlist[["motifPWM"]],
                                rep(1:ceiling(n_motif/k), each=k)[1:n_motif])
        n_group <- length(motif_in_group)

        WriteMotifOrder <- 1
        cl <- makeCluster(.Object@paramlist[["n.cores"]])
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
                                     subject = .Object@paramlist[["genome"]],
                                     min.score = .Object@paramlist[["min.score"]],
                                     with.score = TRUE)

            for(i in seq(thisGroup.motifnum)){
                # processing motif scan
                motif_name <- names(sitesetList[i])
                motif_len <- .Object@paramlist[["motifPWM.len"]][[motif_name]]
                # for case
                output_data <- IRanges::subsetByOverlaps(x = sitesetList[[i]],
                                                         ranges = case.peak,
                                                         ignore.strand = TRUE)
                output_data <- sort(x = output_data, ignore.strand = TRUE)
                output_data <- as.data.frame(output_data)
                output_data <- within(output_data, rm(width))
                output_path <- paste(.Object@paramlist[["scanO.dir"]],
                                     "/", .Object@paramlist[["prefix"]][1], "_",
                                     motif_name, sep = "")
                case_save_info[WriteMotifOrder, 1] <- motif_name
                case_save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output_path)
                case_save_info[WriteMotifOrder, 3] <- motif_len
                write.table(x = output_data, file = output_path, row.names = FALSE,
                            col.names = FALSE, quote = FALSE)
                # for ctrl
                output_data <- IRanges::subsetByOverlaps(x = sitesetList[[i]],
                                                         ranges = ctrl.peak,
                                                         ignore.strand = TRUE)
                output_data <- sort(x = output_data, ignore.strand = TRUE)
                output_data <- as.data.frame(output_data)
                output_data <- within(output_data, rm(width))
                output_path <- paste(.Object@paramlist[["scanO.dir"]],
                                     "/", .Object@paramlist[["prefix"]][2], "_",
                                     motif_name, sep = "")
                ctrl_save_info[WriteMotifOrder, 1] <- motif_name
                ctrl_save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output_path)
                ctrl_save_info[WriteMotifOrder, 3] <- motif_len
                write.table(x = output_data, file = output_path, row.names = FALSE,
                            col.names = FALSE, quote = FALSE)
                # for olap
                output_data <- IRanges::subsetByOverlaps(x = sitesetList[[i]],
                                                         ranges = backg.peak,
                                                         ignore.strand = TRUE)
                output_data <- sort(x = output_data, ignore.strand = TRUE)
                output_data <- as.data.frame(output_data)
                output_data <- within(output_data, rm(width))
                output_path <- paste(.Object@paramlist[["scanO.dir"]],
                                     "/", .Object@paramlist[["prefix"]][3], "_",
                                     motif_name, sep = "")
                backg_save_info[WriteMotifOrder, 1] <- motif_name
                backg_save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output_path)
                backg_save_info[WriteMotifOrder, 3] <- motif_len
                write.table(x = output_data, file = output_path, row.names = FALSE,
                            col.names = FALSE, quote = FALSE)

                # processing motif enrichment
                case_overlap <- GenomicRanges::findOverlaps(query = case_mid.peak,
                                                            subject = sitesetList[[i]],
                                                            ignore.strand = TRUE)
                ctrl_overlap <- GenomicRanges::findOverlaps(query = ctrl_mid.peak,
                                                            subject = sitesetList[[i]],
                                                            ignore.strand = TRUE)
                backg_overlap <- GenomicRanges::findOverlaps(query = backg_mid.peak,
                                                             subject = sitesetList[[i]],
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
        }
        stopCluster(cl)

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

#' @title Search Motif Position in Given Regions
#' @description
#' Search motif position in given genome regions according PWM matrix.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacpeakComp}}.
#' @param peak1 peak file path.
#' @param peak2 peak file path.
#' @param background background peak file path.
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
#' # library(R.utils)
#' # library(BSgenome.Hsapiens.UCSC.hg19)
#' # p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="ATACpipe")
#' # p2bz <- system.file("extdata", "Example_peak2.bed.bz2", package="ATACpipe")
#' # peak1_path <- as.vector(bunzip2(filename = p1bz,
#' # destname = file.path(getwd(), "Example_peak1.bed"),
#' # ext="bz2", FUN=bzfile, overwrite=TRUE , remove = FALSE))
#' # peak2_path <- as.vector(bunzip2(filename = p2bz,
#' # destname = file.path(getwd(), "Example_peak2.bed"),
#' # ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' # peakcom.output <- peakcomp(bedInput1 = peak1_path, bedInput2 = peak2_path,
#' # olap.rate = 0.1)
#'
#' # pwm <- readRDS(system.file("extdata", "motifPWM.rds", package="ATACpipe"))
#' # output <- atacMotifScanPair(atacProc = peakcom.output,
#' # genome = BSgenome.Hsapiens.UCSC.hg19,
#' # motifPWM = pwm)
#'}
#'
#' @seealso
#' \code{\link{atacpeakComp}}



#' @name atacMotifScanPair
#' @export
#' @docType methods
#' @rdname atacMotifScanPair-methods
setGeneric("atacMotifScanPair",
           function(atacProc, peak1 = NULL, peak2 = NULL, background = NULL, genome = NULL,
                    motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                    n.cores = NULL, prefix = NULL, ...) standardGeneric("atacMotifScanPair"))

#' @rdname atacMotifScanPair-methods
#' @aliases atacMotifScanPair
setMethod(
    f = "atacMotifScanPair",
    signature = "ATACProc",
    definition = function(atacProc, peak1 = NULL, peak2 = NULL, background = NULL, genome = NULL,
                          motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                          n.cores = NULL, prefix = NULL, ...){
        atacproc <- new(
            "RMotifScanPair",
            atacProc = atacProc,
            peak1 = peak1,
            peak2 = peak2,
            background = background,
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

#' @rdname atacMotifScanPair-methods
#' @export
motifscanpair <- function(peak1 = NULL, peak2 = NULL, background = NULL, genome = NULL,
                          motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                          n.cores = NULL, prefix = NULL, ...){
    atacproc <- new(
        "RMotifScanPair",
        atacProc = NULL,
        peak1 = peak1,
        peak2 = peak2,
        background = background,
        genome = genome,
        motifPWM = motifPWM,
        min.score = min.score,
        scanO.dir = scanO.dir,
        n.cores = n.cores,
        prefix = prefix)
    atacproc <- process(atacproc)
    invisible(atacproc)
}








