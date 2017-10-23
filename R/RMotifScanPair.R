mid_interval <- function(x, inerval = 100){
    x_mid <- as.integer((start(x) + end(x)) / 2)
    start(x) <- x_mid - inerval
    end(x) <- x_mid + inerval
    return(x)
}


RMotifScanPair <- R6::R6Class(
    classname = "RMotifScanPair",
    inherit = ATACProc,

    public = list(
        initialize = function(atacProc, peak1 = NULL, peak2 = NULL,
                              background = NULL, genome = NULL,
                              motifPWM = NULL, min.score = "85%",
                              scanO.dir = NULL, n.cores = NULL,
                              prefix = NULL, editable = FALSE){
            super$initialize("RMotifScanPair", editable, list(arg1 = atacProc))

            if(!is.null(atacProc)){
                private$paramlist[["peak1"]] <- atacProc$getParam("bedOutput")[1]
                private$paramlist[["peak2"]] <- atacProc$getParam("bedOutput")[2]
                private$paramlist[["background"]] <- atacProc$getParam("bedOutput")[3]
            }else{
                private$paramlist[["peak1"]] <- peak1
                private$paramlist[["peak2"]] <- peak2
                private$paramlist[["background"]] <- background
            }

            if(!is.null(genome)){
                private$paramlist[["genome"]] <- genome
            }else{
                private$paramlist[["genome"]] <- .obtainConfigure("bsgenome")
            }

            private$paramlist[["motifPWM"]] <- motifPWM
            private$paramlist[["motifPWM.len"]] <- lapply(X = private$paramlist[["motifPWM"]], FUN = ncol)
            private$paramlist[["min.score"]] <- min.score

            if(is.null(prefix)){
                peak1.prefix <- tools::file_path_sans_ext(base::basename(private$paramlist[["peak1"]]))
                peak2.prefix <- tools::file_path_sans_ext(base::basename(private$paramlist[["peak2"]]))
                background.prefix <- tools::file_path_sans_ext(base::basename(private$paramlist[["background"]]))
                private$paramlist[["prefix"]] <- c(peak1.prefix, peak2.prefix, background.prefix)
            }else{
                private$paramlist[["prefix"]] <- prefix
            }

            if(is.null(n.cores)){
                private$paramlist[["n.cores"]] <- .obtainConfigure("threads")
            }else{
                private$paramlist[["n.cores"]] <- n.cores
            }

            if(is.null(scanO.dir)){
                private$paramlist[["scanO.dir"]] <- dirname(private$paramlist[["peak1"]])
            }else{
                private$paramlist[["scanO.dir"]] <- scanO.dir
            }

            private$paramlist[["rdsOutput.peak1"]] <- paste(
                private$paramlist[["scanO.dir"]], "/",
                private$paramlist[["prefix"]][1],
                "_", "RMotifScanPair.rds", sep = ""
            )
            private$paramlist[["rdsOutput.peak2"]] <- paste(
                private$paramlist[["scanO.dir"]], "/",
                private$paramlist[["prefix"]][2],
                "_", "RMotifScanPair.rds", sep = ""
            )
            private$paramlist[["rdsOutput.background"]] <- paste(
                private$paramlist[["scanO.dir"]], "/",
                private$paramlist[["prefix"]][3],
                "_", "RMotifScanPair.rds", sep = ""
            )
            # parameter check
            private$paramValidation()
        } # initialization end

    ), # public end

    private = list(
        processing = function(){
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("peak1 file:%s", private$paramlist[["peak1"]]))
            private$writeLog(sprintf("peak2 file:%s", private$paramlist[["peak2"]]))
            private$writeLog(sprintf("background file:%s", private$paramlist[["background"]]))
            private$writeLog(sprintf("Output destination:%s", private$paramlist[["scanO.dir"]]))

            # original peak
            case.peak <- rtracklayer::import(private$paramlist[["peak1"]], format = "bed")
            ctrl.peak <- rtracklayer::import(private$paramlist[["peak2"]], format = "bed")
            backg.peak <- rtracklayer::import(private$paramlist[["background"]], format = "bed")
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
            n_motif <- length(private$paramlist[["motifPWM"]])
            k <- private$paramlist[["n.cores"]] * 2
            motif_in_group <- split(private$paramlist[["motifPWM"]],
                                    rep(1:ceiling(n_motif/k), each=k)[1:n_motif])
            n_group <- length(motif_in_group)

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
                    # processing motif scan
                    motif_name <- names(sitesetList[i])
                    motif_len <- private$paramlist[["motifPWM.len"]][[motif_name]]
                    # for case
                    output_data <- IRanges::subsetByOverlaps(x = sitesetList[[i]],
                                                             ranges = case.peak,
                                                             ignore.strand = TRUE)
                    output_data <- sort(x = output_data, ignore.strand = TRUE)
                    output_data <- as.data.frame(output_data)
                    output_data <- within(output_data, rm(width))
                    output_path <- paste(private$paramlist[["scanO.dir"]],
                                         "/", private$paramlist[["prefix"]][1], "_",
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
                    output_path <- paste(private$paramlist[["scanO.dir"]],
                                         "/", private$paramlist[["prefix"]][2], "_",
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
                    output_path <- paste(private$paramlist[["scanO.dir"]],
                                         "/", private$paramlist[["prefix"]][3], "_",
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

            saveRDS(object = case_save_info, file = private$paramlist[["rdsOutput.peak1"]])
            saveRDS(object = ctrl_save_info, file = private$paramlist[["rdsOutput.peak2"]])
            saveRDS(object = backg_save_info, file = private$paramlist[["rdsOutput.background"]])

        }, # processing end

        checkRequireParam = function(){
            if(is.null(private$paramlist[["peak1"]])){
                stop("Parameter atacProc or peak1 is required!")
            }
            if(is.null(private$paramlist[["peak2"]])){
                stop("Parameter atacProc or peak2 is required!")
            }
            if(is.null(private$paramlist[["background"]])){
                stop("Parameter atacProc or background is required!")
            }
            if(is.null(private$paramlist[["motifPWM"]])){
                stop("Parameter motifPWM is required!")
                if(!is.list(private$paramlist[["motifPWM"]])){
                    stop("Parameter motifPWM must be a list!")
                }
            }
        }, # checkRequireParam end

        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["peak1"]])
            private$checkFileExist(private$paramlist[["peak2"]])
            private$checkFileExist(private$paramlist[["background"]])
            private$checkPathExist(private$paramlist[["scanO.dir"]])
        }, # checkAllPath end

        getReportValImp = function(item){
            if(item == "rdsOutput.peak1"){
                fp <- readRDS(private$paramlist[["rdsOutput.peak1"]])
            }else if(item == "rdsOutput.peak2"){
                fp <- readRDS(private$paramlist[["rdsOutput.peak2"]])
            }else if(item == "rdsOutput.background"){
                fp <- readRDS(private$paramlist[["rdsOutput.background"]])
            }
            return(fp)
        },

        getReportItemsImp = function(){
            return(c("rdsOutput.peak1", "rdsOutput.peak2", "rdsOutput.background"))
        }

    ) # private end

) # R6 class end

#' @name atacMotifScanPair
#' @aliases atacMotifScanPair
#' @aliases motifscanpair
#' @title Search Motif Position in Given Regions
#' @description
#' Search motif position in given genome regions according PWM matrix.
#' @param atacProc \code{\link{ATACProc}} object scalar.
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
#' @param prefix prefix for Output file.
#' @details This function scan motif position in a given genome regions.
#' @return An invisible \code{\link{ATACProc}} object scalar for
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


#' @rdname atacMotifScanPair
#' @export
atacMotifScanPair <- function(atacProc, peak1 = NULL, peak2 = NULL, background = NULL, genome = NULL,
                              motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                              n.cores = NULL, prefix = NULL){
    tmp <- RMotifScanPair$new(atacProc, peak1, peak2, background, genome,
                          motifPWM, min.score, scanO.dir, n.cores, prefix)
    tmp$process()
    invisible(tmp)
}

#' @rdname atacMotifScanPair
#' @export
motifscanpair <- function(peak1 = NULL, peak2 = NULL, background = NULL, genome = NULL,
                          motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                          n.cores = NULL, prefix = NULL){
    tmp <- RMotifScanPair$new(atacProc = NULL, peak1, peak2, background, genome,
                          motifPWM, min.score, scanO.dir, n.cores, prefix)
    tmp$process()
    invisible(tmp)
}
