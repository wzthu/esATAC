setClass(Class = "RPeakComp",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "RPeakComp",
    definition = function(.Object, atacProcPeak1, atacProcPeak2, ..., bedInput1 = NULL,
                          bedInput2 = NULL, bedOutput = NULL, operation = NULL,
                          editable = FALSE){
        .Object <- init(.Object,"RPeakComp",editable,list(arg1 = atacProcPeak1, arg2 = atacProcPeak2))
        if(!is.null(atacProcPeak1)){
            .Object@paramlist[["bedInput1"]] <- atacProcPeak1$getParam("bedOutput")
        }else{
            .Object@paramlist[["bedInput1"]] <- bedInput1
        }
        if(!is.null(atacProcPeak2)){
            .Object@paramlist[["bedInput2"]] <- atacProcPeak2$getParam("bedOutput")
        }else{
            .Object@paramlist[["bedInput2"]] <- bedInput2
        }
        .Object@paramlist[["operation"]] <- operation
        regexProcName <- "(bed)"
        prefix1 <- getBasenamePrefix(.Object, .Object@paramlist[["bedInput1"]], regexProcName)
        prefix2 <- getBasenamePrefix(.Object, .Object@paramlist[["bedInput2"]], regexProcName)

        if(.Object@paramlist[["operation"]] == "overlap"){
            if(is.null(bedOutput)){
                overlap_file <- paste("overlap_", prefix1, "_", prefix2, ".bed", sep = "")
                .Object@paramlist[["bedOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                              overlap_file)
            }else{
                .Object@paramlist[["bedOutput"]] <- bedOutput
            }
        }else if(.Object@paramlist[["operation"]] == "diff"){
            if(is.null(bedOutput)){
                diff_1.name <- paste("diff_", prefix1, ".bed", sep = "")
                diff_1.path <- file.path(.obtainConfigure("tmpdir"), diff_1.name)
                diff_2.name <- paste("diff_", prefix2, ".bed", sep = "")
                diff_2.path <- file.path(.obtainConfigure("tmpdir"), diff_2.name)
                .Object@paramlist[["bedOutput"]] <- c(diff_1.path, diff_2.path)
            }else{
                .Object@paramlist[["bedOutput"]] <- bedOutput
            }
        }

        paramValidation(.Object)
        .Object

    }
)


setMethod(
    f = "processing",
    signature = "RPeakComp",
    definition = function(.Object,...){
        .Object <- writeLog(.Object,paste0("processing file:"))
        .Object <- writeLog(.Object,sprintf("source:%s",.Object@paramlist[["bedInput1"]]))
        .Object <- writeLog(.Object,sprintf("source:%s",.Object@paramlist[["bedInput2"]]))
        .Object <- writeLog(.Object,sprintf("destination:%s",.Object@paramlist[["bedOutput"]]))

        bed1 <- rtracklayer::import(con = .Object@paramlist[["bedInput1"]], format = "bed")
        bed2 <- rtracklayer::import(con = .Object@paramlist[["bedInput2"]], format = "bed")
        if(.Object@paramlist[["operation"]] == "overlap"){
            output_data <- GenomicRanges::intersect(bed1, bed2, ignore.strand = TRUE)
            rtracklayer::export(object = output_data,
                                con = .Object@paramlist[["bedOutput"]],
                                format = "bed")
        }else{
            overlaps <- GenomicRanges::findOverlaps(query = bed1, subject = bed2,
                                                    ignore.strand = TRUE)
            bed1_index <- seq(length(bed1))
            bed2_index <- seq(length(bed2))
            bed1_output <- bed1[setdiff(bed1_index, S4Vectors::queryHits(overlaps))]
            bed2_output <- bed2[setdiff(bed2_index, S4Vectors::subjectHits(overlaps))]
            rtracklayer::export(object = bed1_output,
                                con = .Object@paramlist[["bedOutput"]][1],
                                format = "bed")
            rtracklayer::export(object = bed2_output,
                                con = .Object@paramlist[["bedOutput"]][2],
                                format = "bed")
        }

        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "RPeakComp",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bedInput1"]])){
            stop("Parameter bedInput1 is requied!")
        }
        if(is.null(.Object@paramlist[["bedInput2"]])){
            stop("Parameter bedInput2 is requied!")
        }
        if(is.null(.Object@paramlist[["operation"]])){
            stop("Parameter operation is requied!")
            stopifnot(.Object@paramlist[["operation"]] %in% c("overlap", "diff"))
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "RPeakComp",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["bedInput1"]])
        checkFileExist(.Object,.Object@paramlist[["bedInput2"]])
    }
)



#' @title Find the overlap or differential peaks between two samples.
#' @description
#' This function compares two peak file and report overlap or differential peaks
#' according to the parameter "operation".
#' @param atacProcPeak1 \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}.
#' @param atacProcPeak2 \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}.
#' @param bedInput1 \code{Character} scalar.
#' Input peak file path. UCSC bed file is recommented. Other file should be
#' able to import as \link[GenomicRanges]{GRanges} objects through
#' \link[rtracklayer]{import}.
#' @param bedInput2 \code{Character} scalar.
#' Input peak file path. UCSC bed file is recommented. Other file should be
#' able to import as \link[GenomicRanges]{GRanges} objects through
#' \link[rtracklayer]{import}.
#' @param bedOutput The output file path. If "operation" is "overlap", it
#' should be a \code{Character} scalar. If "operation" is "diff", it should be
#' a \code{vector} contains 2 file path(\code{Character}).
#' @param operation "overlap" or "diff".
#' @return An invisible \code{\link{ATACProc}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="ATACpipe")
#' p2bz <- system.file("extdata", "Example_peak1.bed.bz2", package="ATACpipe")
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE , remove = FALSE))
#' peak2_path <- as.vector(bunzip2(filename = p2bz,
#' destname = file.path(getwd(), "Example_peak2.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' ## overlap peak
#' peakcomp(bedInput1 = peak1_path, bedInput2 = peak2_path,
#' operation = "overlap")
#' ## differential peak
#' peakcomp(bedInput1 = peak1_path, bedInput2 = peak2_path,
#' operation = "diff")
#'
#' @seealso
#' \code{\link{atacPeakCalling}}

#' @name atacpeakComp
#' @export
#' @docType methods
#' @rdname atacpeakComp-methods
setGeneric("atacpeakComp",function(atacProcPeak1, atacProcPeak2, bedInput1 = NULL,
                                  bedInput2 = NULL, bedOutput = NULL, operation = NULL) standardGeneric("atacpeakComp"))

#' @rdname atacpeakComp-methods
#' @aliases atacpeakComp
setMethod(
    f = "atacpeakComp",
    signature = "ATACProc",
    definition = function(atacProcPeak1, atacProcPeak2, bedInput1 = NULL,
                          bedInput2 = NULL, bedOutput = NULL, operation = NULL){
        atacproc <- new(
            "RPeakComp",
            atacProcPeak1 = atacProcPeak1,
            atacProcPeak2 = atacProcPeak2,
            bedInput1 = bedInput1,
            bedInput2 = bedInput2,
            bedOutput = bedOutput,
            operation = operation)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname atacpeakComp-methods
#' @export
peakcomp <- function(bedInput1, bedInput2, bedOutput = NULL, operation = NULL){
    atacproc <- new(
        "RPeakComp",
        atacProcPeak1 = NULL,
        atacProcPeak2 = NULL,
        bedInput1 = bedInput1,
        bedInput2 = bedInput2,
        bedOutput = bedOutput,
        operation = operation)
    atacproc <- process(atacproc)
    invisible(atacproc)
}


