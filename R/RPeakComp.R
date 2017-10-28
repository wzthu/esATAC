setClass(Class = "RPeakComp",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "RPeakComp",
    definition = function(.Object, atacProcPeak1 = NULL, atacProcPeak2 = NULL, ...,
                          bedInput1 = NULL, bedInput2 = NULL, bedOutput = NULL,
                          olap.rate = NULL, editable = FALSE){
        .Object <- init(.Object, "RPeakComp", editable, list(arg1 = atacProcPeak1, arg2 = atacProcPeak2))

        # necessary parameters
        if(!is.null(atacProcPeak1)){
            .Object@paramlist[["bedInput1"]] <- getParam(atacProcPeak1, "bedOutput")
        }else{
            .Object@paramlist[["bedInput1"]] <- bedInput1
        }
        if(!is.null(atacProcPeak2)){
            .Object@paramlist[["bedInput2"]] <- getParam(atacProcPeak2, "bedOutput")
        }else{
            .Object@paramlist[["bedInput2"]] <- bedInput2
        }
        .Object@paramlist[["olap.rate"]] <- olap.rate
        regexProcName <- "(bed)"
        .Object@paramlist[["prefix1"]] <- getBasenamePrefix(.Object, .Object@paramlist[["bedInput1"]], regexProcName)
        .Object@paramlist[["prefix2"]] <- getBasenamePrefix(.Object, .Object@paramlist[["bedInput2"]], regexProcName)

        venn_file <- paste("VennPlot_", .Object@paramlist[["prefix1"]], "_", .Object@paramlist[["prefix2"]], ".pdf", sep = "")
        .Object@paramlist[["venn.plot"]] <- file.path(.obtainConfigure("tmpdir"), venn_file)

        venn_data <- paste("VennData_", .Object@paramlist[["prefix1"]], "_", .Object@paramlist[["prefix2"]], ".rds", sep = "")
        .Object@paramlist[["venn.data"]] <- file.path(.obtainConfigure("tmpdir"), venn_data)

        # unnecessary parameters
        if(is.null(bedOutput)){
            bedInput1_specific <- paste(.Object@paramlist[["prefix1"]], "_specific.bed", sep = "")
            bedInput1_specific_path <- file.path(.obtainConfigure("tmpdir"), bedInput1_specific)
            bedInput2_specific <- paste(.Object@paramlist[["prefix2"]], "_specific.bed", sep = "")
            bedInput2_specific_path <- file.path(.obtainConfigure("tmpdir"), bedInput2_specific)
            overlap_file <- paste("overlap_", .Object@paramlist[["prefix1"]], "_", .Object@paramlist[["prefix2"]], ".bed", sep = "")
            overlap_file_path <- file.path(.obtainConfigure("tmpdir"), overlap_file)
            .Object@paramlist[["bedOutput"]] <- c(bedInput1_specific_path,
                                                  bedInput2_specific_path,
                                                  overlap_file_path)
        }else{
            .Object@paramlist[["bedOutput"]] <- bedOutput
        }
        paramValidation(.Object)
        .Object

    }
)


setMethod(
    f = "processing",
    signature = "RPeakComp",
    definition = function(.Object,...){
        .Object <- writeLog(.Object, paste0("processing file:"))
        .Object <- writeLog(.Object, sprintf("bed1 source:%s", .Object@paramlist[["bedInput1"]]))
        .Object <- writeLog(.Object, sprintf("bed2 source:%s", .Object@paramlist[["bedInput2"]]))
        .Object <- writeLog(.Object, sprintf("bed1 specific peak:%s", .Object@paramlist[["bedOutput"]][1]))
        .Object <- writeLog(.Object, sprintf("bed2 specific peak:%s", .Object@paramlist[["bedOutput"]][2]))
        .Object <- writeLog(.Object, sprintf("overlap peak:%s", .Object@paramlist[["bedOutput"]][3]))

        gr_a <- rtracklayer::import(con = .Object@paramlist[["bedInput1"]], format = "bed")
        gr_b <- rtracklayer::import(con = .Object@paramlist[["bedInput2"]], format = "bed")
        o <- GenomicRanges::findOverlaps(query = gr_a, subject = gr_b, ignore.strand = TRUE)
        o_1 <- gr_a[S4Vectors::queryHits(o)]
        o_2 <- gr_b[S4Vectors::subjectHits(o)]

        peak_intersect <- IRanges::pintersect(o_1, o_2)
        keep <- (IRanges::width(peak_intersect)/pmin(IRanges::width(o_1),IRanges::width(o_2)) >= .Object@paramlist[["olap.rate"]])

        overlap_peak <- IRanges::union(o_1[keep], o_2[keep])
        a_diff <- IRanges::setdiff(gr_a, o_1[keep])
        b_diff <- IRanges::setdiff(gr_b, o_2[keep])

        rtracklayer::export(object = a_diff, con = .Object@paramlist[["bedOutput"]][1], format = "bed")
        rtracklayer::export(object = b_diff, con = .Object@paramlist[["bedOutput"]][2], format = "bed")
        rtracklayer::export(object = overlap_peak, con = .Object@paramlist[["bedOutput"]][3], format = "bed")

        num_a <- length(a_diff)
        num_b <- length(b_diff)
        num_olap <- length(overlap_peak)
        num_gr_a <- length(gr_a)
        num_gr_b <- length(gr_b)
        saveRDS(object = c(num_a, num_b, num_olap, num_gr_a, num_gr_b), file = .Object@paramlist[["venn.data"]])

        venn.plot <- VennDiagram::draw.pairwise.venn(
            area1 = num_a + num_olap,
            area2 = num_b + num_olap,
            cross.area = num_olap,
            category = c(.Object@paramlist[["prefix1"]], .Object@paramlist[["prefix2"]]),
            fill = c("blue", "red"),
            cex = 3,
            cat.cex = 2,
            cat.dist = 0.01,
            cat.just = list(c(-1, -1), c(1, 1)),
            ext.pos = 30,
            ext.dist = -0.05,
            ext.length = 0.85,
            ext.line.lwd = 2,
            ext.line.lty = "dashed"
        )
        pdf(file = .Object@paramlist[["venn.plot"]])
        grid::grid.draw(venn.plot)
        dev.off()

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
    }
)


setMethod(
    f = "checkAllPath",
    signature = "RPeakComp",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["bedInput1"]])
        checkFileExist(.Object,.Object@paramlist[["bedInput2"]])
        checkPathExist(.Object,.Object@paramlist[["bedOutput"]][1])
        checkPathExist(.Object,.Object@paramlist[["bedOutput"]][2])
        checkPathExist(.Object,.Object@paramlist[["bedOutput"]][3])
    }
)


setMethod(
    f = "getReportValImp",
    signature = "RPeakComp",
    definition = function(.Object, item){
        if(item == "venn.data"){
            vd <- readRDS(.Object@paramlist[["venn.data"]])
            return(vd)
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "RPeakComp",
    definition = function(.Object){
        return(c("venn.data"))
    }
)



#' @title Find the overlap or differential peaks between two samples.
#' @description
#' This function compares two peak file and report overlap or differential peaks
#' according to the parameter "operation".
#' @param atacProcPeak1 \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}.
#' @param atacProcPeak2 \code{\link{ATACProc-class}} object scalar.
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
#' @param bedOutput The output file path. File name order: bedInput1 specific
#' peaks, bedInput2 specific peaks, overlap peaks.
#' @param olap.rate Overlap rate, if the overlap region between 2 peak is more
#' than this rate of the short peak, these two peak are considered to be
#' overlap and will be merged to a bigger peak. Default: 0.2. NOTICE: multi-peak will be
#' merged together!
#' @param ... Additional arguments, currently unused.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @importFrom VennDiagram draw.pairwise.venn
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @importFrom GenomicRanges findOverlaps
#' @examples
#'
#' library(R.utils)
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="esATAC")
#' p2bz <- system.file("extdata", "Example_peak1.bed.bz2", package="esATAC")
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
                                   bedInput2 = NULL, bedOutput = NULL, olap.rate = 0.2, ...) standardGeneric("atacpeakComp"))

#' @rdname atacpeakComp-methods
#' @aliases atacpeakComp
setMethod(
    f = "atacpeakComp",
    signature = "ATACProc",
    definition = function(atacProcPeak1, atacProcPeak2, bedInput1 = NULL,
                          bedInput2 = NULL, bedOutput = NULL, olap.rate = 0.2, ...){
        atacproc <- new(
            "RPeakComp",
            atacProcPeak1 = atacProcPeak1,
            atacProcPeak2 = atacProcPeak2,
            bedInput1 = bedInput1,
            bedInput2 = bedInput2,
            bedOutput = bedOutput,
            olap.rate = olap.rate)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname atacpeakComp-methods
#' @export
peakcomp <- function(bedInput1 = NULL, bedInput2 = NULL, bedOutput = NULL, olap.rate = 0.2, ...){
    atacproc <- new(
        "RPeakComp",
        atacProcPeak1 = NULL,
        atacProcPeak2 = NULL,
        bedInput1 = bedInput1,
        bedInput2 = bedInput2,
        bedOutput = bedOutput,
        olap.rate = olap.rate)
    atacproc <- process(atacproc)
    invisible(atacproc)
}


