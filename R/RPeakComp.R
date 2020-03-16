setClass(Class = "RPeakComp",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "RPeakComp",
    definition = function(.Object,prevSteps = list(),... ){
        allparam <- list(...)
        bedInput1 <- allparam[["bedInput1"]]
        bedInput2 <- allparam[["bedInput2"]]
        bedOutput <- allparam[["bedOutput"]]
        olap.rate <- allparam[["olap.rate"]]

        atacProcPeak1 <- NULL
        if(length(prevSteps) > 0){
            atacProcPeak1 <- prevSteps[[1]]
        }
        
        atacProcPeak2 <- NULL
        if(length(prevSteps) > 0){
            atacProcPeak2 <- prevSteps[[1]]
        }
        
        # necessary parameters
        if(!is.null(atacProcPeak1)){
            input(.Object)[["bedInput1"]] <- getParam(atacProcPeak1, "bedOutput")
        }else{
            input(.Object)[["bedInput1"]] <- bedInput1
        }
        if(!is.null(atacProcPeak2)){
            input(.Object)[["bedInput2"]] <- getParam(atacProcPeak2, "bedOutput")
        }else{
            input(.Object)[["bedInput2"]] <- bedInput2
        }
        param(.Object)[["olap.rate"]] <- olap.rate
        regexProcName <- "(bed)"
        
        output(.Object)[["venn.plot"]] <- getStepWorkDir(.Object, filename = "venn.plot.pdf")
        output(.Object)[["venn.data"]] <-  getStepWorkDir(.Object, filename = "venn.data.rds")

        # unnecessary parameters
        output(.Object)[["bedOutput1"]] <- getAutoPath(.Object,input(.Object)[["bedInput1"]], "bed", "peak1.bed")
        output(.Object)[["bedOutput2"]] <- getAutoPath(.Object,input(.Object)[["bedInput2"]], "bed", "peak2.bed")
        if(is.null(bedOutput)){
            output(.Object)[["bedOutput"]] <- getStepWorkDir(.Object, filename = "overlap.bed")
        }else{
            output(.Object)[["bedOutput"]] <- addFileSuffix(bedOutput,"bed")
        }
        
        peak1 <- basename(output(.Object)[["bedOutput1"]])
        peak2 <- basename(output(.Object)[["bedOutput2"]])
        param(.Object)[["prefix1"]] <- substring(peak1,1,nchar(peak1)-4)
        param(.Object)[["prefix2"]] <- substring(peak2,1,nchar(peak2)-4)
        

        .Object

    }
)


setMethod(
    f = "processing",
    signature = "RPeakComp",
    definition = function(.Object,...){
       

        gr_a <- rtracklayer::import(con = input(.Object)[["bedInput1"]], format = "bed")
        gr_b <- rtracklayer::import(con = input(.Object)[["bedInput2"]], format = "bed")
        o <- GenomicRanges::findOverlaps(query = gr_a, subject = gr_b, ignore.strand = TRUE)
        o_1 <- gr_a[S4Vectors::queryHits(o)]
        o_2 <- gr_b[S4Vectors::subjectHits(o)]

        peak_intersect <- IRanges::pintersect(o_1, o_2)
        keep <- (IRanges::width(peak_intersect)/pmin(IRanges::width(o_1),IRanges::width(o_2)) >= param(.Object)[["olap.rate"]])

        overlap_peak <- IRanges::union(o_1[keep], o_2[keep])
        a_diff <- IRanges::setdiff(gr_a, o_1[keep])
        b_diff <- IRanges::setdiff(gr_b, o_2[keep])

        rtracklayer::export(object = a_diff, con = output(.Object)[["bedOutput1"]], format = "bed")
        rtracklayer::export(object = b_diff, con = output(.Object)[["bedOutput2"]], format = "bed")
        rtracklayer::export(object = overlap_peak, con = output(.Object)[["bedOutput"]], format = "bed")

        num_a <- length(a_diff)
        num_b <- length(b_diff)
        num_olap <- length(overlap_peak)
        num_gr_a <- length(gr_a)
        num_gr_b <- length(gr_b)
        saveRDS(object = c(num_a, num_b, num_olap, num_gr_a, num_gr_b), file = output(.Object)[["venn.data"]])

        venn.plot <- VennDiagram::draw.pairwise.venn(
            area1 = num_a + num_olap,
            area2 = num_b + num_olap,
            cross.area = num_olap,
            category = c(param(.Object)[["prefix1"]], param(.Object)[["prefix2"]]),
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
        pdf(file = output(.Object)[["venn.plot"]])
        grid::grid.draw(venn.plot)
        dev.off()

        
        .Object
    }
)

setMethod(
    f = "genReport",
    signature = "RPeakComp",
    definition = function(.Object, ...){
        obj <- readRDS(output(.Object)[["venn.data"]])
        report(.Object)$venn.data <- obj
        .Object
    }
)


#' @name RPeakComp
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
#' able to import as \code{\link{GRanges}} objects through
#' \code{\link{import}}.
#' @param bedInput2 \code{Character} scalar.
#' Input peak file path. UCSC bed file is recommented. Other file should be
#' able to import as \code{\link{GRanges}} objects through
#' \code{\link{import}}.
#' @param bedOutput The output file path for overlap peaks.
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
#' p2bz <- system.file("extdata", "Example_peak2.bed.bz2", package="esATAC")
#' \dontrun{
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE , remove = FALSE))
#' peak2_path <- as.vector(bunzip2(filename = p2bz,
#' destname = file.path(getwd(), "Example_peak2.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' output <- peakcomp(bedInput1 = peak1_path, bedInput2 = peak2_path,
#' olap.rate = 0.1)
#' }
#' @seealso
#' \code{\link{atacPeakCalling}}


setGeneric("atacPeakComp",function(atacProcPeak1, atacProcPeak2, bedInput1 = NULL,
                                   bedInput2 = NULL, bedOutput = NULL, olap.rate = 0.2, ...) standardGeneric("atacPeakComp"))

#' @rdname RPeakComp
#' @aliases atacPeakComp
#' @export
setMethod(
    f = "atacPeakComp",
    signature = "ATACProc",
    definition = function(atacProcPeak1, atacProcPeak2, bedInput1 = NULL,
                          bedInput2 = NULL, bedOutput = NULL, olap.rate = 0.2, ...){
        allpara <- c(list(Class = "RPeakComp", prevSteps = list(atacProcPeak1,atacProcPeak2)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname RPeakComp
#' @aliases peakcomp
#' @export
peakcomp <- function(bedInput1 = NULL, bedInput2 = NULL, bedOutput = NULL, olap.rate = 0.2, ...){
    allpara <- c(list(Class = "RPeakComp", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}


