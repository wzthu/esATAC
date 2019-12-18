setClass(Class = "RSNPs",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "RSNPs",
    definition = function(.Object,prevSteps = list(), ...){
        allparam <- list(...)
        snp.info <- allparam[["snp.info"]]
        region.info <- allparam[["region.info"]]
        annoOutput <- allparam[["annoOutput"]]
        
        withend <- allparam[["withend"]]

        if(!is.null(withend)){
            param(.Object)[["withend"]] <- TRUE
        }else{
            param(.Object)[["withend"]] <- FALSE
        }
        
        
        atacProc <- NULL
        if(length(prevSteps) > 0){
            atacProc <- prevSteps[[1]]
        }

        if((!is.null(atacProc)) && (stepType(atacProc) == "PeakCallingFseq")){
            param(.Object)[["type"]] <- "file"
            input(.Object)[["region.info"]] <- output(atacProc)[["bedOutput"]]
        }else if((!is.null(atacProc)) && (stepType(atacProc) == "RMotifScan")){
            param(.Object)[["type"]] <- "rds"
            input(.Object)[["region.info"]] <- output(atacProc)[["rdsOutput"]]
        }else{
            param(.Object)[["type"]] <- "file"
            input(.Object)[["region.info"]] <- region.info
        }
        if(!is.null(snp.info)){
            input(.Object)[["snp.info"]] <- snp.info
        }else{
            input(.Object)[["snp.info"]] <- getRefFiles("SNP")
            param(.Object)[["withend"]] <- TRUE
        }


        if(is.null(annoOutput)){
            output(.Object)[["annoOutput"]] <- 
                getAutoPath(.Object, input(.Object)[["region.info"]],"bed|rds", ".txt")
               
        }else{
            output(.Object)[["annoOutput"]] <- addFileSuffix(annoOutput,"txt")
        }

        .Object

    }
)


setMethod(
    f = "processing",
    signature = "RSNPs",
    definition = function(.Object,...){
        SNP_info <- read.delim(file = input(.Object)[["snp.info"]],
                               header = FALSE)
        if(param(.Object)[["withend"]]){
            snp_gr <- with(SNP_info, GRanges(SNP_info[, 1], IRanges(SNP_info[, 2] - 1, SNP_info[, 3])))
        }else{
            snp_gr <- with(SNP_info, GRanges(SNP_info[, 1], IRanges(SNP_info[, 2] - 1, SNP_info[, 2])))
        }

        if(param(.Object)[["type"]] == "file"){
            peak_info <- read.table(file = input(.Object)[["region.info"]],
                                    header = FALSE)
            peak_gr <- with(peak_info, GRanges(peak_info[, 1], IRanges(peak_info[, 2], peak_info[, 3])))

            overlaps <- GenomicRanges::findOverlaps(query = peak_gr,
                                                    subject = snp_gr,
                                                    ignore.strand = TRUE)
            output0 <- cbind(peak_info[S4Vectors::queryHits(overlaps), ],
                            SNP_info[S4Vectors::subjectHits(overlaps), ])
            write.table(x = output0, file = output(.Object)[["annoOutput"]],
                        quote = FALSE, sep = "\t", row.names = FALSE,
                        col.names = FALSE)
        }else{
            motif_info <- readRDS(input(.Object)[["region.info"]])
            motif_num <- nrow(motif_info)
            names(motif_info) <- c("name", "path", "length")
            file_name.txt <- data.frame()
            for(i in seq(motif_num)){
                motif_df <- read.table(motif_info$path[i],
                                       col.names = c("chr", "start", "end", "strand", "score", "sequence"))
                peak_gr <- GRanges(seqnames = motif_df$chr,
                                   ranges = IRanges(start = motif_df$start, end = motif_df$end),
                                   strand = motif_df$strand,
                                   score = motif_df$score,
                                   seq = motif_df$sequence)
                overlaps <- GenomicRanges::findOverlaps(query = peak_gr,
                                                        subject = snp_gr,
                                                        ignore.strand = TRUE)
                output0 <- cbind(motif_df[S4Vectors::queryHits(overlaps), ],
                                SNP_info[S4Vectors::subjectHits(overlaps), ])
                file_name <- file.path(dirname(output(.Object)[["annoOutput"]]),
                                       paste(motif_info$name[i], "_snps", sep = ""))
                write.table(x = output0, file = file_name, quote = FALSE,
                            sep = "\t", row.names = FALSE, col.names = FALSE)
                file_name.txt[i, 1] <- file_name
            }
            write.table(x = file_name.txt, file = output(.Object)[["annoOutput"]],
                        quote = FALSE, sep = "\t", row.names = FALSE,
                        col.names = FALSE)
        }
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "RSNPs",
    definition = function(.Object, ...){
        .Object
    }
)

#' @name RSNPs
#' @title Find whether snps are in the given regions.
#' @description
#' Find snps(user providing) in given regions.
#' This function do not consider strand.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacMotifScan}}.
#' If from \code{\link{atacPeakCalling}}, the output file would contain the snps
#' in given region. If from \code{\link{atacMotifScan}}, the output file would
#' contain file path to the output of every motif.
#' @param snp.info \code{Character} scalar.
#' Input snp info path. There are two type of input files(you can specify by
#' parameter withend).
#' 1.The first 2 column must be chr, position.
#' e.g. chr13   39776775    rs7993214.
#' Other columns could be other information about snps.
#' 2.The first 3 column must be chr, start, end.
#' e.g. chr13   39776775    39776775    rs7993214.
#' Other columns could be other information about snps.
#' When genome is hg19, using human disease as default.
#' @param region.info \code{Character} scalar.
#' Input region info path. The first 3 column must be chr, position, end. The
#' standard BED format is recommended.
#' @param annoOutput \code{Character} scalar.
#' Output path.
#' @param ... withend Your snp data has only one position column or 2.
#' @return An invisible \code{\link{ATACProc-class}} object scalar.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="esATAC")
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' snps <- system.file("extdata", "snp_info", package="esATAC")
#' #snpanno(snp.info = snps, region.info = peak1_path)
#'
#' @seealso
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacMotifScan}}
#' @importFrom  utils read.delim

setGeneric("atacSNPAnno",function(atacProc, snp.info = NULL, region.info = NULL,
                                  annoOutput = NULL, ...) standardGeneric("atacSNPAnno"))

#' @rdname RSNPs
#' @aliases atacSNPAnno
#' @export
setMethod(
    f = "atacSNPAnno",
    signature = "ATACProc",
    definition = function(atacProc, snp.info = NULL, region.info = NULL,
                          annoOutput = NULL, ...){
        allpara <- c(list(Class = "RSNPs", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname RSNPs
#' @aliases snpanno
#' @export
snpanno <- function(snp.info = NULL, region.info = NULL, annoOutput = NULL, ...){
    allpara <- c(list(Class = "RSNPs", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
