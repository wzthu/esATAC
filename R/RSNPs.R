setClass(Class = "RSNPs",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "RSNPs",
    definition = function(.Object,atacProc,..., snp.info = NULL, region.info = NULL,
                          annoOutput = NULL, editable = FALSE){
        .Object <- init(.Object,"RSNPs",editable,list(arg1=atacProc))

        if((!is.null(atacProc)) && (class(atacProc)[1] == "PeakCallingFseq")){
            .Object@paramlist[["type"]] <- "file"
            .Object@paramlist[["region.info"]] <- atacProc$getParam("bedOutput")
            regexProcName <- sprintf("(bed|%s)", atacProc$getProcName())
        }else if((!is.null(atacProc)) && (class(atacProc)[1] == "RMotifScan")){
            .Object@paramlist[["type"]] <- "rds"
            .Object@paramlist[["region.info"]] <- atacProc$getParam("rdsOutput")
            regexProcName <- sprintf("(rds|%s)", atacProc$getProcName())
        }else{
            .Object@paramlist[["type"]] <- "file"
            .Object@paramlist[["region.info"]] <- region.info
            regexProcName <- "(bed)"
        }
        .Object@paramlist[["snp.info"]] <- snp.info

        if(is.null(annoOutput)){
            prefix <- getBasenamePrefix(.Object, .Object@paramlist[["region.info"]], regexProcName)
            .Object@paramlist[["annoOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                           paste0(prefix, ".", getProcName(.Object), ".df"))
        }else{
            name_split <- unlist(base::strsplit(x = annoOutput, split = ".", fixed = TRUE))
            suffix <- tail(name_split, 1)
            name_split <- head(name_split, -1)
            if(suffix == "df"){
                .Object@paramlist[["annoOutput"]] <- annoOutput
            }else{
                .Object@paramlist[["annoOutput"]] <- paste(annoOutput, "df", sep = ".")
            }
        }
        paramValidation(.Object)
        .Object

    }
)


setMethod(
    f = "processing",
    signature = "RSNPs",
    definition = function(.Object,...){
        .Object <- writeLog(.Object,paste0("processing file:"))
        .Object <- writeLog(.Object,sprintf("SNP source:%s",.Object@paramlist[["snp.info"]]))
        .Object <- writeLog(.Object,sprintf("Region source:%s",.Object@paramlist[["region.info"]]))
        .Object <- writeLog(.Object,sprintf("Destination:%s",.Object@paramlist[["annoOutput"]]))

        SNP_info <- read.table(file = .Object@paramlist[["snp.info"]],
                               header = FALSE)
        snp_gr <- with(SNP_info, GRanges(V1, IRanges(V2 - 1, V2)))
        if(.Object@paramlist[["type"]] == "file"){
            peak_info <- read.table(file = .Object@paramlist[["region.info"]],
                                    header = FALSE)
            peak_gr <- with(peak_info, GRanges(V1, IRanges(V2, V3)))

            overlaps <- GenomicRanges::findOverlaps(query = peak_gr,
                                                    subject = snp_gr,
                                                    ignore.strand = TRUE)
            output <- cbind(peak_info[S4Vectors::queryHits(overlaps), ],
                            SNP_info[S4Vectors::subjectHits(overlaps), ])
            write.table(x = output, file = .Object@paramlist[["annoOutput"]],
                        quote = FALSE, sep = "\t", row.names = FALSE,
                        col.names = FALSE)
        }else{
            motif_info <- readRDS(.Object@paramlist[["region.info"]])
            motif_num <- nrow(motif_info)
            names(motif_info) <- c("name", "path", "length")
            file_name.df <- data.frame()
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
                output <- cbind(motif_df[S4Vectors::queryHits(overlaps), ],
                                SNP_info[S4Vectors::subjectHits(overlaps), ])
                file_name <- file.path(dirname(.Object@paramlist[["annoOutput"]]),
                                       paste(motif_info$name[i], "_snps", sep = ""))
                write.table(x = output, file = file_name, quote = FALSE,
                            sep = "\t", row.names = FALSE, col.names = FALSE)
                file_name.df[i, 1] <- file_name
            }
            write.table(x = file_name.df, file = .Object@paramlist[["annoOutput"]],
                        quote = FALSE, sep = "\t", row.names = FALSE,
                        col.names = FALSE)
        }
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "RSNPs",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["region.info"]])){
            stop("Parameter region.info is required!")
        }
        if(is.null(.Object@paramlist[["snp.info"]])){
            stop("Parameter snp.info is required!")
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "RSNPs",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["snp.info"]])
        checkPathExist(.Object,.Object@paramlist[["annoOutput"]])
    }
)


#' @title Find whether snps are in the given regions.
#' @description
#' Find snps(user providing) in given regions.
#' This function do not consider strand.
#' @param atacProc \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacMotifScan}}.
#' If from \code{\link{atacPeakCalling}}, the output file would contain the snps
#' in given region. If from \code{\link{atacMotifScan}}, the output file would
#' contain file path to the output of every motif.
#' @param snp.info \code{Character} scalar.
#' Input snp info path. The first 3 column must be chr, position, snp_name,
#' e.g. chr13   39776775    rs7993214. Other columns could be other information
#' about snps.
#' @param region.info \code{Character} scalar.
#' Input region info path. The first 3 column must be chr, position, end.
#' @param annoOutput \code{Character} scalar.
#' Output path.
#' @return An invisible \code{\link{ATACProc}} object scalar.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="ATACpipe")
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' snps <- system.file("extdata", "snp_info", package="ATACpipe")
#' #snpanno(snp.info = snps, region.info = peak1_path)
#'
#' @seealso
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacMotifScan}}

#' @name atacSNPAnno
#' @export
#' @docType methods
#' @rdname atacSNPAnno-methods
setGeneric("atacSNPAnno",function(atacProc = NULL, snp.info = NULL, region.info = NULL,
                                  annoOutput = NULL) standardGeneric("atacSNPAnno"))

#' @rdname atacSNPAnno-methods
#' @aliases atacSNPAnno
setMethod(
    f = "atacSNPAnno",
    signature = "ATACProc",
    definition = function(atacProc = NULL, snp.info = NULL, region.info = NULL,
                          annoOutput = NULL){
        atacproc <- new(
            "RSNPs",
            atacProc = atacProc,
            snp.info = snp.info,
            region.info = region.info,
            annoOutput = annoOutput)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)
#' @rdname atacSNPAnno-methods
#' @export
snpanno <- function(snp.info, region.info = NULL, annoOutput = NULL){
    atacproc <- new(
        "RSNPs",
        atacProc = NULL,
        snp.info = snp.info,
        region.info = region.info,
        annoOutput = annoOutput)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
