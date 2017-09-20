RSNPs <- R6::R6Class(
    classname = "RSNPs",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, snp.info = NULL, region.info = NULL,
                              annoOutput = NULL, editable = FALSE){
            super$initialize("RSNPs", editable, list(arg1 = atacProc))

            # necessary parameters
            if((!is.null(atacProc)) && (class(atacProc)[1] == "PeakCallingFseq")){
                private$paramlist[["type"]] <- "file"
                private$paramlist[["region.info"]] <- atacProc$getParam("bedOutput")
                regexProcName <- sprintf("(bed|%s)", atacProc$getProcName())
            }else if((!is.null(atacProc)) && (class(atacProc)[1] == "RMotifScan")){
                private$paramlist[["type"]] <- "rds"
                private$paramlist[["region.info"]] <- atacProc$getParam("rdsOutput")
                regexProcName <- sprintf("(rds|%s)", atacProc$getProcName())
            }else{
                private$paramlist[["type"]] <- "file"
                private$paramlist[["region.info"]] <- region.info
                regexProcName <- "(bed)"
            }
            private$paramlist[["snp.info"]] <- snp.info
            # unnecessary parameters
            if(is.null(annoOutput)){
                prefix <- private$getBasenamePrefix(private$paramlist[["region.info"]], regexProcName)
                private$paramlist[["annoOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                               paste0(prefix, ".", self$getProcName(), ".df"))
            }else{
                name_split <- unlist(base::strsplit(x = annoOutput, split = ".", fixed = TRUE))
                suffix <- tail(name_split, 1)
                name_split <- head(name_split, -1)
                if(suffix == "df"){
                    private$paramlist[["annoOutput"]] <- annoOutput
                }else{
                    private$paramlist[["annoOutput"]] <- paste(annoOutput, "df", sep = ".")
                }
            }
            # parameter check
            private$paramValidation()
        } # initialization end

    ), # piblic end


    private = list(
        processing = function(){
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("SNP source:%s", private$paramlist[["snp.info"]]))
            private$writeLog(sprintf("Region source:%s", private$paramlist[["region.info"]]))
            private$writeLog(sprintf("Destination:%s", private$paramlist[["annoOutput"]]))
            SNP_info <- read.table(file = private$paramlist[["snp.info"]],
                                   header = FALSE)
            snp_gr <- with(SNP_info, GRanges(V1, IRanges(V2 - 1, V2)))
            if(private$paramlist[["type"]] == "file"){
                peak_info <- read.table(file = private$paramlist[["region.info"]],
                                        header = FALSE)
                peak_gr <- with(peak_info, GRanges(V1, IRanges(V2, V3)))
                # This function do not consider strand
                overlaps <- GenomicRanges::findOverlaps(query = peak_gr,
                                                        subject = snp_gr,
                                                        ignore.strand = TRUE)
                output <- cbind(peak_info[S4Vectors::queryHits(overlaps), ],
                                SNP_info[S4Vectors::subjectHits(overlaps), ])
                write.table(x = output, file = private$paramlist[["annoOutput"]],
                            quote = FALSE, sep = "\t", row.names = FALSE,
                            col.names = FALSE)
            }else{
                motif_info <- readRDS(private$paramlist[["region.info"]])
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
                    file_name <- file.path(dirname(private$paramlist[["annoOutput"]]),
                                           paste(motif_info$name[i], "_snps", sep = ""))
                    write.table(x = output, file = file_name, quote = FALSE,
                                sep = "\t", row.names = FALSE, col.names = FALSE)
                    file_name.df[i, 1] <- file_name
                }
                write.table(x = file_name.df, file = private$paramlist[["annoOutput"]],
                            quote = FALSE, sep = "\t", row.names = FALSE,
                            col.names = FALSE)
            }

        }, # processing end

        checkRequireParam = function(){
            if(is.null(private$paramlist[["region.info"]])){
                stop("Parameter region.info is required!")
            }
            if(is.null(private$paramlist[["snp.info"]])){
                stop("Parameter snp.info is required!")
            }
        }, # checkRequireParam end

        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["snp.info"]])
            private$checkPathExist(private$paramlist[["annoOutput"]])
        } # checkAllPath end

    ) # private end

) # class end


#' @name atacMotifScan
#' @aliases atacSNPAnno
#' @aliases snpanno
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
#' @seealso
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacMotifScan}}

#' @rdname atacSNPAnno
#' @export
atacSNPAnno <- function(atacProc = NULL, snp.info = NULL, region.info = NULL,
                    annoOutput = NULL){
    tmp <- RSNPs$new(atacProc, snp.info, region.info, annoOutput)
    tmp$process()
    invisible(tmp)
}

#' @rdname atacSNPAnno
#' @export
snpanno <- function(snp.info, region.info = NULL, annoOutput = NULL){
    tmp <- RSNPs$new(atacProc = NULL, snp.info, region.info, annoOutput)
    tmp$process()
    invisible(tmp)
}
