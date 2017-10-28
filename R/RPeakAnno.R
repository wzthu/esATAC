setClass(Class = "RPeakAnno",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "RPeakAnno",
    definition = function(.Object, atacProc, ..., peakInput = NULL, tssRegion = NULL, TxDb = NULL, level = NULL,
                          genomicAnnotationPriority = NULL, annoDb = NULL, addFlankGeneInfo = NULL,
                          flankDistance = NULL, sameStrand = NULL, ignoreOverlap = NULL, ignoreUpstream = NULL,
                          ignoreDownstream = NULL, overlap = NULL, annoOutput = NULL, editable = FALSE){
        .Object <- init(.Object, "RPeakAnno", editable, list(arg1 = atacProc))

        if(!is.null(atacProc)){ # class from PeakCallingFseq
            .Object@paramlist[["peakInput"]] <- getParam(atacProc, "bedOutput")
            regexProcName <- sprintf("(bed|%s)", getProcName(atacProc))
        }else{
            .Object@paramlist[["peakInput"]] <- peakInput
            regexProcName <- "(bed|%s)"
        }
        .Object@paramlist[["tssRegion"]] <- tssRegion
        if(!is.null(TxDb)){
            .Object@paramlist[["TxDb"]] <- TxDb
        }else{
            .Object@paramlist[["TxDb"]] <- .obtainConfigure("knownGene")
        }

        .Object@paramlist[["level"]] <- level
        .Object@paramlist[["genomicAnnotationPriority"]] <- genomicAnnotationPriority
        if(is.null(annoDb)){
            .Object@paramlist[["annoDb"]] <- .obtainConfigure("annoDb")
        }else{
            .Object@paramlist[["annoDb"]] <- annoDb
        }
        .Object@paramlist[["addFlankGeneInfo"]] <- addFlankGeneInfo
        .Object@paramlist[["flankDistance"]] <- flankDistance
        .Object@paramlist[["sameStrand"]] <- sameStrand
        .Object@paramlist[["ignoreOverlap"]] <- ignoreOverlap
        .Object@paramlist[["ignoreUpstream"]] <- ignoreUpstream
        .Object@paramlist[["ignoreDownstream"]] <- ignoreDownstream
        .Object@paramlist[["overlap"]] <- overlap

        # unnecessary parameters
        if(is.null(annoOutput)){
            prefix <- getBasenamePrefix(.Object, .Object@paramlist[["peakInput"]], regexProcName)
            annoOutput.dir <- file.path(.obtainConfigure("tmpdir"),
                                        paste0(prefix, ".", getProcName(.Object)))
            .Object@paramlist[["annoOutput.pdf"]] <- paste(annoOutput.dir,
                                                           ".pdf", sep = "")
            .Object@paramlist[["annoOutput.df"]] <- paste(annoOutput.dir,
                                                          ".df", sep = "")
            .Object@paramlist[["annoOutput.rds"]] <- paste(annoOutput.dir,
                                                           ".rds", sep = "")
        }else{
            name_split <- unlist(base::strsplit(x = annoOutput, split = ".", fixed = TRUE))
            suffix <- tail(name_split, 1)
            name_split <- head(name_split, -1)
            if(suffix == "df"){
                .Object@paramlist[["annoOutput.df"]] <- annoOutput
                .Object@paramlist[["annoOutput.pdf"]] <- paste(name_split, "pdf", sep = ".")
                .Object@paramlist[["annoOutput.rds"]] <- paste(name_split, "rds", sep = ".")
            }else{
                .Object@paramlist[["annoOutput.df"]] <- paste(annoOutput, "df", sep = ".")
                .Object@paramlist[["annoOutput.pdf"]] <- paste(annoOutput, "pdf", sep = ".")
                .Object@paramlist[["annoOutput.rds"]] <- paste(annoOutput, "rds", sep = ".")
            }
        }

        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RPeakAnno",
    definition = function(.Object,...){
        .Object <- writeLog(.Object,paste0("processing file:"))
        .Object <- writeLog(.Object,sprintf("source:%s",.Object@paramlist[["peakInput"]]))
        .Object <- writeLog(.Object,sprintf("dataframe destination:%s",.Object@paramlist[["annoOutput.df"]]))
        .Object <- writeLog(.Object,sprintf("Image destination:%s",.Object@paramlist[["annoOutput.pdf"]]))

        peakGRange <- rtracklayer::import(con = .Object@paramlist[["peakInput"]], format = "bed")
        peakAn <- ChIPseeker::annotatePeak(peak = peakGRange,
                                           tssRegion = .Object@paramlist[["tssRegion"]],
                                           TxDb = .Object@paramlist[["TxDb"]],
                                           level = .Object@paramlist[["level"]],
                                           genomicAnnotationPriority = .Object@paramlist[["genomicAnnotationPriority"]],
                                           annoDb = .Object@paramlist[["annoDb"]],
                                           addFlankGeneInfo = .Object@paramlist[["addFlankGeneInfo"]],
                                           flankDistance = .Object@paramlist[["flankDistance"]],
                                           sameStrand = .Object@paramlist[["sameStrand"]],
                                           ignoreOverlap = .Object@paramlist[["ignoreOverlap"]],
                                           ignoreUpstream = .Object@paramlist[["ignoreUpstream"]],
                                           ignoreDownstream = .Object@paramlist[["ignoreDownstream"]],
                                           overlap = .Object@paramlist[["overlap"]])
        saveRDS(peakAn, .Object@paramlist[["annoOutput.rds"]])
        pdf(file = .Object@paramlist[["annoOutput.pdf"]])
        print(.Object@paramlist[["annoOutput.pdf"]])
        ChIPseeker::plotAnnoPie(x = peakAn)
        dev.off()
        tmp_file <- as.data.frame(peakAn)
        colnames(tmp_file)[1] <- "chromatin"
        write.table(x = tmp_file, file = .Object@paramlist[["annoOutput.df"]],
                    quote = FALSE, row.names = FALSE, sep = "\t",
                    col.names = TRUE, append = FALSE)
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "RPeakAnno",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["peakInput"]])){
            stop("peakInput is required.")
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "RPeakAnno",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["peakInput"]])
        checkPathExist(.Object,.Object@paramlist[["annoOutput.df"]])
        checkPathExist(.Object,.Object@paramlist[["annoOutput.pdf"]])
        checkPathExist(.Object,.Object@paramlist[["annoOutput.rds"]])
    }
)


setMethod(
    f = "getReportValImp",
    signature = "RPeakAnno",
    definition = function(.Object, item, ...){
        if(item == "annoOutput.rds"){
            peakAnno <- readRDS(.Object@paramlist[["annoOutput.rds"]])
            return(peakAnno)
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "RPeakAnno",
    definition = function(.Object, ...){
        return(c("annoOutput.rds"))
    }
)




#' @title Annotate ATAC-seq Peak
#' @description
#' This function annotates ATAC-seq peak by a given annotation database.
#' For more information, please see \link[ChIPseeker]{annotatePeak}.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}.
#' @param peakInput \code{Character} scalar.
#' Input peak file path. UCSC bed file is recommented. Other file should be
#' able to import as \link[GenomicRanges]{GRanges} objects through
#' \link[rtracklayer]{import}.
#' @param tssRegion Region range of TSS, default:c(-1000, 1000).
#' @param TxDb TxDb object, annotation database.
#' @param level "transcript" or "gene".
#' @param genomicAnnotationPriority genomic annotation priority.
#' @param annoDb Gene annotation database.
#' @param addFlankGeneInfo logical, add flanking gene information from
#' the peaks.
#' @param flankDistance distance of flanking sequence.
#' @param sameStrand logical, whether find nearest/overlap gene in the
#' same strand.
#' @param ignoreOverlap logical, whether ignore overlap of TSS with peak.
#' @param ignoreUpstream logical, if True only annotate gene at the 3'
#' of the peak.
#' @param ignoreDownstream logical, if True only annotate gene at the
#' 5' of the peak.
#' @param overlap one of 'TSS' or 'all', if overlap="all", then gene overlap
#' with peak will be reported as nearest gene, no matter the overlap is at TSS
#' region or not.
#' @param annoOutput \code{Character} scalar.
#' the output file path.
#' @param ... Additional arguments, currently unused.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#'
#' @examples
#'
#' library(R.utils)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="esATAC")
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' #peakanno(peakInput = peak1_path, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#' #annoDb = 'org.Hs.eg.db')
#'
#' @references Guangchuang Yu, Li-Gen Wang, Qing-Yu He. ChIPseeker: an
#' R/Bioconductor package for ChIP peak annotation, comparison and
#' visualization. Bioinformatics 2015, 31(14):2382-2383
#' @seealso
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacGOAnalysis}}
#' @importFrom ChIPseeker annotatePeak
#' @importFrom ChIPseeker plotAnnoPie


#' @name atacPeakAnno
#' @export
#' @docType methods
#' @rdname atacPeakAnno-methods
setGeneric("atacPeakAnno",function(atacProc, peakInput = NULL, tssRegion = c(-1000, 1000),
                                  TxDb = NULL, level = "transcript",
                                  genomicAnnotationPriority = c("Promoter", "5UTR",
                                                                "3UTR", "Exon", "Intron",
                                                                "Downstream", "Intergenic"),
                                  annoDb = NULL, addFlankGeneInfo = FALSE,
                                  flankDistance = 5000, sameStrand = FALSE,
                                  ignoreOverlap = FALSE, ignoreUpstream = FALSE,
                                  ignoreDownstream = FALSE, overlap = "TSS",
                                  annoOutput = NULL, ...) standardGeneric("atacPeakAnno"))


#' @rdname atacPeakAnno-methods
#' @aliases atacPeakAnno
setMethod(
    f = "atacPeakAnno",
    signature = "ATACProc",
    definition = function(atacProc, peakInput = NULL, tssRegion = c(-1000, 1000),
                          TxDb = NULL, level = "transcript",
                          genomicAnnotationPriority = c("Promoter", "5UTR",
                                                        "3UTR", "Exon", "Intron",
                                                        "Downstream", "Intergenic"),
                          annoDb = NULL, addFlankGeneInfo = FALSE,
                          flankDistance = 5000, sameStrand = FALSE,
                          ignoreOverlap = FALSE, ignoreUpstream = FALSE,
                          ignoreDownstream = FALSE, overlap = "TSS",
                          annoOutput = NULL, ...){
        atacproc <- new(
            "RPeakAnno",
            atacProc = atacProc,
            peakInput = peakInput,
            tssRegion = tssRegion,
            TxDb = TxDb,
            level = level,
            genomicAnnotationPriority = genomicAnnotationPriority,
            annoDb = annoDb,
            addFlankGeneInfo = addFlankGeneInfo,
            flankDistance = flankDistance,
            sameStrand = sameStrand,
            ignoreOverlap = ignoreOverlap,
            ignoreUpstream = ignoreUpstream,
            ignoreDownstream = ignoreDownstream,
            overlap = overlap,
            annoOutput = annoOutput)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname atacPeakAnno-methods
#' @export
peakanno <- function(peakInput, tssRegion = c(-1000, 1000),
                     TxDb = NULL, level = "transcript",
                     genomicAnnotationPriority = c("Promoter", "5UTR",
                                                   "3UTR", "Exon", "Intron",
                                                   "Downstream", "Intergenic"),
                     annoDb = NULL, addFlankGeneInfo = FALSE,
                     flankDistance = 5000, sameStrand = FALSE,
                     ignoreOverlap = FALSE, ignoreUpstream = FALSE,
                     ignoreDownstream = FALSE, overlap = "TSS",
                     annoOutput = NULL, ...){
    atacproc <- new(
        "RPeakAnno",
        atacProc = NULL,
        peakInput = peakInput,
        tssRegion = tssRegion,
        TxDb = TxDb,
        level = level,
        genomicAnnotationPriority = genomicAnnotationPriority,
        annoDb = annoDb,
        addFlankGeneInfo = addFlankGeneInfo,
        flankDistance = flankDistance,
        sameStrand = sameStrand,
        ignoreOverlap = ignoreOverlap,
        ignoreUpstream = ignoreUpstream,
        ignoreDownstream = ignoreDownstream,
        overlap = overlap,
        annoOutput = annoOutput)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
