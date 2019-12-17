setClass(Class = "RPeakAnno",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "RPeakAnno",
    definition = function(.Object, prevSteps = list(), ...){
        allparam <- list(...)
        peakInput <- allparam[["peakInput"]]
        tssRegion <- allparam[["tssRegion"]]
        TxDb <- allparam[["TxDb"]]
        level <- allparam[["level"]]
        genomicAnnotationPriority <- allparam[["genomicAnnotationPriority"]]
        annoDb <- allparam[["annoDb"]]
        addFlankGeneInfo <- allparam[["addFlankGeneInfo"]]
        flankDistance <- allparam[["flankDistance"]]
        sameStrand <- allparam[["sameStrand"]]
        ignoreOverlap <- allparam[["ignoreOverlap"]]
        ignoreUpstream <- allparam[["ignoreUpstream"]]
        ignoreDownstream <- allparam[["ignoreDownstream"]]
        overlap <- allparam[["overlap"]]
        annoOutput <- allparam[["annoOutput"]]
      
        atacProc <- NULL
        if(length(prevSteps) > 0){
            atacProc <- prevSteps[[1]]
        }
        if(!is.null(atacProc)){ # class from PeakCallingFseq
            input(.Object)[["peakInput"]] <- getParam(atacProc, "bedOutput")
        }else{
            input(.Object)[["peakInput"]] <- peakInput
        }
        param(.Object)[["tssRegion"]] <- tssRegion
        if(!is.null(TxDb)){
            param(.Object)[["TxDb"]] <- TxDb
        }else{
            param(.Object)[["TxDb"]] <- getRefRc("knownGene")
        }

        param(.Object)[["level"]] <- level
        param(.Object)[["genomicAnnotationPriority"]] <- genomicAnnotationPriority
        if(is.null(annoDb)){
            param(.Object)[["annoDb"]] <- getRefRc("knownGene")
        }else{
            param(.Object)[["annoDb"]] <- annoDb
        }
        param(.Object)[["addFlankGeneInfo"]] <- addFlankGeneInfo
        param(.Object)[["flankDistance"]] <- flankDistance
        param(.Object)[["sameStrand"]] <- sameStrand
        param(.Object)[["ignoreOverlap"]] <- ignoreOverlap
        param(.Object)[["ignoreUpstream"]] <- ignoreUpstream
        param(.Object)[["ignoreDownstream"]] <- ignoreDownstream
        param(.Object)[["overlap"]] <- overlap

        # unnecessary parameters
        if(is.null(annoOutput)){
            output(.Object)[["annoOutput.pdf"]] <- getAutoPath(.Object, input(.Object)[["peakInput"]], "bed", ".pdf")
            output(.Object)[["annoOutput.txt"]] <- getAutoPath(.Object, input(.Object)[["peakInput"]], "bed", ".txt")
            output(.Object)[["annoOutput.rds"]] <- getAutoPath(.Object, input(.Object)[["peakInput"]], "bed", ".rds")
        }else{
            output(.Object)[["annoOutput.txt"]] <- addFileSuffix(annoOutput,".txt")
            output(.Object)[["annoOutput.pdf"]] <-  addFileSuffix(annoOutput,".pdf")
            output(.Object)[["annoOutput.rds"]] <-  addFileSuffix(annoOutput,".rds")
        }
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RPeakAnno",
    definition = function(.Object,...){
        peakGRange <- rtracklayer::import(con = input(.Object)[["peakInput"]], format = "bed")
        txdb<-param(.Object)[["TxDb"]]
        if(is.character(txdb)){
          library(txdb,character.only = TRUE)
          txdb <- get0(txdb)
        }
        peakAn <- ChIPseeker::annotatePeak(peak = peakGRange,
                                           tssRegion = param(.Object)[["tssRegion"]],
                                           TxDb = txdb,
                                           level = param(.Object)[["level"]],
                                           genomicAnnotationPriority = param(.Object)[["genomicAnnotationPriority"]],
                                           annoDb = param(.Object)[["annoDb"]],
                                           addFlankGeneInfo = param(.Object)[["addFlankGeneInfo"]],
                                           flankDistance = param(.Object)[["flankDistance"]],
                                           sameStrand = param(.Object)[["sameStrand"]],
                                           ignoreOverlap = param(.Object)[["ignoreOverlap"]],
                                           ignoreUpstream = param(.Object)[["ignoreUpstream"]],
                                           ignoreDownstream = param(.Object)[["ignoreDownstream"]],
                                           overlap = param(.Object)[["overlap"]])
        saveRDS(peakAn, output(.Object)[["annoOutput.rds"]])
        pdf(file = output(.Object)[["annoOutput.pdf"]])
        ChIPseeker::plotAnnoPie(x = peakAn)
        dev.off()
        tmp_file <- as.data.frame(peakAn)
        colnames(tmp_file)[1] <- "chromatin"
        write.table(x = tmp_file, file = output(.Object)[["annoOutput.txt"]],
                    quote = FALSE, row.names = FALSE, sep = "\t",
                    col.names = TRUE, append = FALSE)
        

        
        .Object
    }
)


setMethod(
  f = "genReport",
  signature = "RPeakAnno",
  definition = function(.Object, ...){
    report(.Object)$annoOutput.rds <- readRDS(output(.Object)[["annoOutput.rds"]])
    report(.Object)$annoOutput <- output(.Object)[["annoOutput.txt"]]
    .Object
  }
)

#' @name RPeakAnno
#' @title Annotate ATAC-seq Peak
#' @description
#' This function annotates ATAC-seq peak by a given annotation database.
#' For more information, please see \code{\link{annotatePeak}}.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}.
#' @param peakInput \code{Character} scalar.
#' Input peak file path. UCSC bed file is recommented. Other file should be
#' able to import as \code{\link{GRanges}} objects through
#' \code{\link{import}}.
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


#' @rdname RPeakAnno
#' @aliases atacPeakAnno
#' @export
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
        allpara <- c(list(Class = "RPeakAnno", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname RPeakAnno
#' @aliases peakanno
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
    allpara <- c(list(Class = "RPeakAnno", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
