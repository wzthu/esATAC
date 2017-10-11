RPeakAnno <- R6::R6Class(
  classname = "RPeakAnno",
  inherit = ATACProc,

  public = list(
    initialize = function(atacProc = NULL, peakInput = NULL, tssRegion = NULL, TxDb = NULL, level = NULL,
                          genomicAnnotationPriority = NULL, annoDb = NULL, addFlankGeneInfo = NULL,
                          flankDistance = NULL, sameStrand = NULL, ignoreOverlap = NULL, ignoreUpstream = NULL,
                          ignoreDownstream = NULL, overlap = NULL, annoOutput = NULL, editable = FALSE){
      super$initialize("RPeakAnno", editable, list(arg1 = atacProc))

      # necessary parameters
      if(!is.null(atacProc)){ # class from PeakCallingFseq
        private$paramlist[["peakInput"]] <- atacProc$getParam("bedOutput")
        regexProcName <- sprintf("(bed|%s)", atacProc$getProcName())
      }else{
        private$paramlist[["peakInput"]] <- peakInput
        regexProcName <- "(bed|%s)"
      }
      private$paramlist[["tssRegion"]] <- tssRegion
      if(!is.null(TxDb)){
        private$paramlist[["TxDb"]] <- TxDb
      }else{
        private$paramlist[["TxDb"]] <- .obtainConfigure("knownGene")
      }

      private$paramlist[["level"]] <- level
      private$paramlist[["genomicAnnotationPriority"]] <- genomicAnnotationPriority
      if(is.null(annoDb)){
          private$paramlist[["annoDb"]] <- .obtainConfigure("annoDb")
      }else{
          private$paramlist[["annoDb"]] <- annoDb
      }
      private$paramlist[["addFlankGeneInfo"]] <- addFlankGeneInfo
      private$paramlist[["flankDistance"]] <- flankDistance
      private$paramlist[["sameStrand"]] <- sameStrand
      private$paramlist[["ignoreOverlap"]] <- ignoreOverlap
      private$paramlist[["ignoreUpstream"]] <- ignoreUpstream
      private$paramlist[["ignoreDownstream"]] <- ignoreDownstream
      private$paramlist[["overlap"]] <- overlap

      # unnecessary parameters
      if(is.null(annoOutput)){
        prefix <- private$getBasenamePrefix(private$paramlist[["peakInput"]], regexProcName)
        annoOutput.dir <- file.path(.obtainConfigure("tmpdir"),
                                    paste0(prefix, ".", self$getProcName()))
        private$paramlist[["annoOutput.img"]] <- paste(annoOutput.dir,
                                                       ".pdf", sep = "")
        private$paramlist[["annoOutput.df"]] <- paste(annoOutput.dir,
                                                      ".df", sep = "")
        private$paramlist[["annoOutput.rds"]] <- paste(annoOutput.dir,
                                                       ".rds", sep = "")
      }else{
        name_split <- unlist(base::strsplit(x = annoOutput, split = ".", fixed = TRUE))
        suffix <- tail(name_split, 1)
        name_split <- head(name_split, -1)
        if(suffix == "df"){
          private$paramlist[["annoOutput.df"]] <- annoOutput
          private$paramlist[["annoOutput.img"]] <- paste(name_split, "pdf", sep = ".")
          private$paramlist[["annoOutput.rds"]] <- paste(name_split, "rds", sep = ".")
        }else{
          private$paramlist[["annoOutput.df"]] <- paste(annoOutput, "df", sep = ".")
          private$paramlist[["annoOutput.img"]] <- paste(annoOutput, "pdf", sep = ".")
          private$paramlist[["annoOutput.rds"]] <- paste(annoOutput, "rds", sep = ".")
        }
      }

      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end

  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s", private$paramlist[["peakInput"]]))
      private$writeLog(sprintf("dataframe destination:%s", private$paramlist[["annoOutput.df"]]))
      private$writeLog(sprintf("Image destination:%s", private$paramlist[["annoOutput.img"]]))
      peakGRange <- rtracklayer::import(con = private$paramlist[["peakInput"]], format = "bed")
      peakAn <- ChIPseeker::annotatePeak(peak = peakGRange,
                                         tssRegion = private$paramlist[["tssRegion"]],
                                         TxDb = private$paramlist[["TxDb"]],
                                         level = private$paramlist[["level"]],
                                         genomicAnnotationPriority = private$paramlist[["genomicAnnotationPriority"]],
                                         annoDb = private$paramlist[["annoDb"]],
                                         addFlankGeneInfo = private$paramlist[["addFlankGeneInfo"]],
                                         flankDistance = private$paramlist[["flankDistance"]],
                                         sameStrand = private$paramlist[["sameStrand"]],
                                         ignoreOverlap = private$paramlist[["ignoreOverlap"]],
                                         ignoreUpstream = private$paramlist[["ignoreUpstream"]],
                                         ignoreDownstream = private$paramlist[["ignoreDownstream"]],
                                         overlap = private$paramlist[["overlap"]])
      saveRDS(peakAn, private$paramlist[["annoOutput.rds"]])
      pdf(file = private$paramlist[["annoOutput.img"]])
      ChIPseeker::plotAnnoPie(x = peakAn)
      dev.off()
      tmp_file <- as.data.frame(peakAn)
      colnames(tmp_file)[1] <- "chromatin"
      write.table(x = tmp_file, file = private$paramlist[["annoOutput.df"]],
                  quote = FALSE, row.names = FALSE, sep = "\t",
                  col.names = TRUE, append = FALSE)
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["peakInput"]])){
        stop("Parameter peakInput is required!")
      }
    },

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["peakInput"]])
      private$checkPathExist(private$paramlist[["annoOutput.df"]])
      private$checkPathExist(private$paramlist[["annoOutput.img"]])
      private$checkPathExist(private$paramlist[["annoOutput.rds"]])
    }, # checkAllPath end

    getReportValImp = function(item){
      if(item == "annoOutput.rds"){
        peakAnno <- readRDS(private$paramlist[["annoOutput.rds"]])
        return(peakAnno)
      }
    },

    getReportItemsImp = function(){
      return(c("annoOutput.rds"))
    }

  ) # private end

) # R6 class end

#' @name atacPeakAnno
#' @aliases atacPeakAnno
#' @aliases peakanno
#' @title Annotate ATAC-seq Peak
#' @description
#' This function annotates ATAC-seq peak by a given annotation database.
#' For more information, please see \link[ChIPseeker]{annotatePeak}.
#' @param atacProc \code{\link{ATACProc}} object scalar.
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
#' @return An invisible \code{\link{ATACProc}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#'
#' @examples
#'
#' library(R.utils)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="ATACpipe")
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
#' @rdname atacPeakAnno
#' @export
atacPeakAnno <- function(atacProc, peakInput = NULL, tssRegion = c(-1000, 1000),
                         TxDb = NULL, level = "transcript",
                         genomicAnnotationPriority = c("Promoter", "5UTR",
                                                       "3UTR", "Exon", "Intron",
                                                       "Downstream", "Intergenic"),
                         annoDb = NULL, addFlankGeneInfo = FALSE,
                         flankDistance = 5000, sameStrand = FALSE,
                         ignoreOverlap = FALSE, ignoreUpstream = FALSE,
                         ignoreDownstream = FALSE, overlap = "TSS",
                         annoOutput = NULL){
  tmp <- RPeakAnno$new(atacProc, peakInput, tssRegion, TxDb, level,
                       genomicAnnotationPriority, annoDb, addFlankGeneInfo,
                       flankDistance, sameStrand, ignoreOverlap, ignoreUpstream,
                       ignoreDownstream, overlap, annoOutput)
  tmp$process()
  invisible(tmp)
}

#' @rdname atacPeakAnno
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
                     annoOutput = NULL){
  tmp <- RPeakAnno$new(atacProc = NULL, peakInput, tssRegion, TxDb, level,
                       genomicAnnotationPriority, annoDb, addFlankGeneInfo,
                       flankDistance, sameStrand, ignoreOverlap, ignoreUpstream,
                       ignoreDownstream, overlap, annoOutput)
  tmp$process()
  invisible(tmp)
}
