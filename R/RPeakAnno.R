RPeakAnno <- R6::R6Class(
  classname = "RPeakAnno",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, peakInput = NULL, tssRegion = NULL, TxDb = NULL,
                          level = NULL, assignGenomicAnnotation = NULL,
                          genomicAnnotationPriority = NULL,
                          annoDb = NULL, addFlankGeneInfo = NULL,
                          flankDistance = NULL, sameStrand = NULL, ignoreOverlap = NULL,
                          ignoreUpstream = NULL, ignoreDownstream = NULL, overlap = NULL,
                          verbose = NULL, annoOutput = NULL, editable = FALSE){
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
      private$paramlist[["TxDb"]] <- TxDb
      private$paramlist[["level"]] <- level
      private$paramlist[["assignGenomicAnnotation"]] <- assignGenomicAnnotation
      private$paramlist[["genomicAnnotationPriority"]] <- genomicAnnotationPriority
      private$paramlist[["annoDb"]] <- annoDb
      private$paramlist[["addFlankGeneInfo"]] <- addFlankGeneInfo
      private$paramlist[["flankDistance"]] <- flankDistance
      private$paramlist[["sameStrand"]] <- sameStrand
      private$paramlist[["ignoreOverlap"]] <- ignoreOverlap
      private$paramlist[["ignoreUpstream"]] <- ignoreUpstream
      private$paramlist[["ignoreDownstream"]] <- ignoreDownstream
      private$paramlist[["overlap"]] <- overlap
      private$paramlist[["verbose"]] <- verbose

      # unnecessary parameters
      if(is.null(annoOutput)){
        prefix <- private$getBasenamePrefix(private$paramlist[["peakInput"]], regexProcName)
        private$paramlist[["annoOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                          paste0(prefix, ".", self$getProcName()))
        private$paramlist[["annoOutput.pdf"]] <- paste(private$paramlist[["annoOutput"]],
                                                        ".pdf", sep = "")
        private$paramlist[["annoOutput.df"]] <- paste(private$paramlist[["annoOutput"]],
                                                        ".df", sep = "")
      }else{
        name_split <- unlist(base::strsplit(x = annoOutput, split = ".", fixed = TRUE))
        suffix <- tail(name_split, 1)
        name_split <- head(name_split, -1)
        if(suffix == "df"){
        private$paramlist[["annoOutput.df"]] <- annoOutput
        private$paramlist[["annoOutput.pdf"]] <- paste(name_split, "pdf", sep = ".")
        }else{
          private$paramlist[["annoOutput.df"]] <- paste(annoOutput, "df", sep = ".")
          private$paramlist[["annoOutput.pdf"]] <- paste(annoOutput, "pdf", sep = ".")
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
      private$writeLog(sprintf("df destination:%s", private$paramlist[["annoOutput.df"]]))
      private$writeLog(sprintf("pdf destination:%s", private$paramlist[["annoOutput.pdf"]]))
      peakGRange <- rtracklayer::import(con = private$paramlist[["peakInput"]], format = "bed")
      peakAn <- ChIPseeker::annotatePeak(peak = peakGRange,
                                         tssRegion = private$paramlist[["tssRegion"]],
                                         TxDb = private$paramlist[["TxDb"]],
                                         level = private$paramlist[["level"]],
                                         assignGenomicAnnotation = private$paramlist[["assignGenomicAnnotation"]],
                                         genomicAnnotationPriority = private$paramlist[["genomicAnnotationPriority"]],
                                         annoDb = private$paramlist[["annoDb"]],
                                         addFlankGeneInfo = private$paramlist[["addFlankGeneInfo"]],
                                         flankDistance = private$paramlist[["flankDistance"]],
                                         sameStrand = private$paramlist[["sameStrand"]],
                                         ignoreOverlap = private$paramlist[["ignoreOverlap"]],
                                         ignoreUpstream = private$paramlist[["ignoreUpstream"]],
                                         ignoreDownstream = private$paramlist[["ignoreDownstream"]],
                                         overlap = private$paramlist[["overlap"]],
                                         verbose = private$paramlist[["verbose"]])
      pdf(file = private$paramlist[["annoOutput.pdf"]])
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
      private$checkPathExist(private$paramlist[["annoOutput"]])
    } # checkAllPath end

  ) # private end

) # R6 class end

#' Using ChIPseeker to annotate peak file.
#'
#' Just for test, other parameters will be add.
#' @param Input peak file, bed format.
#' @param Output Output file, default Input_path/output.csv.
PeakAnno <- function(atacProc = NULL, peakInput = NULL, tssRegion = c(-3000, 3000), TxDb = NULL,
                     level = "transcript", assignGenomicAnnotation = TRUE,
                     genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                   "Downstream", "Intergenic"),
                     annoDb = NULL, addFlankGeneInfo = FALSE,
                     flankDistance = 5000, sameStrand = FALSE, ignoreOverlap = FALSE,
                     ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "TSS",
                     verbose = TRUE, annoOutput = NULL){
  tmp <- RPeakAnno$new(atacProc, peakInput, tssRegion, TxDb,
                       level, assignGenomicAnnotation,
                       genomicAnnotationPriority,
                       annoDb, addFlankGeneInfo,
                       flankDistance, sameStrand, ignoreOverlap,
                       ignoreUpstream, ignoreDownstream, overlap,
                       verbose, annoOutput)
  tmp$process()
  return(tmp)
}
