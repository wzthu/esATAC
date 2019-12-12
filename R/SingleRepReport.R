
setClass(Class = "SingleRepReport",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "SingleRepReport",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        htmlOutput <- allparam[["htmlOutput"]]
        
        
        if(is.null(htmlOutput)){
            output(.Object)$htmlOutput <- getStepWorkDir(.Object, filename = "Report.html")
        }else{
            output(.Object)$htmlOutput <- getStepWorkDir(.Object, filename = htmlOutput)
        }
        
        .Object
    }
)

#' @importFrom rmarkdown render
#' @importFrom pipeFrame stepType

setMethod(
    f = "processing",
    signature = "SingleRepReport",
    definition = function(.Object, ...){
        htmlOutput <- output(.Object)[["htmlOutput"]]
        prevSteps <- list(...)[["prevSteps"]]
        prevStepsType <- lapply(prevSteps, function(step){
            return(stepType(step))
        })
        names(prevSteps) <- unlist(prevStepsType)
        
        wholesummary <- NULL
        filtstat <- NULL
        
        
        unzipAndMerge <- prevSteps[["UnzipAndMerge"]]
        renamer <- prevSteps[["Renamer"]]
        removeAdapter <- prevSteps[["RemoveAdapter"]]
        bowtie2Mapping <- prevSteps[["Bowtie2Mapping"]]
        libComplexQC <- prevSteps[["LibComplexQC"]]
        sam2Bed <- prevSteps[["SamToBed"]]
        bedToBigWig <- prevSteps[["BedToBigWig"]]
        tssqc100 <- prevSteps[["TSSQCNFR"]]
        tssqc180_247 <- prevSteps[["TSSQCneucleosome"]]
        fragLenDistr <- prevSteps[["FragLenDistr"]]
        peakCalling <- prevSteps[["PeakCallingFseq"]]
        DHSQC <- prevSteps[["PeakQCDHS"]]
        blacklistQC <- prevSteps[["PeakQCblacklist"]]
        fripQC <- prevSteps[["FRiPQC"]]
        shortBed <- prevSteps[["BedUtils"]]
        Peakanno <- prevSteps[["RPeakAnno"]]
        goAna <- prevSteps[["RGo"]]
        output_motifscan <- prevSteps[["RMotifScan"]]
        cs_output <- prevSteps[["CutSitePre"]]
        footprint <- prevSteps[["CutSiteCountR"]]
        atacQC <- prevSteps[["Renamer"]]
        
        
        if(property(.Object)$singleEnd && !property(.Object)$interleave){
            wholesummary <- data.frame(Item=c("Sequence Files Type",
                                             "Original total reads",
                                             "-- Reads after adapter removing (ratio)",
                                             "-- -- Total mapped reads (ratio of original reads)",
                                             "-- -- -- Unique locations mapped uniquely by reads",
                                             "-- -- -- Uniquely mappable reads",
                                             "-- -- -- Non-Redundant Fraction (NRF)",
                                             "-- -- -- Locations with only 1 reads mapping uniquely",
                                             "-- -- -- Locations with only 2 reads mapping uniquely",
                                             "-- -- -- PCR Bottlenecking Coefficients 1 (PBC1)",
                                             "-- -- -- PCR Bottlenecking Coefficients 2 (PBC2)",
                                             "-- -- -- Non-mitochondrial reads (ratio)",
                                             "-- -- -- -- Unique mapped reads (ratio)",
                                             "-- -- -- -- -- Duplicate removed reads (final for use)",
                                             "-- -- -- -- -- -- -- Total peaks",
                                             "-- -- -- -- -- -- -- Peaks overlaped with union DHS ratio",
                                             "-- -- -- -- -- -- -- Peaks overlaped with blacklist ratio",
                                             "-- -- -- -- -- -- Fraction of reads in peaks (FRiP)"),
                                      
                                      Value=c(seqtype,
                                              getVMShow(getReportVal(removeAdapter,"statisticslist")[[1]]),
                                              getVMRShow(as.integer(getReportVal(removeAdapter,"statisticslist")[["Number of retained reads"]])/frag,
                                                         getReportVal(removeAdapter,"statisticslist")[[1]]),
                                              getVMRShow(getReportVal(sam2Bed,"total"),
                                                         getReportVal(removeAdapter,"statisticslist")[[1]]),
                                              getVMShow(getReportVal(libComplexQC,"total")),
                                              #sam2Bed$getReportVal("non-mitochondrial-multimap")),
                                              getVMShow(getReportVal(libComplexQC,"nonMultimap")),
                                              #getf(libComplexQC$getReportVal("NRF")),
                                              getRshow(getReportVal(libComplexQC,"total"),
                                                       getReportVal(libComplexQC,"nonMultimap")),
                                              getVMShow(getReportVal(libComplexQC,"one")),
                                              #sam2Bed$getReportVal("non-mitochondrial-multimap")),
                                              getVMShow(getReportVal(libComplexQC,"two")),
                                              #sam2Bed$getReportVal("non-mitochondrial-multimap")),
                                              #getf(libComplexQC$getReportVal("PBC1")),
                                              getRshow(getReportVal(libComplexQC,"one"),
                                                       getReportVal(libComplexQC,"total")),
                                              #getf(libComplexQC$getReportVal("PBC2")),
                                              getRshow(getReportVal(libComplexQC,"one"),
                                                       getReportVal(libComplexQC,"two")),
                                              getVMRShow(getReportVal(sam2Bed,"non-mitochondrial"),
                                                         getReportVal(sam2Bed,"total")),
                                              getVMRShow(getReportVal(sam2Bed,"non-mitochondrial-multimap"),
                                                         getReportVal(sam2Bed,"total")),
                                              getVMRShow(getReportVal(sam2Bed,"save"),
                                                         getReportVal(sam2Bed,"total")),
                                              sprintf("%d",as.numeric(getReportVal(fripQC,"totalPeaks"))),
                                              ifelse(is.null(DHSQC),"NA",getPer(getReportVal(DHSQC,"qcbedRate"))),
                                              ifelse(is.null(blacklistQC),"NA",getPer(getReportVal(blacklistQC,"qcbedRate"))),
                                              getPer(getReportVal(fripQC,"FRiP")))
                                      ,
                                      `Reference`=c("SE / PE",
                                                    "",
                                                    ">99%",
                                                    ">95%",
                                                    "",
                                                    "",
                                                    ">0.7",
                                                    "",
                                                    "",
                                                    ">0.7",
                                                    ">3",
                                                    ">70%",
                                                    "",
                                                    ">25M",
                                                    "",
                                                    "",
                                                    "",
                                                    ""
                                      )
                                      #`Annotation`=c()
            )
            filtstat = data.frame(
                Item=c("Original total reads",
                       "Reads after adapter removing (ratio)",
                       "Total mapped reads (ratio of original reads)",
                       "-- Non-mitochondrial reads (ratio)",
                       "-- -- Unique mapped reads (ratio)",
                       "-- -- -- Duplicate removed reads (ratio final for use)"
                ),
                Value=c(getVMShow(getReportVal(removeAdapter,"statisticslist")[[1]],TRUE),
                        getVMRShow(as.integer(getReportVal(removeAdapter,"statisticslist")[["Number of retained reads"]])/frag,
                                   getReportVal(removeAdapter,"statisticslist")[[1]],TRUE),
                        getVMRShow(getReportVal(sam2Bed,"total"),
                                   getReportVal(removeAdapter,"statisticslist")[[1]],TRUE),
                        getVMRShow(getReportVal(sam2Bed,"non-mitochondrial"),
                                   getReportVal(sam2Bed,"total"),TRUE),
                        getVMRShow(getReportVal(sam2Bed,"non-mitochondrial-multimap"),
                                   getReportVal(sam2Bed,"total"),TRUE),
                        getVMRShow(getReportVal(sam2Bed,"save"),
                                   getReportVal(sam2Bed,"total"),TRUE)
                ),
                `Reference`=c("",
                              ">99%",
                              ">95%",
                              ">70%",
                              ">60%",
                              ">25M,>60%"
                )
            )
        }else{
            wholesummary = data.frame(Item=c("Sequence Files Type",
                                             "Original total reads",
                                             "-- Reads after adapter removing (ratio)",
                                             "-- -- Total mapped reads (ratio of original reads)",
                                             "-- -- -- Unique locations mapped uniquely by reads",
                                             "-- -- -- Uniquely mappable reads",
                                             "-- -- -- Non-Redundant Fraction (NRF)",
                                             "-- -- -- Locations with only 1 reads mapping uniquely",
                                             "-- -- -- Locations with only 2 reads mapping uniquely",
                                             "-- -- -- PCR Bottlenecking Coefficients 1 (PBC1)",
                                             "-- -- -- PCR Bottlenecking Coefficients 2 (PBC2)",
                                             "-- -- -- Non-mitochondrial reads (ratio)",
                                             "-- -- -- -- Unique mapped reads (ratio)",
                                             "-- -- -- -- -- Duplicate removed reads (final for use)",
                                             #"---- -- -- -- -- Transcription start site (TSS) enrichment",
                                             "-- -- -- -- -- -- Nucleosome free reads (<100bp)",
                                             "-- -- -- -- -- -- -- Total peaks",
                                             "-- -- -- -- -- -- -- Peaks overlaped with union DHS ratio",
                                             "-- -- -- -- -- -- -- Peaks overlaped with blacklist ratio",
                                             "-- -- -- -- -- -- Fraction of reads in peaks (FRiP)"),
                                      
                                      Value=c(seqtype,
                                              getVMShow(getReportVal(removeAdapter,"statisticslist")[[1]]),
                                              getVMRShow(as.integer(getReportVal(removeAdapter,"statisticslist")[["Number of retained reads"]])/frag,
                                                         getReportVal(removeAdapter,"statisticslist")[[1]]),
                                              getVMRShow(getReportVal(sam2Bed,"total"),
                                                         getReportVal(removeAdapter,"statisticslist")[[1]]),
                                              getVMShow(getReportVal(libComplexQC,"total")),
                                              #sam2Bed$getReportVal("non-mitochondrial-multimap")),
                                              getVMShow(getReportVal(libComplexQC,"nonMultimap")),
                                              #getf(libComplexQC$getReportVal("NRF")),
                                              getRshow(getReportVal(libComplexQC,"total"),
                                                       getReportVal(libComplexQC,"nonMultimap")),
                                              getVMShow(getReportVal(libComplexQC,"one")),
                                              #sam2Bed$getReportVal("non-mitochondrial-multimap")),
                                              getVMShow(getReportVal(libComplexQC,"two")),
                                              #sam2Bed$getReportVal("non-mitochondrial-multimap")),
                                              #getf(libComplexQC$getReportVal("PBC1")),
                                              getRshow(getReportVal(libComplexQC,"one"),
                                                       getReportVal(libComplexQC,"total")),
                                              #getf(libComplexQC$getReportVal("PBC2")),
                                              getRshow(getReportVal(libComplexQC,"one"),
                                                       getReportVal(libComplexQC,"two")),
                                              getVMRShow(getReportVal(sam2Bed,"non-mitochondrial"),
                                                         getReportVal(sam2Bed,"total")),
                                              getVMRShow(getReportVal(sam2Bed,"non-mitochondrial-multimap"),
                                                         getReportVal(sam2Bed,"total")),
                                              getVMRShow(getReportVal(sam2Bed,"save"),
                                                         getReportVal(sam2Bed,"total")),
                                              #"",
                                              getVMRShow(getReportVal(shortBed,"save"),
                                                         getReportVal(shortBed,"total")),
                                              sprintf("%d",as.numeric(getReportVal(fripQC,"totalPeaks"))),
                                              ifelse(is.null(DHSQC),"NA",getPer(getReportVal(DHSQC,"qcbedRate"))),
                                              ifelse(is.null(blacklistQC),"NA",getPer(getReportVal(blacklistQC,"qcbedRate"))),
                                              getPer(getReportVal(fripQC,"FRiP")))
                                      ,
                                      `Reference`=c("SE / PE",
                                                    "",
                                                    ">99%",
                                                    ">95%",
                                                    "",
                                                    "",
                                                    ">0.7",
                                                    "",
                                                    "",
                                                    ">0.7",
                                                    ">3",
                                                    ">70%",
                                                    "",
                                                    ">25M",
                                                    #"",
                                                    "",
                                                    "",
                                                    "",
                                                    "",
                                                    ""
                                      )
                                      #`Annotation`=c()
            )
            filtstat = data.frame(
                Item=c("Original total reads",
                       "Reads after adapter removing (ratio)",
                       "Total mapped reads (ratio of original reads)",
                       "-- Non-mitochondrial reads (ratio)",
                       "-- -- Unique mapped reads (ratio)",
                       "-- -- -- Duplicate removed reads (ratio final for use)"
                ),
                Value=c(getVMShow(getReportVal(removeAdapter,"statisticslist")[[1]],TRUE),
                        getVMRShow(as.integer(getReportVal(removeAdapter,"statisticslist")[["Number of retained reads"]])/frag,
                                   getReportVal(removeAdapter,"statisticslist")[[1]],TRUE),
                        getVMRShow(getReportVal(sam2Bed,"total"),
                                   getReportVal(removeAdapter,"statisticslist")[[1]],TRUE),
                        getVMRShow(getReportVal(sam2Bed,"non-mitochondrial"),
                                   getReportVal(sam2Bed,"total"),TRUE),
                        getVMRShow(getReportVal(sam2Bed,"non-mitochondrial-multimap"),
                                   getReportVal(sam2Bed,"total"),TRUE),
                        getVMRShow(getReportVal(sam2Bed,"save"),
                                   getReportVal(sam2Bed,"total"),TRUE)
                ),
                `Reference`=c("",
                              ">99%",
                              ">95%",
                              ">70%",
                              ">60%",
                              ">25M,>60%"
                )
            )
        }
        
        
        
        save(prevSteps, wholesummary, filtstat, file = getStepWorkDir(.Object = .Object, filename = "PrevSteps.Rdata"))
        
        reportmkd <- getStepWorkDir(.Object = .Object, filename = "Report.Rmd")
        
        reportmkd1 <- getStepWorkDir(.Object = .Object, filename = "Report.code.Rmd")
        
        file.copy(from = system.file(package = "esATAC", "extdata","Report.Rmd"),
                  to = reportmkd,overwrite = TRUE)
        
  #      file.copy(from = system.file(package = "enrichTF", "extdata","Report.code.Rmd"),
  #                to = reportmkd1,overwrite = TRUE)
        
        render(reportmkd)
        
  #      render(reportmkd1, quiet = TRUE)
        
        
        .Object
    }
)


#' @name SingleRepReport
#' @importFrom rtracklayer import
#' @importFrom rtracklayer import.bed
#' @title Final report for single group of regions
#' @description
#' When user call all steps in the pipeline, the final report can be generated.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' Any steps object in this package is acceptable when the pipeline is ready.
#' @param htmlOutput \code{Character} scalar.
#' HTML report file directory
#' Default: NULL ("Report.html")
#' @param ... Additional arguments, currently unused.
#' @details
#' The report is HTML format. All link in HTML file is the relative directory
#' in report step folder and other step folder
#' If user want to move HTML file and keep all link access available,
#' they should move the whole pipeline folder at the same time.
#' @return An invisible \code{\link{ATACProc-class}}
#' object (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacPipe}}



setGeneric("atacSingleRepReport",
           function(prevStep, htmlOutput = NULL,...)
               standardGeneric("atacSingleRepReport"))



#' @rdname SingleRepReport
#' @aliases atacSingleRepReport
#' @export
setMethod(
    f = "atacSingleRepReport",
    signature = "Step",
    definition = function(prevStep, htmlOutput = NULL, ...){
        allpara <- c(list(Class = "SingleRepReport",
                          prevSteps = list(prevStep), isReportStep = TRUE),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
