getPer<- function(per){
    per <- as.numeric(per)*100
    return(paste0(sprintf("%.1f",per),"%"))
}
getf<- function(per){
    per <- as.numeric(per)
    return(sprintf("%.2f",per))
}
getVM <- function(readSize){
    readSize <- as.numeric(readSize)
    return(c(readSize,readSize/1e6))
}
getVMR <- function(readSize,total){
    readSize <- as.character(readSize)
    total <- as.integer(total)
    vm <- getVM(readSize)
    return(c(vm, 100*vm[1]/total,total))
}
getVMShow <- function(readSize,detail = FALSE){
    vm <- getVM(readSize)
    if(detail){
        return(sprintf("%.1fM (%d)",vm[2],vm[1]))
    }else{
        return(sprintf("%.1fM",vm[2]))
    }

}
getVMRShow <- function(readSize,total,detail = FALSE){
    vmr <- getVMR(readSize,total)
    total <- as.numeric(total)
    if(detail){
        return(sprintf("%.1fM (%.2f%s, %d / %d)",vmr[2],vmr[3],"%",vmr[1],vmr[4]))
    }else{
        return(sprintf("%.1fM (%.2f%s)",vmr[2],vmr[3],"%"))
    }

}
getRshow <- function(readSize,total,detail = FALSE){
    vmr <- getVMR(readSize,total)
    total <- as.numeric(total)
    if(detail){
        return(sprintf("%.2f, %d / %d",vmr[3]/100,vmr[1],vmr[4]))
    }else{
        return(sprintf("%.2f",vmr[3]/100))
    }

}


getSuffix0 <- function(filePath){
    filename<-basename(filePath)
    lst=strsplit(filename,"\\.")[[1]]
    if(length(lst)==1){
        return(NULL)
    }else{
        return(lst[length(lst)])
    }
}
getSuffixlessFileName0 <- function(filePath){
    sfx=getSuffix0(filePath)
    if(is.null(sfx)){
        return(basename(filePath))
    }else {
        return(strsplit(basename(filePath),paste0(".",sfx)))
    }
}


#' @docType package
#' @name esATAC-package
#' @details
#' See packageDescription('esATAC') for package details.
#'
#' @title An Easy-to-use Systematic pipeline for ATACseq data analysis
#' @description
#' This package provides a framework and complete preset pipeline for
#' the quantification and analysis of ATAC-seq Reads. It covers raw sequencing
#' reads preprocessing (FASTQ files), reads alignment (Rbowtie2), aligned reads
#' file operation (SAM, BAM, and BED files), peak calling (fseq), genome
#' annotations (Motif, GO, SNP analysis) and quality control report. The package
#' is managed by dataflow graph. It is easy for user to pass variables seamlessly
#' between processes and understand the workflow. Users can process FASTQ files
#' through end-to-end preset pipeline which produces a pretty HTML report for
#' quality control and preliminary statistical results, or customize workflow
#' starting from any intermediate stages with esATAC functions easily and flexibly.
#'
#' Preset pipeline for case study is shown below.
#' For case-control study, see \code{\link{atacPipe2}}.
#'
#'
#' NOTE:
#' Build bowtie index in the function may take some time.
#' If you already have bowtie2 index files or
#' you want to download(\url{ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes})
#' instead of building,
#' you can let esATAC skip the steps by renaming them following the format
#' (genome+suffix) and put them in reference installation path (refdir).
#' Example: hg19 bowtie2 index files
#'
#' \itemize{
#' \item hg19.1.bt2
#' \item hg19.2.bt2
#' \item hg19.3.bt2
#' \item hg19.4.bt2
#' \item hg19.rev.1.bt2
#' \item hg19.rev.2.bt2
#' }
#'
#' For single end reads FASTQ files,
#' The required parameters are fastqInput1 and adapter1.
#' For paired end reads non-interleaved FASTQ files (interleave=FALSE,defualt),
#' The required parameters are fastqInput1 and fastqInput2.
#' Otherwise, parameter fastqInput2 is not required (interleave=TRUE)
#'
#' The paths of sequencing data replicates can be a \code{Character} vector.
#' For example:
#'
#' fastqInput1=c("file_1.rep1.fastq","file_1.rep2.fastq")
#'
#' fastqInput2=c("file_2.rep1.fastq","file_2.rep2.fastq")
#'
#' The result will be return by the function.
#' An HTML report file will be created for paired end reads.
#' Intermediate files will be save at tmpdir path (default is ./)
#'
#' @param fastqInput1 \code{Character} vector. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates paired
#' with file paths in fastqInput2
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}
#' @param fastqInput2 \code{Character} vector. It contains file paths with #2
#' mates paired with file paths in fastqInput1.
#' For single-end sequencing files and interleaved paired-end sequencing
#' files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' @param adapter1 \code{Character} scalar. It is an adapter sequence for file1.
#' For single end data, it is requied.
#' @param adapter2 \code{Character} scalar. It is an adapter sequence for file2.
#' @param interleave \code{Logical} scalar. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param createReport \code{Logical} scalar. If the HTML report file will be created.
#' @param motifPWM \code{List} scalar. Motif PWM, a list.
#' @param prefix \code{Character} scalar. Temporary file prefix for identifying files
#' when multiple pipeline generating file in the same tempdir.
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param ... Configure "refdir", "genome", "threads", "tmpdir" for this function.
#' They will overwrite global configuration.
#' If you need to set globally, see \link{configureValue}.
#' \describe{
#'   \item{refdir}{\code{Character} scalar, the path for reference data being installed to and storage.}
#'   \item{genome}{\code{Character} scalar, the genome(like hg19, mm10, etc.) reference data in "refdir" to be used in the pipeline.}
#'   \item{tmpdir}{\code{Character} scalar, the temporary file storage path}
#'   \item{threads}{\code{Integer} scalar, the max threads allowed to be created}
#' }
#' Parameter "chr", which chromatin the program will processing.
#' \describe{
#'    \item{chr}{\code{Character} scalar, It must be identical with the filename of cut site information files or a subset. Default:c(1:22, "X", "Y").}
#' }
#' @return \code{List} scalar. It is a list that save the result of the pipeline.
#' Slot "filelist": the input file paths.
#' Slot "wholesummary": a dataframe that for quality control summary
#' Slot "atacProcs": \code{\link{ATACProc-class}} objects generated by each process in the pipeline.
#' Slot "filtstat": a dataframe that summary the reads filted in each process.
#'
#' @author Zheng Wei and Wei Zhang
#' @seealso
#' \code{\link{printMap}},
#' \code{\link{atacPipe2}},
#' \code{\link{atacRenamer}},
#' \code{\link{atacRemoveAdapter}},
#' \code{\link{atacBowtie2Mapping}},
#' \code{\link{atacPeakCalling}},
#' \code{\link{atacMotifScan}}
#' @examples
#' \dontrun{
#' ## These codes are time consuming so they will not be run and
#' ## checked by bioconductor checker.
#'
#'
#' # call pipeline
#' # for a quick example(only CTCF will be processing)
#' conclusion <-
#'   atacPipe(
#'        # MODIFY: Change these paths to your own case files!
#'        # e.g. fastqInput1 = "your/own/data/path.fastq"
#'        fastqInput1 = system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz"),
#'        fastqInput2 = system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz"),
#'        # MODIFY: Set the genome for your data
#'        genome = "hg19",
#'        motifPWM = getMotifPWM(motif.file = system.file("extdata", "CTCF.txt", package="esATAC"), is.PWM = FALSE))
#'
#' # call pipeline
#' # for overall example(all human motif in JASPAR will be processed)
#' conclusion <-
#'   atacPipe(
#'        # MODIFY: Change these paths to your own case files!
#'        # e.g. fastqInput1 = "your/own/data/path.fastq"
#'        fastqInput1 = system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz"),
#'        fastqInput2 = system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz"),
#'        # MODIFY: Set the genome for your data
#'        genome = "hg19")
#' }
#' @export

atacPipe <- function(fastqInput1,fastqInput2=NULL, adapter1 = NULL, adapter2 = NULL,
                     interleave = FALSE,  createReport = TRUE, motifPWM = NULL, prefix = NULL,
                     chr = c(1:22, "X", "Y"), ...){ #saveTmp = TRUE,



    if(is.null(fastqInput2)&&!interleave&&is.null(adapter1)){
        stop("adapter1 should not be NULL for single end sequencing data")
    }
    if(!is.null(fastqInput1)){
        fastqInput1 = normalizePath(fastqInput1)
    }
    if(!is.null(fastqInput2)){
        fastqInput2 = normalizePath(fastqInput2)
    }

    param.tmp <- list(...)
    if(!(!is.null(param.tmp[["dontSet"]])&&param.tmp[["dontSet"]])){
        if(!is.null(param.tmp[["refdir"]])){
            options(atacConf=setConfigure("refdir",param.tmp[["refdir"]]))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","refdir"))){
                dir.create(file.path("esATAC_pipeline","refdir"))
            }
            options(atacConf=setConfigure("refdir",file.path("esATAC_pipeline","refdir")))
        }
        if(!is.null(param.tmp[["threads"]])){
            options(atacConf=setConfigure("threads",as.integer(param.tmp[["threads"]])))
            message(getConfigure("threads"))
        }else{
            options(atacConf=setConfigure("threads",2L))
        }
        if(!is.null(param.tmp[["tmpdir"]])){
            options(atacConf=setConfigure("tmpdir",param.tmp[["tmpdir"]]))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","esATAC_result"))){
                dir.create(file.path("esATAC_pipeline","esATAC_result"))
            }
            options(atacConf=setConfigure("tmpdir",file.path("esATAC_pipeline","esATAC_result")))
        }
        if(!is.null(param.tmp[["genome"]])){
            options(atacConf=setConfigure("genome",param.tmp[["genome"]]))
        }else{
            stop("parameter genome is required")
        }
    }
    


    unzipAndMerge <- atacUnzipAndMerge(fastqInput1 = fastqInput1,fastqInput2 = fastqInput2,interleave = interleave)
    atacQC <- atacQCReport(atacProc = unzipAndMerge)
    renamer <- atacRenamer(unzipAndMerge)
    removeAdapter <- atacRemoveAdapter(renamer, adapter1 = adapter1, adapter2 = adapter2)
    bowtie2Mapping <- atacBowtie2Mapping(removeAdapter)
    libComplexQC <- atacLibComplexQC(bowtie2Mapping)
    sam2Bed <-atacSamToBed(bowtie2Mapping,maxFragLen = 2000)
    bedToBigWig <- atacBedToBigWig(sam2Bed)
    tssqc100 <-atacTSSQC(sam2Bed,reportPrefix = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(fastqInput1[1]),".tssqc100")),fragLenRange = c(0,100))

    if(is.null(fastqInput2)&&!interleave){
        peakCalling <- atacPeakCalling(sam2Bed)
        DHSQC <- atacPeakQC(peakCalling,qcbedInput = "DHS",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(fastqInput1[1]),".DHSQC")))
        blacklistQC <- atacPeakQC(peakCalling,qcbedInput = "blacklist",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(fastqInput1[1]),".blacklistQC")))
        fripQC <- atacFripQC(atacProcReads = sam2Bed,atacProcPeak = peakCalling)
    }else{
        tssqc180_247 <-atacTSSQC(sam2Bed,reportPrefix = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(fastqInput1[1]),".tssqc180_247")),fragLenRange = c(180,247))
        fragLenDistr <- atacFragLenDistr(sam2Bed)
        shortBed <- atacBedUtils(sam2Bed,maxFragLen = 100, chrFilterList = NULL)
        peakCalling <- atacPeakCalling(shortBed)
        DHSQC <- atacPeakQC(peakCalling,qcbedInput = "DHS",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(fastqInput1[1]),".DHSQC")))
        blacklistQC <- atacPeakQC(peakCalling,qcbedInput = "blacklist",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(fastqInput1[1]),".blacklistQC")))
        fripQC <- atacFripQC(atacProcReads = shortBed,atacProcPeak = peakCalling)



        if(is.null(motifPWM)){
            pwm <- getMotifPWM(JASPARdb = TRUE, Species = "9606")
        }else{
            pwm <- motifPWM
        }

        Peakanno <- atacPeakAnno(atacProc = peakCalling)
        goAna <- atacGOAnalysis(atacProc = Peakanno, ont = "BP", pvalueCutoff = 0.01)
        output_motifscan <- atacMotifScan(atacProc = peakCalling, motifPWM = pwm, min.score = "85%", prefix = prefix)
        cs_output <- atacExtractCutSite(atacProc = sam2Bed, prefix = prefix)
        footprint <- atacCutSiteCount(atacProcCutSite = cs_output, atacProcMotifScan = output_motifscan,
                                      strandLength = 100, prefix = prefix, chr = chr)
    }

    if(interleave){
        seqtype <- "paired end (PE,interleave)"
        frag <- 2
    }else if(is.null(fastqInput2)){
        seqtype <- "single end (SE)"
        frag <- 1
    }else{
        seqtype <- "paired end (PE)"
        frag <- 2
    }
    if(is.null(fastqInput2)){
        filelist <- data.frame(`File(s)`=fastqInput1)
    }else{
        filelist <- data.frame(`Mate1 files`=fastqInput1,
                          `Mate2 files`=fastqInput2)
    }
    if(is.null(fastqInput2)&&!interleave){
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
                                          getPer(getReportVal(DHSQC,"qcbedRate")),
                                          getPer(getReportVal(blacklistQC,"qcbedRate")),
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
        atacProcs=list(unzipAndMerge = unzipAndMerge,
                       renamer = renamer,
                       removeAdapter = removeAdapter,
                       bowtie2Mapping = bowtie2Mapping,
                       libComplexQC = libComplexQC,
                       sam2Bed = sam2Bed,
                       bedToBigWig = bedToBigWig,
                       tssqc100 = tssqc100,
                       peakCalling = peakCalling,
                       DHSQC = DHSQC,
                       blacklistQC = blacklistQC,
                       fripQC = fripQC,
                       atacQC = atacQC
        )
        conclusion <- list(filelist=filelist,
                           wholesummary = wholesummary,
                           atacProcs = atacProcs,
                           filtstat = filtstat
        )
        return(conclusion)
    }

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
                                      getPer(getReportVal(DHSQC,"qcbedRate")),
                                      getPer(getReportVal(blacklistQC,"qcbedRate")),
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
    atacProcs=list(unzipAndMerge = unzipAndMerge,
                   renamer = renamer,
                   removeAdapter = removeAdapter,
                   bowtie2Mapping = bowtie2Mapping,
                   libComplexQC = libComplexQC,
                   sam2Bed = sam2Bed,
                   bedToBigWig = bedToBigWig,
                   tssqc100 = tssqc100,
                   tssqc180_247 = tssqc180_247,
                   fragLenDistr = fragLenDistr,
                   peakCalling = peakCalling,
                   DHSQC = DHSQC,
                   blacklistQC = blacklistQC,
                   fripQC = fripQC,
                   shortBed = shortBed,
                   Peakanno = Peakanno,
                   goAna = goAna,
                   output_motifscan = output_motifscan,
                   cs_output = cs_output,
                   footprint = footprint,
                   atacQC = atacQC
    )
    conclusion <- list(filelist=filelist,
                       wholesummary = wholesummary,
                       atacProcs = atacProcs,
                       filtstat = filtstat
                       )

    if(createReport){
        filename <- strsplit(fastqInput1,".fastq|.FASTQ|.FQ|.fq")[[1]][1]
        filename <- basename(filename)

        rmdfile<-system.file(package="esATAC", "extdata", "Report.Rmd")
        rmdtext<-readChar(rmdfile,nchars=file.info(rmdfile)$size,useBytes = TRUE)
        #rmdtext<-sprintf(rmdtext,filename)

        workdir <- getwd()
        save(filelist,wholesummary,filtstat,atacProcs,workdir,file = file.path(.obtainConfigure("tmpdir"),"Report.Rdata"))

        writeChar(rmdtext,con = file.path(.obtainConfigure("tmpdir"),"Report.Rmd"),useBytes = TRUE)
        render(file.path(.obtainConfigure("tmpdir"),"Report.Rmd"))
        #knit(file.path(.obtainConfigure("tmpdir"),"Report.Rmd"), file.path(.obtainConfigure("tmpdir"),"Report.md"))
        #markdownToHTML(file.path(.obtainConfigure("tmpdir"),"Report.md"), file.path(.obtainConfigure("tmpdir"),"Report.html"))
        #browseURL(paste0('file://', file.path(.obtainConfigure("tmpdir"),"Report.html")))
    }

    invisible(conclusion)



}



#' @name atacPipe2
#' @title Pipeline for case-control paired-end sequencing data
#' @description
#' The preset pipeline to process case control study sequencing data.
#' An HTML report file, result files(e.g. BED, BAM files) and
#' conclusion list will generated. See detail for usage.
#' @param case \code{List} scalar. Input for case sample. \code{fastqInput1},
#' the path(s) of the mate 1 fastq file(s), is required. \code{fastqInput2},
#' the path(s) of the mate 2 fastq file(s), is required, when \code{interleave=FALSE}.
#' \code{adapter1} and \code{adapter2} are optional.
#' @param control \code{List} scalar. Input for control sample. \code{fastqInput1},
#' the path(s) of the mate 1 fastq file(s), is required. \code{fastqInput2},
#' the path(s) of the mate 2 fastq file(s), is required, when \code{interleave=FALSE}.
#' \code{adapter1} and \code{adapter2} are optional.
#' @param interleave \code{Logical} scalar. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param createReport \code{Logical} scalar. If the HTML report file will be created.
#' @param motifPWM \code{List} scalar. Motif PWM list.
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param ... Configure "refdir", "genome", "threads", "tmpdir" for this function.
#' They will overwrite global configuration.
#' If you need to set globally, see \link{configureValue}.
#' \describe{
#'   \item{refdir}{\code{Character} scalar, the path for reference data being installed to and storage.}
#'   \item{genome}{\code{Character} scalar, the genome(like hg19, mm10, etc.) reference data in "refdir" to be used in the pipeline.}
#'   \item{tmpdir}{\code{Character} scalar, the temporary file storage path}
#'   \item{threads}{\code{Integer} scalar, the max threads allowed to be created}
#' }
#' Parameter "chr", which chromatin the program will processing.
#' \describe{
#'    \item{chr}{\code{Character} scalar, It must be identical with the filename of cut site information files or a subset. Default:c(1:22, "X", "Y").}
#' }
#' @return \code{List} scalar. It is a list that save the result of the pipeline.
#' Slot "wholesummary": a dataframe for quality control summary of  case and control data
#' Slot "caselist" and "ctrlist": Each of them is a list that save the result for case or control data.
#' Slots of "caselist" and "ctrllist":
#' Slot "filelist": the input file paths.
#' Slot "wholesummary": a dataframe for quality control summary of case or control data
#' Slot "atacProcs": \code{\link{ATACProc-class}} objects generated by each process in the pipeline.
#' Slot "filtstat": a dataframe that summary the reads filted in each process.
#' @details
#' NOTE:
#' Build bowtie index in this function may take some time. If you already have bowtie2 index files or you want to download(ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes) instead of building, you can let esATAC skip the steps by renaming them following the format (genome+suffix) and put them in reference installation path (refdir).
#' Example: hg19 bowtie2 index files
#'
#' \itemize{
#' \item hg19.1.bt2
#' \item hg19.2.bt2
#' \item hg19.3.bt2
#' \item hg19.4.bt2
#' \item hg19.rev.1.bt2
#' \item hg19.rev.2.bt2
#' }
#'
#' For single end reads FASTQ files,
#' The required parameters are fastqInput1 and adapter1.
#' For paired end reads non-interleaved FASTQ files (interleave=FALSE,defualt),
#' The required parameters are fastqInput1 and fastqInput2.
#' Otherwise, parameter fastqInput2 is not required (interleave=TRUE)
#'
#' The paths of sequencing data replicates can be a \code{Character} vector.
#' For example:
#'
#' fastqInput1=c("file_1.rep1.fastq","file_1.rep2.fastq")
#'
#' fastqInput2=c("file_2.rep1.fastq","file_2.rep2.fastq")
#'
#' The result will be return by the function.
#' An HTML report file will be created for paired end reads.
#' Intermediate files will be save at tmpdir path (default is ./)
#' @author Zheng Wei and Wei Zhang
#' @seealso
#' \code{\link{atacPipe}}
#' @import JASPAR2016
#' @importFrom TFBSTools getMatrixSet
#' @importFrom TFBSTools toPWM
#' @importFrom TFBSTools name
#' @examples
#' \dontrun{
#' ## These codes are time consuming so they will not be run and
#' ## checked by bioconductor checker.
#'
#'
#' # call pipeline
#' # for a quick example(only CTCF will be processed)
#' conclusion <-
#'    atacPipe2(
#'        # MODIFY: Change these paths to your own case files!
#'        # e.g. fastqInput1 = "your/own/data/path.fastq"
#'        case=list(fastqInput1 = system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz"),
#'                 fastqInput2 = system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz")),
#'        # MODIFY: Change these paths to your own control files!
#'        # e.g. fastqInput1 = "your/own/data/path.fastq"
#'        control=list(fastqInput1 = system.file(package="esATAC", "extdata", "chr20_1.2.fq.bz2"),
#'                     fastqInput2 = system.file(package="esATAC", "extdata", "chr20_2.2.fq.bz2")),
#'        # MODIFY: Set the genome for your data
#'        genome = "hg19",
#'        motifPWM = getMotifPWM(motif.file = system.file("extdata", "CTCF.txt", package="esATAC"), is.PWM = FALSE))
#'
#' # call pipeline
#' # for overall example(all human motif in JASPAR will be processed)
#' conclusion <-
#'    atacPipe2(
#'        # MODIFY: Change these paths to your own case files!
#'        # e.g. fastqInput1 = "your/own/data/path.fastq"
#'        case=list(fastqInput1 = system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz"),
#'                 fastqInput2 = system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz")),
#'        # MODIFY: Change these paths to your own control files!
#'        # e.g. fastqInput1 = "your/own/data/path.fastq"
#'        control=list(fastqInput1 = system.file(package="esATAC", "extdata", "chr20_1.2.fq.bz2"),
#'                     fastqInput2 = system.file(package="esATAC", "extdata", "chr20_2.2.fq.bz2")),
#'        # MODIFY: Set the genome for your data
#'        genome = "hg19")
#'}
#' @export
#'
atacPipe2 <- function(case = list(fastqInput1="paths/To/fastq1",fastqInput2="paths/To/fastq2", adapter1 = NULL, adapter2 = NULL),
                      control =list(fastqInput1="paths/To/fastq1",fastqInput2="paths/To/fastq2", adapter1 = NULL, adapter2 = NULL),
                      interleave = FALSE, createReport = TRUE, motifPWM = NULL, chr = c(1:22, "X", "Y"), ...){ #saveTmp = TRUE,
    if(case[["fastqInput1"]]=="paths/To/fastq1"||is.null(case[["fastqInput1"]])){
        stop("fastqInput1 for case can not be NULL")
    }
    if(!interleave && (case[["fastqInput2"]]=="paths/To/fastq2"|| is.null(case[["fastqInput2"]]))){
        stop("fastqInput2 for case can not be NULL")
    }
    if(control[["fastqInput1"]]=="paths/To/fastq1"||is.null(control[["fastqInput1"]])){
        stop("fastqInput1 for control can not be NULL")
    }
    if(!interleave && (control[["fastqInput2"]]=="paths/To/fastq2"||is.null(control[["fastqInput2"]]))){
        stop("fastqInput2 for control can not be NULL")
    }

    param.tmp <- list(...)
    if(!(!is.null(param.tmp[["dontSet"]])&&param.tmp[["dontSet"]])){
        if(!is.null(param.tmp[["refdir"]])){
            options(atacConf=setConfigure("refdir",param.tmp[["refdir"]]))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","refdir"))){
                dir.create(file.path("esATAC_pipeline","refdir"))
            }
            options(atacConf=setConfigure("refdir",file.path("esATAC_pipeline","refdir")))
        }
        if(!is.null(param.tmp[["threads"]])){
            options(atacConf=setConfigure("threads",as.integer(param.tmp[["threads"]])))
            message(getConfigure("threads"))
        }else{
            options(atacConf=setConfigure("threads",2L))
        }
        if(!is.null(param.tmp[["tmpdir"]])){
            options(atacConf=setConfigure("tmpdir",param.tmp[["tmpdir"]]))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","esATAC_result"))){
                dir.create(file.path("esATAC_pipeline","esATAC_result"))
            }
            options(atacConf=setConfigure("tmpdir",file.path("esATAC_pipeline","esATAC_result")))
        }
        if(!is.null(param.tmp[["genome"]])){
            options(atacConf=setConfigure("genome",param.tmp[["genome"]]))
        }else{
            stop("parameter genome is required")
        }
    }

    if(is.null(motifPWM)){
        pwm <- getMotifPWM(JASPARdb = TRUE, Species = "9606")
    }else{
        pwm <- motifPWM
    }

    caselist <- atacPipe(fastqInput1 = case[["fastqInput1"]],fastqInput2 = case[["fastqInput2"]],
               adapter1 = case[["adapter1"]], adapter2 = case[["adapter2"]],interleave = interleave,
                createReport = FALSE, motifPWM =pwm, prefix = "CASE_all_data", chr = chr, dontSet=TRUE) #saveTmp = TRUE,
    ctrllist <- atacPipe(fastqInput1 = control[["fastqInput1"]],fastqInput2 = control[["fastqInput2"]],
               adapter1 = control[["adapter1"]], adapter2 = control[["adapter2"]],interleave = interleave,
                createReport = FALSE, motifPWM =pwm, prefix = "CTRL_all_data", chr = chr, dontSet=TRUE) #saveTmp = TRUE,


    bed.case <- getParam(caselist$atacProcs$sam2Bed, "bedOutput")
    bed.ctrl <- getParam(ctrllist$atacProcs$sam2Bed, "bedOutput")

    case.peak <- getParam(caselist$atacProcs$peakCalling, "bedOutput")
    ctrl.peak <- getParam(ctrllist$atacProcs$peakCalling,"bedOutput")

    peakCom <- peakcomp(bedInput1 = case.peak, bedInput2 = ctrl.peak)
    case_specific.peak <- getParam(peakCom, "bedOutput")[1]
    ctrl_specific.peak <- getParam(peakCom, "bedOutput")[2]
    overlap.peak <- getParam(peakCom, "bedOutput")[3]

    # for case
    Peakanno.case <- peakanno(peakInput = case_specific.peak)
    goAna.case <- atacGOAnalysis(atacProc = Peakanno.case, ont = "BP", pvalueCutoff = 0.01)

    # for ctrl
    Peakanno.ctrl <- peakanno(peakInput = ctrl_specific.peak)
    goAna.ctrl <- atacGOAnalysis(atacProc = Peakanno.ctrl, ont = "BP", pvalueCutoff = 0.01)

    mout <- atacMotifScanPair(atacProc = peakCom, motifPWM = pwm, min.score = "90%")
    cs_case <- extractcutsite(bedInput = bed.case, prefix = "CASE")
    cs_ctrl <- extractcutsite(bedInput = bed.ctrl, prefix = "CTRL")

    footprint.case <- atacCutSiteCount(atacProcCutSite = cs_case,
                                       motif_info = getParam(mout, "rdsOutput.peak1"),
                                       strandLength = 100, prefix = "Case", chr = chr)

    footprint.ctrl <- atacCutSiteCount(atacProcCutSite = cs_ctrl,
                                       motif_info = getParam(mout, "rdsOutput.peak2"),
                                       strandLength = 100, prefix = "Ctrl", chr = chr)

    comp_result <- list(
        peakCom = peakCom,
        goAna.case = goAna.case,
        goAna.ctrl = goAna.ctrl,
        mout = mout,
        footprint.case = footprint.case,
        footprint.ctrl = footprint.ctrl
    )

    wholesummary <- data.frame(Item = caselist[["wholesummary"]][["Item"]],
                          Case = caselist[["wholesummary"]][["Value"]],
                          Control = ctrllist[["wholesummary"]][["Value"]],
                          Reference = ctrllist[["wholesummary"]][["Reference"]])

    conclusion <- list(caselist = caselist,
                       ctrllist = ctrllist,
                       wholesummary = wholesummary
                )
    casefilelist <- caselist[["filelist"]]
    ctrlfilelist <- ctrllist[["filelist"]]

    filtstat <- data.frame(Item = caselist[["filtstat"]][["Item"]],
                           Case = caselist[["filtstat"]][["Value"]],
                           Control = ctrllist[["filtstat"]][["Value"]],
                           Reference = ctrllist[["filtstat"]][["Reference"]])

    if(createReport){
        #filename <- strsplit(case[["fastqInput1"]],".fastq|.FASTQ|.FQ|.fq")[[1]][1]
        #filename <- basename(filename)

        rmdfile<-system.file(package="esATAC", "extdata", "Report2.Rmd")
        rmdtext<-readChar(rmdfile,nchars=file.info(rmdfile)$size,useBytes = TRUE)
        #rmdtext<-sprintf(rmdtext,filename)

        workdir <- getwd()
        save(casefilelist,ctrlfilelist,wholesummary,filtstat,caselist,ctrllist,workdir,file = file.path(.obtainConfigure("tmpdir"),"Report2.Rdata"))

        writeChar(rmdtext,con = file.path(.obtainConfigure("tmpdir"),"Report2.Rmd"),useBytes = TRUE)
        render(file.path(.obtainConfigure("tmpdir"),"Report2.Rmd"))
        #knit(file.path(.obtainConfigure("tmpdir"),"Report.Rmd"), file.path(.obtainConfigure("tmpdir"),"Report.md"))
        #markdownToHTML(file.path(.obtainConfigure("tmpdir"),"Report.md"), file.path(.obtainConfigure("tmpdir"),"Report.html"))
        #browseURL(paste0('file://', file.path(.obtainConfigure("tmpdir"),"Report.html")))
    }

    invisible(conclusion)

}

