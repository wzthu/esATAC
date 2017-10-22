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


getSuffix = function(filePath){
    filename<-basename(filePath)
    lst=strsplit(filename,"\\.")[[1]]
    if(length(lst)==1){
        return(NULL)
    }else{
        return(lst[length(lst)])
    }
}
getSuffixlessFileName = function(filePath){
    sfx=getSuffix(filePath)
    if(is.null(sfx)){
        return(basename(filePath))
    }else {
        return(strsplit(basename(filePath),paste0(".",sfx)))
    }
}



#' @name atacPipe
#' @title Pipeline for single replicate
#' @description
#' The pipeline to process sequencing data into destination files including
#' a HTML report file reads storage files (BED BAM)
#' and various quality control report files.
#' @param fastqInput1 \code{Character} vector. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates paired
#' with file paths in fastqInput2
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}
#' @param fastqInput2 \code{Character} vector. It contains file paths with #2
#' mates paired with file paths in fastqInput1
#' For single-end sequencing files and interleaved paired-end sequencing
#' files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' @param adapter1 \code{Character}. It is an adapter sequence for file1.
#' For single end data, it is requied.
#' @param adapter2 \code{Character}. It is an adapter sequence for file2.
#' @param interleave \code{Logical}. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param createReport \code{Logical}. If the HTML report file will be created.
#' @param prefix For identifying files.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{atacBedUtils}}
#' @examples
#' \dontrun{
#' td<-tempdir()
#' options(atacConf=setConfigure("threads",4))
#' dir.create(file.path(td,"ref"))
#' options(atacConf=setConfigure("refdir",file.path(td,"ref")))
#' options(atacConf=setConfigure("genome","hg19"))
#' bedbzfile11 <- c(
#'     system.file(package="ATACpipe", "extdata", "chr20_1.1.fq.gz"),
#'     system.file(package="ATACpipe", "extdata", "chr20_1.2.fq.bz2")
#' )
#' bedbzfile12 <- c(
#'     system.file(package="ATACpipe", "extdata", "chr20_2.1.fq.gz"),
#'     system.file(package="ATACpipe", "extdata", "chr20_2.2.fq.bz2")
#' )
#' # for single end
#' dir.create(file.path(td,"single"))
#' options(atacConf=setConfigure("tmpdir",file.path(td,"single")))
#' atacPipe(fastqInput1 = bedbzfile11,adapter1 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGACTGAAG")
#' # for paired end
#' dir.create(file.path(td,"paired"))
#' options(atacConf=setConfigure("tmpdir",file.path(td,"paired")))
#' atacPipe(fastqInput1 = bedbzfile11,fastqInput2 = bedbzfile12)
#' }
#' @export

atacPipe <- function(fastqInput1,fastqInput2=NULL, adapter1 = NULL, adapter2 = NULL,
                     interleave = FALSE,  createReport = TRUE, prefix = NULL){ #saveTmp = TRUE,


    if(is.null(fastqInput2)&&!interleave&&is.null(adapter1)){
        stop("adapter1 should not be NULL for single end sequencing data")
    }
    if(!is.null(fastqInput1)){
        fastqInput1 = normalizePath(fastqInput1)
    }
    if(!is.null(fastqInput2)){
        fastqInput2 = normalizePath(fastqInput2)
    }
    unzipAndMerge <- atacUnzipAndMerge(fastqInput1 = fastqInput1,fastqInput2 = fastqInput2,interleave = interleave)
    atacQC <- atacQCReport(atacProc = unzipAndMerge)
    renamer <- atacRenamer(unzipAndMerge)
    removeAdapter <- atacRemoveAdapter(renamer, adapter1 = adapter1, adapter2 = adapter2)
    bowtie2Mapping <- atacBowtie2Mapping(removeAdapter)
    libComplexQC <- atacLibComplexQC(bowtie2Mapping)
    sam2Bed <-atacSamToBed(bowtie2Mapping,maxFregLen = 2000)
    bedToBigWig <- atacBedToBigWig(sam2Bed)
    tssqc100 <-atacTSSQC(sam2Bed,reportPrefix = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName(fastqInput1[1]),".tssqc100")),fregLenRange = c(0,100))

    if(is.null(fastqInput2)&&!interleave){
        peakCalling <- atacPeakCalling(sam2Bed)
        DHSQC <- atacPeakQC(peakCalling,qcbedInput = "DHS",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName(fastqInput1[1]),".DHSQC")))
        blacklistQC <- atacPeakQC(peakCalling,qcbedInput = "blacklist",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName(fastqInput1[1]),".blacklistQC")))
        fripQC <- atacFripQC(atacProcReads = sam2Bed,atacProcPeak = peakCalling)
    }else{
        tssqc180_247 <-atacTSSQC(sam2Bed,reportPrefix = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName(fastqInput1[1]),".tssqc180_247")),fregLenRange = c(180,247))
        fregLenDistr <- atacFregLenDistr(sam2Bed)
        shortBed <- atacBedUtils(sam2Bed,maxFregLen = 100, chrFilterList = NULL)
        peakCalling <- atacPeakCalling(shortBed)
        DHSQC <- atacPeakQC(peakCalling,qcbedInput = "DHS",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName(fastqInput1[1]),".DHSQC")))
        blacklistQC <- atacPeakQC(peakCalling,qcbedInput = "blacklist",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName(fastqInput1[1]),".blacklistQC")))
        fripQC <- atacFripQC(atacProcReads = shortBed,atacProcPeak = peakCalling)

        Peakanno <- atacPeakAnno(atacProc = peakCalling)
        goAna <- atacGOAnalysis(atacProc = Peakanno, ont = "BP", pvalueCutoff = 0.01)
        pwm <- readRDS(system.file("extdata", "motifPWM.rds", package="ATACpipe"))
        output_motifscan <- atacMotifScan(atacProc = peakCalling, motifPWM = pwm, min.score = "85%", prefix = prefix)
        cs_output <- atacExtractCutSite(atacProc = sam2Bed, prefix = prefix)
        footprint <- atacCutSiteCount(atacProcCutSite = cs_output, atacProcMotifScan = output_motifscan,
                                      strandLength = 100, prefix = prefix)
    }

    if(interleave){
        seqtype <- "paired end (PE,interleave)"
        freg <- 2
    }else if(is.null(fastqInput2)){
        seqtype <- "single end (SE)"
        freg <- 1
    }else{
        seqtype <- "paired end (PE)"
        freg <- 2
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
                                          getVMRShow(as.integer(getReportVal(removeAdapter,"statisticslist")[["Number of retained reads"]])/freg,
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
                    getVMRShow(as.integer(getReportVal(removeAdapter,"statisticslist")[["Number of retained reads"]])/freg,
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
                                      getVMRShow(as.integer(getReportVal(removeAdapter,"statisticslist")[["Number of retained reads"]])/freg,
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
                 getVMRShow(as.integer(getReportVal(removeAdapter,"statisticslist")[["Number of retained reads"]])/freg,
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
                   fregLenDistr = fregLenDistr,
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

        rmdfile<-system.file(package="ATACpipe", "extdata", "Report.Rmd")
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
#' @title Pipeline for single replicate
#' @description
#' The pipeline to process case control study sequencing data
#' into destination files including
#' a HTML report file reads storage files (BED BAM)
#' and various quality control report files.
#' @param case \code{List} scalar. Input for case sample. \code{fastqInput1},
#' the path of the mate 1 fastq file, is required. \code{fastqInput2},
#' the path of the mate 2 fastq file, is required, when \code{interleave=FALSE}.
#' \code{adapter1} and \code{adapter2} are optional.
#' @param control \code{List} scalar. Input for control sample. \code{fastqInput1},
#' the path of the mate 1 fastq file, is required. \code{fastqInput2},
#' the path of the mate 2 fastq file, is required, when \code{interleave=FALSE}.
#' \code{adapter1} and \code{adapter2} are optional.
#' @param interleave \code{Logical}. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param createReport \code{Logical}. If the HTML report file will be created.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacPipe}}
#' @export
#'
atacPipe2 <- function(case = list(fastqInput1="paths/To/fastq1",fastqInput2="paths/To/fastq2", adapter1 = NULL, adapter2 = NULL),
                      control =list(fastqInput1="paths/To/fastq1",fastqInput2="paths/To/fastq2", adapter1 = NULL, adapter2 = NULL),
                      interleave = FALSE, createReport = TRUE){ #saveTmp = TRUE,
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
    caselist <- atacPipe(fastqInput1 = case[["fastqInput1"]],fastqInput2 = case[["fastqInput2"]],
               adapter1 = case[["adapter1"]], adapter2 = case[["adapter2"]],interleave = interleave,
                createReport = FALSE, prefix = "CASE_all_data") #saveTmp = TRUE,
    ctrllist <- atacPipe(fastqInput1 = control[["fastqInput1"]],fastqInput2 = control[["fastqInput2"]],
               adapter1 = control[["adapter1"]], adapter2 = control[["adapter2"]],interleave = interleave,
                createReport = FALSE, prefix = "CTRL_all_data") #saveTmp = TRUE,

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

        rmdfile<-system.file(package="ATACpipe", "extdata", "Report2.Rmd")
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

