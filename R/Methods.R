#' @name atacPipe
#' @title Unzip and merge fastq files
#' @description
#' Unzip and merge fastq files that are in format of bzip, gzip or fastq
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
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{atacBedUtils}}
#' @export
atacPipe <- function(fastqInput1,fastqInput2=NULL, adapter1 = NULL, adapter2 = NULL,interleave = FALSE, saveTmp = TRUE, createReport = TRUE){

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
    tssqc100 <-atacTSSQC(sam2Bed,reportPrefix = file.path(.obtainConfigure("tmpdir"),"tssqc100"),fregLenRange = c(0,100))

    if(is.null(fastqInput2)&&!interleave){
        peakCalling <- atacPeakCalling(sam2Bed)
        DHSQC <- atacPeakQC(peakCalling,qcbedInput = "DHS",reportOutput = file.path(.obtainConfigure("tmpdir"),"DHSQC"))
        blacklistQC <- atacPeakQC(peakCalling,qcbedInput = "blacklist",reportOutput = file.path(.obtainConfigure("tmpdir"),"blacklistQC"))
        fripQC <- atacFripQC(atacProcReads = sam2Bed,atacProcPeak = peakCalling)
    }else{
        tssqc180_247 <-atacTSSQC(sam2Bed,reportPrefix = file.path(.obtainConfigure("tmpdir"),"tssqc180_247"),fregLenRange = c(180,247))
        fregLenDistr <- atacFregLenDistr(sam2Bed)
        shortBed <- atacBedUtils(sam2Bed,maxFregLen = 100, chrFilterList = NULL)
        peakCalling <- atacPeakCalling(shortBed)
        DHSQC <- atacPeakQC(peakCalling,qcbedInput = "DHS",reportOutput = file.path(.obtainConfigure("tmpdir"),"DHSQC"))
        blacklistQC <- atacPeakQC(peakCalling,qcbedInput = "blacklist",reportOutput = file.path(.obtainConfigure("tmpdir"),"blacklistQC"))
        fripQC <- atacFripQC(atacProcReads = shortBed,atacProcPeak = peakCalling)
        Peakanno <- atacPeakAnno(atacProc = peakCalling, annoDb = "org.Hs.eg.db")
        goAna <- atacGOAnalysis(atacProc = Peakanno, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01)
        pwm <- readRDS(system.file("extdata", "motifPWM.rds", package="ATACFlow"))
        output_motifscan <- atacMotifScan(atacProc = peakCalling, motifPWM = pwm, min.score = "90%")
        cs_output <- atacExtractCutSite(atacProc = sam2Bed, prefix = "ATAC")
        footprint <- atacCutSiteCount(atacProcCutSite = cs_output, atacProcMotifScan = output_motifscan, strandLength = 100)
    }

    if(interleave){
        seqtype <- "paired end(interleave)"
        freg <- 2
    }else if(is.null(fastqInput2)){
        seqtype <- "single end"
        freg <- 1
    }else{
        seqtype <- "paired end"
        freg <- 2
    }
    if(is.null(fastqInput2)){
        filelist <- data.frame(`File(s)`=fastqInput1)
    }else{
        filelist <- data.frame(`Mate1 files`=fastqInput1,
                          `Mate2 files`=fastqInput2)
    }
    wholesummary = data.frame(Item=c("Sequence Type",
                                      "Original total reads",
                                      "Adapter removed reads for mapping",
                                      "Locations mapped once / twice / total",
                                      "Non-Redundant Fraction (NRF)",
                                      "PCR Bottlenecking Coefficients 1 (PBC1)",
                                      "PCR Bottlenecking Coefficients 2 (PBC2)",
                                      "Chromosome mitochondria reads ratio",
                                      "Total peaks",
                                      "Fraction of reads in peaks(FRiP)",
                                      "Peaks ratio overlaped with union DHS",
                                      "Peaks ratio overlaped with blacklist"),
                              Value=c(seqtype,
                                       removeAdapter$getReportVal("statisticslist")[[1]],
                                       as.integer(removeAdapter$getReportVal("statisticslist")[["Number of retained reads"]])/freg,
                                       sprintf("%s / %s / %s",libComplexQC$getReportVal("one"),libComplexQC$getReportVal("two"),libComplexQC$getReportVal("total")),
                                       sprintf("%.2f",as.numeric(libComplexQC$getReportVal("NRF"))),
                                       sprintf("%.2f",as.numeric(libComplexQC$getReportVal("PBC1"))),
                                       sprintf("%.2f",as.numeric(libComplexQC$getReportVal("PBC2"))),
                                       sprintf("%.2f",as.numeric(sam2Bed$getReportVal("filted"))/as.numeric(sam2Bed$getReportVal("total"))),
                                       sprintf("%d",as.numeric(fripQC$getReportVal("totalPeaks"))),
                                       sprintf("%.2f",as.numeric(fripQC$getReportVal("FRiP"))),
                                       sprintf("%.2f",as.numeric(DHSQC$getReportVal("qcbedRate"))),
                                       sprintf("%.2f",as.numeric(blacklistQC$getReportVal("qcbedRate"))))
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
                       atacProcs = atacProcs
                       )

    if(createReport){
        filename <- strsplit(fastqInput1,".fastq|.FASTQ|.FQ|.fq")[[1]][1]
        filename <- basename(filename)

        rmdfile<-system.file(package="ATACFlow", "extdata", "Report.Rmd")
        rmdtext<-readChar(rmdfile,nchars=file.info(rmdfile)$size,useBytes = TRUE)
        rmdtext<-sprintf(rmdtext,filename)

        workdir <- getwd()
        save(filelist,wholesummary,atacProcs,workdir,file = file.path(.obtainConfigure("tmpdir"),"Report.Rdata"))

        writeChar(rmdtext,con = file.path(.obtainConfigure("tmpdir"),"Report.Rmd"),useBytes = TRUE)
        render(file.path(.obtainConfigure("tmpdir"),"Report.Rmd"))
        #knit(file.path(.obtainConfigure("tmpdir"),"Report.Rmd"), file.path(.obtainConfigure("tmpdir"),"Report.md"))
        #markdownToHTML(file.path(.obtainConfigure("tmpdir"),"Report.md"), file.path(.obtainConfigure("tmpdir"),"Report.html"))
        #browseURL(paste0('file://', file.path(.obtainConfigure("tmpdir"),"Report.html")))
    }

    invisible(conclusion)



}




