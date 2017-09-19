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
    renamer <- atacRenamer(unzipAndMerge)
    removeAdapter <- atacRemoveAdapter(renamer, adapter1 = adapter1, adapter2 = adapter2)
    bowtie2Mapping <- atacBowtie2Mapping(removeAdapter)
    libComplexQC <- atacLibComplexQC(bowtie2Mapping)
    sam2Bed <-atacSam2Bed(bowtie2Mapping,maxFregLen = 2000)
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
                   fripQC = fripQC
    )
    conclusion <- list(filelist=filelist,
                       wholesummary = wholesummary,
                       atacProcs = atacProcs
                       )
                       
    if(createReport){
        filename <- strsplit(fastqInput1,".fastq|.FASTQ|.FQ|.fq")[[1]][1]
        filename <- basename(filename)
        
        rmdfile<-system.file(package="atacpipe", "extdata", "Report.Rmd")
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
   
    
    
    
}

atacPrintMap <-function(atacProc=NULL,preProc=FALSE,nextProc=TRUE,curProc=TRUE,display=TRUE){
    if(is.null(atacProc)){
        GraphMng$new()$printMap(display=display)
    }else if(class(atacProc)=="character"){
        GraphMng$new()$printMap(atacProc,preProc,nextProc,curProc,display=display)
    }else{
        atacProc$printMap(preProc,nextProc,curProc,display=display)
        invisible(atacProc)
    }
}

atacPrintNextList<-function(atacProc){
    if(class(atacProc)=="character"){
        GraphMng$new()$getNextProcs(atacProc)
    }else{
        GraphMng$new()$getNextProcs(atacProc$getProcName())
        invisible(atacProc)
    }
}

atacPrintPrevList<-function(atacProc){
    if(class(atacProc)=="character"){
        GraphMng$new()$getPrevProcs(atacProc)
    }else{
        GraphMng$new()$getPrevProcs(atacProc$getProcName())
        invisible(atacProc)
    }
}

atacClearCache <- function(atacProc){
    atacProc$clearCache()
    invisible(atacProc)
}

atacProcessing <- function(atacProc){
    atacProc$process()
    invisible(atacProc)
}


