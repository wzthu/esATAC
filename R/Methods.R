atacPipe <- function(fastqInput1,fastqInput2=NULL, adapter1 = NULL, adapter2 = NULL,interleave = FALSE, saveTmp = TRUE){
    if(is.null(fastqInput2)&&!interleave&&is.null(adapter1)){
        stop("adapter1 should not be NULL for single end sequencing data")
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
        atacPeakQC(peakCalling,qcbedInput = "DHS",reportOutput = file.path(.obtainConfigure("tmpdir"),"DHSQC"))
        atacPeakQC(peakCalling,qcbedInput = "blacklist",reportOutput = file.path(.obtainConfigure("tmpdir"),"blacklistQC"))
        atacFripQC(atacProcReads = sam2Bed,atacProcPeak = peakCalling)
    }else{
        tssqc180_247 <-atacTSSQC(sam2Bed,reportPrefix = file.path(.obtainConfigure("tmpdir"),"tssqc180_247"),fregLenRange = c(180,247))
        fregLenDistr <- atacFregLenDistr(sam2Bed)
        shortBed <- atacBedUtils(sam2Bed,maxFregLen = 100, chrFilterList = NULL)
        peakCalling <- atacPeakCalling(shortBed)
        atacPeakQC(peakCalling,qcbedInput = "DHS",reportOutput = file.path(.obtainConfigure("tmpdir"),"DHSQC"))
        atacPeakQC(peakCalling,qcbedInput = "blacklist",reportOutput = file.path(.obtainConfigure("tmpdir"),"blacklistQC"))
        atacFripQC(atacProcReads = shortBed,atacProcPeak = peakCalling)
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


