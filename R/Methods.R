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
#' Preset pipeline for single replicate case study is shown below.
#'
#' For multi-replicates case study, see \code{\link{atacRepsPipe}}.
#'
#' For single replicate case-control study, see \code{\link{atacPipe2}}.
#'
#' For multi-replicates case-control study, see \code{\link{atacRepsPipe2}}.
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
#' @param genome \code{Character} scalar. The genome(like hg19, mm10, etc.) reference data in "refdir" to be used in the pipeline.
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
#' @param refdir \code{Character} scalar. The path for reference data being installed to and storage.
#' @param tmpdir \code{Character} scalar. The temporary file storage path.
#' @param threads \code{Integer} scalar. The max threads allowed to be created.
#' @param adapter1 \code{Character} scalar. It is an adapter sequence for file1.
#' For single end data, it is requied.
#' @param adapter2 \code{Character} scalar. It is an adapter sequence for file2.
#' @param interleave \code{Logical} scalar. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param basicAnalysis \code{Logical} scalar. If it is TRUE, the pipeline will skip the time consuming steps
#' like GO annoation and motif analysis
#' @param createReport \code{Logical} scalar. If the HTML report file will be created.
#' @param motifs either\code{\link{PFMatrix}}, \code{\link{PFMatrixList}},
#' \code{\link{PWMatrix}}, \code{\link{PWMatrixList}}, default: vertebrates motif from JASPAR.
#' @param pipelineName \code{Character} scalar. Temporary file prefix for identifying files
#' when multiple pipeline generating file in the same tempdir.
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param p.cutoff p-value cutoff for returning motifs, default: 1e-6.
#' @param ... Additional arguments, currently unused.
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
#' \code{\link{atacMotifScan}},
#' \code{\link{atacRepsPipe}},
#' \code{\link{atacRepsPipe2}}
#' @examples
#' \dontrun{
#' ## These codes are time consuming so they will not be run and
#' ## checked by bioconductor checker.
#'
#'
#' # call pipeline
#' # for a quick example(only CTCF and BATF3 will be processing)
#' conclusion <-
#'   atacPipe(
#'        # MODIFY: Change these paths to your own case files!
#'        # e.g. fastqInput1 = "your/own/data/path.fastq"
#'        fastqInput1 = system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz"),
#'        fastqInput2 = system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz"),
#'        # MODIFY: Set the genome for your data
#'        genome = "hg19",
#'        motifs = getMotifInfo(motif.file = system.file("extdata", "CustomizedMotif.txt", package="esATAC")))
#'
#' # call pipeline
#' # for overall example(all vertebrates motif in JASPAR will be processed)
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

atacPipe <- function(genome, 
                     fastqInput1, 
                     fastqInput2=NULL, 
                     tmpdir = file.path(getwd(),"esATAC-pipeline"), 
                     refdir = file.path(tmpdir,"refdir"), 
                     threads = 2, 
                     adapter1 = NULL, 
                     adapter2 = NULL,
                     interleave = FALSE,  
                     basicAnalysis = FALSE, 
                     createReport = TRUE, 
                     motifs = NULL, 
                     pipelineName = "pipe",
                     chr = c(1:22, "X", "Y"), 
                     p.cutoff = 1e-6, ...){ #saveTmp = TRUE,
    
    prefix = pipelineName
    
    print(fastqInput1)
    print(fastqInput2)
        
    if(is.null(fastqInput2)&&!interleave&&is.null(adapter1)){
        stop("adapter1 should not be NULL for single end sequencing data")
    }
    if(!is.null(fastqInput1)){
        fastqInput1 = normalizePath(fastqInput1)
    }
    if(!is.null(fastqInput2)){
        fastqInput2 = normalizePath(fastqInput2)
    }
    
    dir.create(tmpdir,recursive = TRUE)
    setTmpDir(tmpdir)
    
    setJobName(pipelineName)
    
    setPipeName(pipelineName)
    
    setRefDir(refdir)
    
    setGenome(genome)
    setThreads(threads)
    unzipAndMerge <- atacUnzipAndMerge(fastqInput1 = fastqInput1,fastqInput2 = fastqInput2,interleave = interleave, ...)
    atacQC <- atacQCReport(atacProc = unzipAndMerge, ...)
    renamer <- atacRenamer(unzipAndMerge, ...)
    if(is.null(adapter1) || is.null(adapter2)){
        findAdapterObj <-atacFindAdapter(renamer, ...)
        removeAdapter <- atacRemoveAdapter(findAdapterObj, ...)
    }else{
        removeAdapter <- atacRemoveAdapter(renamer, adapter1 = adapter1, adapter2 = adapter2, ...)
    }
    bowtie2Mapping <- atacBowtie2Mapping(removeAdapter, ...)
    libComplexQC <- atacLibComplexQC(bowtie2Mapping, ...)
    samToBam <- atacSam2Bam(bowtie2Mapping,...)
#    sam2Bed <-atacSamToBed(bowtie2Mapping,maxFragLen = 2000, ...)
    sam2Bed <- atacBam2Bed(samToBam, ...)
    bedToBigWig <- atacBedToBigWig(sam2Bed, ...)
    tssqc100 <-atacTSSQC(sam2Bed,fragLenRange = c(0,100), newStepType = "TSSQCNFR", ...)
    
    if(is.null(fastqInput2)&&!interleave){
        if(testPeakCallingMACS2()){
            peakCalling <- atacPeakCallingMACS2(sam2Bed, ...)
        }else{
            peakCalling <- atacPeakCalling(sam2Bed, ...)
        }
        DHSQC <- atacPeakQC(peakCalling,qcbedInput = "DHS", ...)
        fripQC <- atacFripQC(atacProc = sam2Bed,atacProcPeak = peakCalling, ...)
    }else{
        tssqc180_247 <- atacTSSQC(sam2Bed,fragLenRange = c(180,247),  newStepType = "TSSQCneucleosome", ...)
        fragLenDistr <- atacFragLenDistr(sam2Bed, ...)
        shortBed <- atacBedUtils(sam2Bed,maxFragLen = 100, chrFilterList = NULL, ...)
        if(testPeakCallingMACS2()){
            peakCalling <- atacPeakCallingMACS2(sam2Bed, ...)
        }else{
            peakCalling <- atacPeakCalling(sam2Bed, ...)
        }
        DHSQC <- atacPeakQC(peakCalling,qcbedInput = "DHS",newStepType = "PeakQCDHS", ...)
        blacklistQC <- atacPeakQC(peakCalling,qcbedInput = "blacklist", newStepType = "PeakQCblacklist", ...)
        fripQC <- atacFripQC(atacProc = shortBed,atacProcPeak = peakCalling, ...)
        
        Peakanno <- atacPeakAnno(atacProc = peakCalling, ...)
        if(!basicAnalysis){
            # GO term
            goAna <- atacGOAnalysis(atacProc = Peakanno, ont = "BP", pvalueCutoff = 0.01, ...)
            
            # set default motif
            if(is.null(motifs)){
                opts <- list()
                opts[["tax_group"]] <- "vertebrates"
                pfm <- getMatrixSet(JASPAR2018::JASPAR2018, opts)
                names(pfm) <- TFBSTools::name(pfm)
                names(pfm) <- gsub(pattern = "[^a-zA-Z0-9]", replacement = "", x = TFBSTools::name(pfm), perl = TRUE)
            }else{
                pfm <- motifs
            }
            output_motifscan <- atacMotifScan(atacProc = peakCalling, motifs = pfm, p.cutoff = p.cutoff, prefix = prefix, ...)
            
            cs_output <- atacExtractCutSite(atacProc = sam2Bed, prefix = prefix, ...)
            footprint <- atacCutSiteCount(atacProcCutSite = cs_output, atacProcMotifScan = output_motifscan,
                                          strandLength = 100, prefix = prefix, chr = chr, ...)
            robj <- atacSingleRepReport(footprint, ...)
        }
    }

    reportDir <- file.path(getTmpDir(), paste0(pipelineName,"_report"))
    dir.create(reportDir)
    file.copy(from = output(peakCalling)$bedOutput, to=file.path(reportDir,"peak.bed"))
    file.copy(from = output(sam2Bed)$bedOutput, to=file.path(reportDir, "frag.bed"))
    file.copy(from = output(samToBam)$bamOutput, to=file.path(reportDir, "reads.bam"))
    file.copy(from = paste0(output(samToBam)$bamOutput,'.bai'), to=file.path(reportDir, "reads.bam.bai"))
}



#' @name atacPipe2
#' @title Pipeline for single replicate case-control paired-end sequencing data
#' @description
#' The preset pipeline to process case control study sequencing data.
#' An HTML report file, result files(e.g. BED, BAM files) and
#' conclusion list will generated. See detail for usage.
#' @param genome \code{Character} scalar. The genome(like hg19, mm10, etc.) reference data in "refdir" to be used in the pipeline.
#' @param case \code{List} scalar. Input for case sample. \code{fastqInput1},
#' the path(s) of the mate 1 fastq file(s), is required. \code{fastqInput2},
#' the path(s) of the mate 2 fastq file(s), is required, when \code{interleave=FALSE}.
#' \code{adapter1} and \code{adapter2} are optional.
#' @param control \code{List} scalar. Input for control sample. \code{fastqInput1},
#' the path(s) of the mate 1 fastq file(s), is required. \code{fastqInput2},
#' the path(s) of the mate 2 fastq file(s), is required, when \code{interleave=FALSE}.
#' \code{adapter1} and \code{adapter2} are optional.
#' @param refdir \code{Character} scalar. The path for reference data being installed to and storage.
#' @param tmpdir \code{Character} scalar. The temporary file storage path.
#' @param threads \code{Integer} scalar. The max threads allowed to be created.
#' @param interleave \code{Logical} scalar. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param createReport \code{Logical} scalar. If the HTML report file will be created.
#' @param motifs either\code{\link{PFMatrix}}, \code{\link{PFMatrixList}},
#' \code{\link{PWMatrix}}, \code{\link{PWMatrixList}}, default: vertebrates motif from JASPAR.
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param p.cutoff p-value cutoff for returning motifs, default: 1e-6.
#' @param ... Additional arguments, currently unused.
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
#' @import JASPAR2018
#' @importFrom TFBSTools getMatrixSet
#' @importFrom TFBSTools name
#' @examples
#' \dontrun{
#' ## These codes are time consuming so they will not be run and
#' ## checked by bioconductor checker.
#'
#'
#' # call pipeline
#' # for a quick example(only CTCF and BATF3 will be processed)
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
#'        motifs = getMotifInfo(motif.file = system.file("extdata", "CustomizedMotif.txt", package="esATAC")))
#'
#' # call pipeline
#' # for overall example(all vertebrates motif in JASPAR will be processed)
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
atacPipe2 <- function(genome, case = list(fastqInput1="paths/To/fastq1",fastqInput2="paths/To/fastq2", adapter1 = NULL, adapter2 = NULL),
                      control =list(fastqInput1="paths/To/fastq1",fastqInput2="paths/To/fastq2", adapter1 = NULL, adapter2 = NULL),
                      tmpdir = file.path(getwd(),"esATAC-pipeline"), 
                      refdir = file.path(tmpdir,"refdir"),  threads=2, interleave = FALSE, createReport = TRUE, motifs = NULL,
                      chr = c(1:22, "X", "Y"), p.cutoff = 1e-6, ...){ #saveTmp = TRUE,
    #stop("this function is still under developing")
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
    
    casePipe <- "CasePipe"
    # case
    atacPipe(genome=genome, 
             fastqInput1=case[["fastqInput1"]], 
             fastqInput2=case[["fastqInput2"]], 
             tmpdir = tmpdir, 
             refdir = refdir, 
             threads = threads, 
             adapter1 = case[["adapter1"]], 
             adapter2 = case[["adapter2"]],
             interleave = FALSE,  
             basicAnalysis = FALSE, 
             createReport = TRUE, 
             motifs = motifs, 
             pipelineName = casePipe,
             chr = chr, 
             p.cutoff = p.cutoff, ...)
    
    controlPipe <- "ControlPipe"
    #control
    atacPipe(genome=genome, 
             fastqInput1=control[["fastqInput1"]], 
             fastqInput2=control[["fastqInput2"]], 
             tmpdir = tmpdir, 
             refdir = refdir, 
             threads = threads, 
             adapter1 = control[["adapter1"]], 
             adapter2 = control[["adapter2"]],
             interleave = FALSE,  
             basicAnalysis = FALSE, 
             createReport = TRUE, 
             motifs = motifs, 
             pipelineName = controlPipe,
             chr = chr, 
             p.cutoff = p.cutoff, ...)
    
    casePipeObj <-getObjsInPipe(casePipe)
    controlPipeObj <- getObjsInPipe(controlPipe)
    
    setJobName("CaseControlCmp")
    setPipeName("CaseControlCmp")
    


    
    message("Begin to generate summary")
    case.frag.obj <- casePipeObj[[paste0(casePipe,'_BamToBed')]]
    ctrl.frag.obj <- controlPipeObj[[paste0(controlPipe,'_BamToBed')]]
    
    case.peak.obj <- casePipeObj[[paste0(casePipe,'_PeakCallingMACS2')]]
    ctrl.peak.obj <- controlPipeObj[[paste0(controlPipe,'_PeakCallingMACS2')]]
    
    
    
    peak.cmp.obj <- atacPeakComp(case.peak.obj,ctrl.peak.obj)
    
    if(testPeakCallingMACS2()){
        atacPeakCallingMACS2(case.peak.obj)
    }else{
        atacPeakCalling()
    }
    
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
    
    # motif analysis
    mout <- atacMotifScanPair(atacProc = peakCom, motifs = pwm, p.cutoff = p.cutoff)
    cs_case <- extractcutsite(bedInput = bed.case, prefix = "CASE")
    cs_ctrl <- extractcutsite(bedInput = bed.ctrl, prefix = "CONTROL")
    footprint.case <- atacCutSiteCount(atacProcCutSite = cs_case,
                                       motif_info = getParam(mout, "rdsOutput.peak1"),
                                       strandLength = 100, prefix = "Case_specific", chr = chr)
    footprint.ctrl <- atacCutSiteCount(atacProcCutSite = cs_ctrl,
                                       motif_info = getParam(mout, "rdsOutput.peak2"),
                                       strandLength = 100, prefix = "Control_specific", chr = chr)
    
    comp_result <- list(
        Peakanno.case = Peakanno.case,
        Peakanno.ctrl = Peakanno.ctrl,
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
    message("Generate summary done")
    if(createReport){
        message("Begin to generate Report")
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
        message("Generate report done")
        if(!(!is.null(param.tmp[["dontSet"]])&&param.tmp[["dontSet"]])){
            file.copy(file.path(.obtainConfigure("tmpdir"),"Report2.html"),esATAC_report, overwrite = TRUE)
            file.copy(getReportVal(caselist$atacProcs$atacQC,"pdf"),esATAC_report, overwrite = TRUE)
            file.copy(getReportVal(ctrllist$atacProcs$atacQC,"pdf"),esATAC_report, overwrite = TRUE)
            file.copy(getReportVal(goAna.case,"goOutput"),esATAC_report, overwrite = TRUE)
            file.copy(getReportVal(goAna.ctrl,"goOutput"),esATAC_report, overwrite = TRUE)
            dir.create(file.path(esATAC_result,"peak"))
            file.copy(getParam(caselist$atacProcs$peakCalling,"bedOutput"),file.path(esATAC_result,"peak"), overwrite = TRUE)
            file.copy(getParam(ctrllist$atacProcs$peakCalling,"bedOutput"),file.path(esATAC_result,"peak"), overwrite = TRUE)
            
            file.copy(getReportVal(caselist$atacProcs$Peakanno,"annoOutput"), esATAC_result, overwrite = TRUE)
            file.copy(getReportVal(ctrllist$atacProcs$Peakanno,"annoOutput"), esATAC_result, overwrite = TRUE)
            file.copy(getReportVal(comp_result$Peakanno.case,"annoOutput"), esATAC_result, overwrite = TRUE)
            file.copy(getReportVal(comp_result$Peakanno.ctrl,"annoOutput"), esATAC_result, overwrite = TRUE)
            
            file.copy(from = getReportVal(caselist$atacProcs$footprint, "pdf.dir"),
                      to = esATAC_result, overwrite = TRUE, recursive = TRUE)
            file.copy(from = getReportVal(ctrllist$atacProcs$footprint, "pdf.dir"),
                      to = esATAC_result, overwrite = TRUE, recursive = TRUE)
            file.copy(from = getReportVal(comp_result$footprint.case, "pdf.dir"),
                      to = esATAC_result, overwrite = TRUE, recursive = TRUE)
            file.copy(from = getReportVal(comp_result$footprint.ctrl, "pdf.dir"),
                      to = esATAC_result, overwrite = TRUE, recursive = TRUE)
            # generate motif enrichment file
            motif_enrich.case <- getReportVal(comp_result$mout, "rdsOutput.peak1")
            motif_enrich.case <- motif_enrich.case[, c(1, 3, 4)]
            colnames(motif_enrich.case) <- c("motif", "motif length", "p_value")
            motif_enrich.case <- motif_enrich.case[order(motif_enrich.case$p_value), ]
            rownames(motif_enrich.case) <- seq(nrow(motif_enrich.case))
            motif_enrich.case_file <- paste(esATAC_result, "/motif_enrichment_Case.txt", sep = "")
            write.table(x = motif_enrich.case, file = motif_enrich.case_file, sep = "\t",
                        row.names = TRUE, col.names = TRUE)
            
            motif_enrich.ctrl <- getReportVal(comp_result$mout, "rdsOutput.peak2")
            motif_enrich.ctrl <- motif_enrich.ctrl[, c(1, 3, 4)]
            colnames(motif_enrich.ctrl) <- c("motif", "motif length", "p_value")
            motif_enrich.ctrl <- motif_enrich.ctrl[order(motif_enrich.ctrl$p_value), ]
            rownames(motif_enrich.ctrl) <- seq(nrow(motif_enrich.ctrl))
            motif_enrich.ctrl_file <- paste(esATAC_result, "/motif_enrichment_Control.txt", sep = "")
            write.table(x = motif_enrich.ctrl, file = motif_enrich.ctrl_file, sep = "\t",
                        row.names = TRUE, col.names = TRUE)
            message(sprintf("type `browseURL(\"%s\")` to view Report in web browser",file.path(esATAC_report,"Report2.html")))
        }else{
            message(sprintf("type `browseURL(\"%s\")` to view Report in web browser",file.path(.obtainConfigure("tmpdir"),"Report2.html")))
        }
    }
    
    

}


checkFilePathExist <- function(afilePaths){
    if(FALSE %in% file.exists(afilePaths)){
        stop(sprintf("file '%s' does not exist",afilePaths))
    }
}

#' @name atacRepsPipe
#' @title Pipeline for multi-replicates case paired-end sequencing data
#' @description
#' The preset pipeline to process multi-replicates case study sequencing data.
#' HTML report files, result files(e.g. BED, BAM files) and
#' conclusion list will generated. See detail for usage.
#' @param genome \code{Character} scalar. The genome(like hg19, mm10, etc.) reference data in "refdir" to be used in the pipeline.
#' @param fastqInput1 \code{List} scalar. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates paired
#' with file paths in fastqInput2
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}.
#' Each element in the fastqInput1 \code{List} is for a replicate
#' It can be a \code{Character} vector of FASTQ files paths to be merged.
#' @param fastqInput2 \code{List} scalar. It contains file paths with #2
#' mates paired with file paths in fastqInput1.
#' For single-end sequencing files and interleaved paired-end sequencing
#' files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' Each element in the fastqInput1 \code{List} is for a replicate
#' It can be a \code{Character} vector of FASTQ files paths to be merged.
#' @param refdir \code{Character} scalar. The path for reference data being installed to and storage.
#' @param tmpdir \code{Character} scalar. The temporary file storage path.
#' @param threads \code{Integer} scalar. The max threads allowed to be created.
#' @param adapter1 \code{Character} scalar. It is an adapter sequence for file1.
#' For single end data, it is requied.
#' @param adapter2 \code{Character} scalar. It is an adapter sequence for file2.
#' @param interleave \code{Logical} scalar. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param createReport \code{Logical} scalar. If the HTML report file will be created.
#' @param motifs either\code{\link{PFMatrix}}, \code{\link{PFMatrixList}},
#' \code{\link{PWMatrix}}, \code{\link{PWMatrixList}}, default: vertebrates motif from JASPAR.
#' @param prefix \code{Character} scalar. Temporary file prefix for identifying files
#' when multiple pipeline generating file in the same tempdir.
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param p.cutoff p-value cutoff for returning motifs, default: 1e-6.
#' @param ... Additional arguments, currently unused.
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
#' # for a quick example(only CTCF and BATF3 will be processing)
#' conclusion <-
#'   atacRepsPipe(
#'        # MODIFY: Change these paths to your own case files!
#'        # e.g. fastqInput1 = "your/own/data/path.fastq"
#'        fastqInput1 = list(system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz"),
#'                           system.file(package="esATAC", "extdata", "chr20_1.2.fq.bz2")),
#'        fastqInput2 = list(system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz"),
#'                           system.file(package="esATAC", "extdata", "chr20_2.2.fq.bz2")),
#'        # MODIFY: Set the genome for your data
#'        genome = "hg19",
#'        motifs = getMotifInfo(motif.file = system.file("extdata", "CustomizedMotif.txt", package="esATAC")))
#'
#' # call pipeline
#' # for overall example(all vertebrates motif in JASPAR will be processed)
#' conclusion <-
#'   atacRepsPipe(
#'        # MODIFY: Change these paths to your own case files!
#'        # e.g. fastqInput1 = "your/own/data/path.fastq"
#'        fastqInput1 = list(system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz"),
#'                           system.file(package="esATAC", "extdata", "chr20_1.2.fq.bz2")),
#'        fastqInput2 = list(system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz"),
#'                           system.file(package="esATAC", "extdata", "chr20_2.2.fq.bz2")),
#'        # MODIFY: Set the genome for your data
#'        genome = "hg19")
#' }
#' @importFrom VennDiagram venn.diagram
#' @importFrom corrplot corrplot
#' @importFrom grDevices colorRampPalette
#' @export
#'
atacRepsPipe <- function(genome, fastqInput1, fastqInput2 = NULL, refdir = NULL, tmpdir = NULL, threads = 2,
                         adapter1 = NULL, adapter2 = NULL, interleave = FALSE,  createReport = TRUE,
                         motifs = NULL, prefix = NULL, chr = c(1:22, "X", "Y"), p.cutoff = 1e-6, ...){
    stop("this function is still under developing")

}

#' @importFrom GenomeInfoDb seqlengths
#' @importFrom rtracklayer import.bed

getBinsReadsCount <- function(bedInput,bsgenome,binsize = 1000){
    abedfile <- import.bed(con = bedInput)
    chrominfo<-seqinfo(bsgenome)
    chroms <- seqnames(chrominfo)
    chromsize <- GenomeInfoDb::seqlengths(chrominfo)
    binnumb <- as.integer(chromsize/binsize)+1
    readsCountsList <- list()
    pos<-as.integer((end(ranges(abedfile))+start(ranges(abedfile)))/2/binsize) + 1
    for(i in 1:length(chroms)){
        readsCountsList[[chroms[i]]]<-rep(0,binnumb[i])
        posadd <- pos[as.character(seqnames(abedfile))==chroms[i]]
        readsCountsNumber<-table(posadd)
        readsCountsList[[chroms[i]]][as.integer(names(readsCountsNumber))] <- readsCountsList[[chroms[i]]][as.integer(names(readsCountsNumber))] + as.numeric(readsCountsNumber)
    }
    return(AnnotationDbi::unlist2(readsCountsList))
    #return(readsCountsList)
}



#' @name atacRepsPipe2
#' @title Pipeline for multi-replicates case-control paired-end sequencing data
#' @description
#' The preset pipeline to process multi-replicates case control study sequencing data.
#' HTML report files, result files(e.g. BED, BAM files) and
#' conclusion list will generated. See detail for usage.
#' @param genome \code{Character} scalar. The genome(like hg19, mm10, etc.) reference data in "refdir" to be used in the pipeline.
#' @param caseFastqInput1 \code{List} scalar. Input for case samples.
#' For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates paired
#' with file paths in fastqInput2
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}.
#' Each element in the caseFastqInput1 \code{List} is for a replicate
#' It can be a \code{Character} vector of FASTQ files paths to be merged.
#' @param caseFastqInput2 \code{List} scalar. Input for case samples.
#' It contains file paths with #2
#' mates paired with file paths in caseFastqInput1
#' For single-end sequencing files and interleaved paired-end sequencing
#' files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' Each element in the caseFastqInput2 \code{List} is for a replicate
#' @param ctrlFastqInput1 \code{List} scalar. Input for control samples.
#' For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates paired
#' with file paths in ctrlFastqInput2
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}.
#' Each element in the ctrlFastqInput1 \code{List} is for a replicate
#' It can be a \code{Character} vector of FASTQ files paths to be merged.
#' @param ctrlFastqInput2 \code{List} scalar. Input for control samples.
#' It contains file paths with #2
#' mates paired with file paths in fastqInput1.
#' For single-end sequencing files and interleaved paired-end sequencing
#' files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' Each element in the ctrlFastqInput1 \code{List} is for a replicate
#' @param caseAdapter1 \code{Character} scalar. Adapter for caseFastqInput1.
#' @param caseAdapter2 \code{Character} scalar. Adapter for caseFastqInput2.
#' @param ctrlAdapter1 \code{Character} scalar. Adapter for ctrlFastqInput1.
#' @param ctrlAdapter2 \code{Character} scalar. Adapter for ctrlFastqInput2.
#' @param refdir \code{Character} scalar. The path for reference data being installed to and storage.
#' @param tmpdir \code{Character} scalar. The temporary file storage path.
#' @param threads \code{Integer} scalar. The max threads allowed to be created.
#' @param interleave \code{Logical} scalar. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param createReport \code{Logical} scalar. If the HTML report file will be created.
#' @param motifs either\code{\link{PFMatrix}}, \code{\link{PFMatrixList}},
#' \code{\link{PWMatrix}}, \code{\link{PWMatrixList}}, default: vertebrates motif from JASPAR.
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param p.cutoff p-value cutoff for returning motifs, default: 1e-6.
#' @param ... Additional arguments, currently unused.
#' @return \code{List} scalar. It is a list that save the result of the pipeline.
#' Slot "caselist" and "ctrlist": Each of them is a list that save the result for case or control data.
#' Slot "comp_result": compare analysis result for case and control data
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
#' @import JASPAR2018
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
#'     atacRepsPipe2(
#'         # MODIFY: Change these paths to your own case files!
#'         # e.g. fastqInput1 = "your/own/data/path.fastq"
#'      caseFastqInput1=list(system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz"),
#'                           system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz")),
#'      # MODIFY: Change these paths to your own case files!
#'      # e.g. fastqInput1 = "your/own/data/path.fastq"
#'      caseFastqInput2=list(system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz"),
#'                           system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz")),
#'      # MODIFY: Change these paths to your own control files!
#'      # e.g. fastqInput1 = "your/own/data/path.fastq"
#'      ctrlFastqInput1=list(system.file(package="esATAC", "extdata", "chr20_1.2.fq.bz2"),
#'                           system.file(package="esATAC", "extdata", "chr20_1.2.fq.bz2")),
#'      # MODIFY: Change these paths to your own control files!
#'      # e.g. fastqInput1 = "your/own/data/path.fastq"
#'      ctrlFastqInput2=list(system.file(package="esATAC", "extdata", "chr20_2.2.fq.bz2"),
#'                           system.file(package="esATAC", "extdata", "chr20_2.2.fq.bz2")),
#'      # MODIFY: Set the genome for your data
#'      genome = "hg19",
#'      motifs = getMotifInfo(motif.file = system.file("extdata", "CustomizedMotif.txt", package="esATAC")))
#'
#'
#' # call pipeline
#' # for overall example(all human motif in JASPAR will be processed)
#' conclusion <-
#'     atacRepsPipe2(
#'         # MODIFY: Change these paths to your own case files!
#'         # e.g. fastqInput1 = "your/own/data/path.fastq"
#'      caseFastqInput1=list(system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz"),
#'                           system.file(package="esATAC", "extdata", "chr20_1.1.fq.gz")),
#'      # MODIFY: Change these paths to your own case files!
#'      # e.g. fastqInput1 = "your/own/data/path.fastq"
#'      caseFastqInput2=list(system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz"),
#'                           system.file(package="esATAC", "extdata", "chr20_2.1.fq.gz")),
#'      # MODIFY: Change these paths to your own control files!
#'      # e.g. fastqInput1 = "your/own/data/path.fastq"
#'      ctrlFastqInput1=list(system.file(package="esATAC", "extdata", "chr20_1.2.fq.bz2"),
#'                           system.file(package="esATAC", "extdata", "chr20_1.2.fq.bz2")),
#'      # MODIFY: Change these paths to your own control files!
#'      # e.g. fastqInput1 = "your/own/data/path.fastq"
#'      ctrlFastqInput2=list(system.file(package="esATAC", "extdata", "chr20_2.2.fq.bz2"),
#'                           system.file(package="esATAC", "extdata", "chr20_2.2.fq.bz2")),
#'      # MODIFY: Set the genome for your data
#'      genome = "hg19"
#'      )
#'}
#' @export
#'
atacRepsPipe2 <- function(genome, caseFastqInput1,caseFastqInput2, ctrlFastqInput1, ctrlFastqInput2,
                          caseAdapter1 = NULL, caseAdapter2 = NULL, ctrlAdapter1 = NULL, ctrlAdapter2 = NULL,
                          refdir=NULL, tmpdir=NULL, threads=2, interleave = FALSE, createReport = TRUE, motifs = NULL,
                          chr = c(1:22, "X", "Y"), p.cutoff = 1e-6, ...){ #saveTmp = TRUE,

    stop("this function is still under developing")

}


