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
#' @param motifPWM \code{List} scalar. Motif PWM, a list, default:vertebrates(JASPAR).
#' @param prefix \code{Character} scalar. Temporary file prefix for identifying files
#' when multiple pipeline generating file in the same tempdir.
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param min.score The minimum score for counting a match. Can be given as a
#' character string containing a percentage (default: "90%") of the highest
#' possible score or as a single number.
#' @param use.SavedPWM use local motif position information. This data is
#' download or generate by users. it must be a rds file and the information
#' saved as GRanges. Once this parameter is used, parameters "motifPWM" will be set to NULL.
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

atacPipe <- function(genome, fastqInput1, fastqInput2=NULL, refdir=NULL, tmpdir=NULL, threads=2, adapter1 = NULL, adapter2 = NULL,
                     interleave = FALSE,  basicAnalysis = FALSE, createReport = TRUE, motifPWM = NULL, prefix = NULL,
                     chr = c(1:22, "X", "Y"), min.score = "90%", use.SavedPWM = NULL, ...){ #saveTmp = TRUE,
    
    if(is.null(fastqInput2)&&!interleave&&is.null(adapter1)){
        stop("adapter1 should not be NULL for single end sequencing data")
    }
    if(!is.null(fastqInput1)){
        fastqInput1 = normalizePath(fastqInput1)
    }
    if(!is.null(fastqInput2)){
        fastqInput2 = normalizePath(fastqInput2)
    }
    
    esATAC_result <- NULL
    esATAC_report <- NULL
    param.tmp <- list(...)
    if(!(!is.null(param.tmp[["dontSet"]])&&param.tmp[["dontSet"]])){
        if(!is.null(refdir)){
            options(atacConf=setConfigure("refdir",refdir))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","refdir"))){
                dir.create(file.path("esATAC_pipeline","refdir"))
            }
            options(atacConf=setConfigure("refdir",file.path("esATAC_pipeline","refdir")))
        }
        if(!is.null(threads)){
            options(atacConf=setConfigure("threads",as.numeric(threads)))
            message(getConfigure("threads"))
        }else{
            options(atacConf=setConfigure("threads",2))
        }
        if(!is.null(genome)){
            options(atacConf=setConfigure("genome",genome))
        }else{
            stop("parameter genome is required")
        }
        if(!is.null(tmpdir)){
            options(atacConf=setConfigure("tmpdir",tmpdir))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","intermediate_results"))){
                dir.create(file.path("esATAC_pipeline","intermediate_results"))
            }else{
                warning(sprintf("path '%s' is exist",file.path("esATAC_pipeline","intermediate_results")))
            }
            options(atacConf=setConfigure("tmpdir",file.path("esATAC_pipeline","intermediate_results")))
            tmpdir <- "esATAC_pipeline"
        }
        if(is.null(param.tmp[["esATAC_result"]])){
            esATAC_result<-file.path(dirname(.obtainConfigure("tmpdir")),"esATAC_result")
            dir.create(esATAC_result)
        }else{
            esATAC_result <- param.tmp[["esATAC_result"]]
        }
        if(is.null(param.tmp[["esATAC_report"]])){
            esATAC_report<-file.path(dirname(.obtainConfigure("tmpdir")),"esATAC_report")
            dir.create(esATAC_report)
        }else{
            esATAC_report <- param.tmp[["esATAC_report"]]
        }
    }else{
        if(is.null(param.tmp[["esATAC_result"]])){
            esATAC_result<-NULL
        }else{
            esATAC_result <- param.tmp[["esATAC_result"]]
        }
        if(is.null(param.tmp[["esATAC_report"]])){
            esATAC_report<-NULL
        }else{
            esATAC_report <- param.tmp[["esATAC_report"]]
        }
    }
    message("final esATAC_result")
    message(esATAC_result)
    message("final esATAC_report")
    message(esATAC_report)
    
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
        
        Peakanno <- atacPeakAnno(atacProc = peakCalling)
        if(!basicAnalysis){
            goAna <- atacGOAnalysis(atacProc = Peakanno, ont = "BP", pvalueCutoff = 0.01)
            if(is.null(use.SavedPWM)){
                if(is.null(motifPWM)){
                    opts <- list()
                    opts[["tax_group"]] <- "vertebrates"
                    pwm <- getMatrixSet(JASPAR2016::JASPAR2016, opts)
                    pwm <- TFBSTools::toPWM(pwm)
                    names(pwm) <- TFBSTools::name(pwm)
                    pwm <- lapply(X = pwm, FUN = TFBSTools::as.matrix)
                    names(pwm) <- gsub(pattern = "[^a-zA-Z0-9]", replacement = "", x = names(pwm), perl = TRUE)
                }else{
                    pwm <- motifPWM
                }
                output_motifscan <- atacMotifScan(atacProc = peakCalling, motifPWM = pwm, min.score = min.score, prefix = prefix)
            }else{
                output_motifscan <- atacMotifScan(atacProc = peakCalling, use.SavedPWM = use.SavedPWM, prefix = prefix)
            }
            
            
            
            cs_output <- atacExtractCutSite(atacProc = sam2Bed, prefix = prefix)
            footprint <- atacCutSiteCount(atacProcCutSite = cs_output, atacProcMotifScan = output_motifscan,
                                          strandLength = 100, prefix = prefix, chr = chr)
        }
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
    if(basicAnalysis){
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
                       atacQC = atacQC
        )
    }else{
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
    }
    
    conclusion <- list(filelist=filelist,
                       wholesummary = wholesummary,
                       atacProcs = atacProcs,
                       filtstat = filtstat
    )
    
    if(createReport){
        message("Begin to generate report")
        filename <- strsplit(fastqInput1,".fastq|.FASTQ|.FQ|.fq")[[1]][1]
        filename <- basename(filename)
        
        if(basicAnalysis){
            rmdfile<-system.file(package="esATAC", "extdata", "basicReport.Rmd")
        }else{
            rmdfile<-system.file(package="esATAC", "extdata", "Report.Rmd")
        }
        
        rmdtext<-readChar(rmdfile,nchars=file.info(rmdfile)$size,useBytes = TRUE)
        #rmdtext<-sprintf(rmdtext,filename)
        
        workdir <- getwd()
        save(filelist,wholesummary,filtstat,atacProcs,workdir,file = file.path(.obtainConfigure("tmpdir"),"Report.Rdata"))
        
        writeChar(rmdtext,con = file.path(.obtainConfigure("tmpdir"),"Report.Rmd"),useBytes = TRUE)
        render(file.path(.obtainConfigure("tmpdir"),"Report.Rmd"))
        #knit(file.path(.obtainConfigure("tmpdir"),"Report.Rmd"), file.path(.obtainConfigure("tmpdir"),"Report.md"))
        #markdownToHTML(file.path(.obtainConfigure("tmpdir"),"Report.md"), file.path(.obtainConfigure("tmpdir"),"Report.html"))
        #browseURL(paste0('file://', file.path(.obtainConfigure("tmpdir"),"Report.html")))
        message("Generate report done")
        
        if(!is.null(esATAC_report)&&!is.null(esATAC_result)){
            file.copy(file.path(.obtainConfigure("tmpdir"),"Report.html"),esATAC_report, overwrite = TRUE)
            file.copy(getReportVal(atacQC,"pdf"),esATAC_report, overwrite = TRUE)
            dir.create(file.path(esATAC_result,"peak"))
            file.copy(getParam(peakCalling,"bedOutput"),file.path(esATAC_result,"peak"), overwrite = TRUE)
            file.copy(getReportVal(Peakanno,"annoOutput"),esATAC_result, overwrite = TRUE)
            if(!basicAnalysis){
                file.copy(getReportVal(goAna,"goOutput"),esATAC_report, overwrite = TRUE)
                file.copy(from = getReportVal(atacProcs$footprint,"pdf.dir"),
                          to = esATAC_result, overwrite = TRUE, recursive = TRUE)
            }
            message(sprintf("type `browseURL(\"%s\")` to view Report in web browser",file.path(esATAC_report,"Report.html")))
        }else{
            message(sprintf("type `browseURL(\"%s\")` to view Report in web browser",file.path(.obtainConfigure("tmpdir"),"Report.html")))
        }
    }
    
    invisible(conclusion)
    
    
    
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
#' @param motifPWM \code{List} scalar. Motif PWM, a list, default:vertebrates(JASPAR).
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param min.score The minimum score for counting a match. Can be given as a
#' character string containing a percentage (default: "90%") of the highest
#' possible score or as a single number.
#' @param use.SavedPWM use local motif position information. This data is
#' download or generate by users. it must be a rds file and the information
#' saved as GRanges. Once this parameter is used, parameters "motifPWM" will be set to NULL.
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
atacPipe2 <- function(genome, case = list(fastqInput1="paths/To/fastq1",fastqInput2="paths/To/fastq2", adapter1 = NULL, adapter2 = NULL),
                      control =list(fastqInput1="paths/To/fastq1",fastqInput2="paths/To/fastq2", adapter1 = NULL, adapter2 = NULL),
                      refdir=NULL, tmpdir=NULL, threads=2,
                      interleave = FALSE, createReport = TRUE, motifPWM = NULL, chr = c(1:22, "X", "Y"), min.score = "90%",
                      use.SavedPWM = NULL,  ...){ #saveTmp = TRUE,
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
        if(!is.null(refdir)){
            options(atacConf=setConfigure("refdir",refdir))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","refdir"))){
                dir.create(file.path("esATAC_pipeline","refdir"))
            }
            options(atacConf=setConfigure("refdir",file.path("esATAC_pipeline","refdir")))
        }
        if(!is.null(threads)){
            options(atacConf=setConfigure("threads",as.numeric(threads)))
            message(getConfigure("threads"))
        }else{
            options(atacConf=setConfigure("threads",2))
        }
        if(!is.null(genome)){
            options(atacConf=setConfigure("genome",genome))
        }else{
            stop("parameter genome is required")
        }
        if(!is.null(tmpdir)){
            options(atacConf=setConfigure("tmpdir",tmpdir))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","intermediate_results"))){
                dir.create(file.path("esATAC_pipeline","intermediate_results"))
            }else{
                warning(sprintf("path '%s' is exist",file.path("esATAC_pipeline","intermediate_results")))
            }
            options(atacConf=setConfigure("tmpdir",file.path("esATAC_pipeline","intermediate_results")))
            tmpdir <- "esATAC_pipeline"
        }
        if(is.null(param.tmp[["esATAC_result"]])){
            esATAC_result<-file.path(dirname(.obtainConfigure("tmpdir")),"esATAC_result")
            dir.create(esATAC_result)
        }else{
            esATAC_result <- param.tmp[["esATAC_result"]]
        }
        if(is.null(param.tmp[["esATAC_report"]])){
            esATAC_report<-file.path(dirname(.obtainConfigure("tmpdir")),"esATAC_report")
            dir.create(esATAC_report)
        }else{
            esATAC_report <- param.tmp[["esATAC_report"]]
        }
    }else{
        if(is.null(param.tmp[["esATAC_result"]])){
            esATAC_result<-NULL
        }else{
            esATAC_result <- param.tmp[["esATAC_result"]]
        }
        if(is.null(param.tmp[["esATAC_report"]])){
            esATAC_report<-NULL
        }else{
            esATAC_report <- param.tmp[["esATAC_report"]]
        }
    }
    
    if(is.null(use.SavedPWM)){
        if(is.null(motifPWM)){
            opts <- list()
            opts[["tax_group"]] <- "vertebrates"
            pwm <- getMatrixSet(JASPAR2016::JASPAR2016, opts)
            pwm <- TFBSTools::toPWM(pwm)
            names(pwm) <- TFBSTools::name(pwm)
            pwm <- lapply(X = pwm, FUN = TFBSTools::as.matrix)
            names(pwm) <- gsub(pattern = "[^a-zA-Z0-9]", replacement = "", x = names(pwm), perl = TRUE)
        }else{
            pwm <- motifPWM
        }
    }
    
    message("Begin to process case sample...")
    if(is.null(use.SavedPWM)){
        caselist <- atacPipe(fastqInput1 = case[["fastqInput1"]],fastqInput2 = case[["fastqInput2"]],
                             adapter1 = case[["adapter1"]], adapter2 = case[["adapter2"]],interleave = interleave,
                             createReport = FALSE, motifPWM =pwm, prefix = "Case_data", chr = chr,
                             min.score = min.score, dontSet=TRUE) #saveTmp = TRUE,
    }else{
        caselist <- atacPipe(fastqInput1 = case[["fastqInput1"]],fastqInput2 = case[["fastqInput2"]],
                             adapter1 = case[["adapter1"]], adapter2 = case[["adapter2"]],interleave = interleave,
                             createReport = FALSE, prefix = "Case_data", chr = chr, use.SavedPWM = use.SavedPWM,
                             dontSet=TRUE) #saveTmp = TRUE,
    }
    
    message("Case sample process done")
    message(" ")
    message("Begin to process control sample")
    if(is.null(use.SavedPWM)){
        ctrllist <- atacPipe(fastqInput1 = control[["fastqInput1"]],fastqInput2 = control[["fastqInput2"]],
                             adapter1 = control[["adapter1"]], adapter2 = control[["adapter2"]],interleave = interleave,
                             createReport = FALSE, motifPWM =pwm, prefix = "Control_data", chr = chr,
                             min.score = min.score, dontSet=TRUE) #saveTmp = TRUE,
    }else{
        ctrllist <- atacPipe(fastqInput1 = control[["fastqInput1"]],fastqInput2 = control[["fastqInput2"]],
                             adapter1 = control[["adapter1"]], adapter2 = control[["adapter2"]],interleave = interleave,
                             createReport = FALSE, prefix = "Control_data", chr = chr, use.SavedPWM = use.SavedPWM,
                             dontSet=TRUE) #saveTmp = TRUE,
    }
    
    message("control sample process done")
    message(" ")
    message("Begin to generate summary")
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
    if(is.null(use.SavedPWM)){
        mout <- atacMotifScanPair(atacProc = peakCom, motifPWM = pwm, min.score = min.score)
    }else{
        mout <- atacMotifScanPair(atacProc = peakCom, use.SavedPWM = use.SavedPWM)
    }
    
    
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
            file.copy(getReportVal(caselist$atacProcs$goAna,"goOutput"),esATAC_report, overwrite = TRUE)
            file.copy(getReportVal(ctrllist$atacProcs$goAna,"goOutput"),esATAC_report, overwrite = TRUE)
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
    
    invisible(conclusion)
    
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
#' @param motifPWM \code{List} scalar. Motif PWM, a list.
#' @param prefix \code{Character} scalar. Temporary file prefix for identifying files
#' when multiple pipeline generating file in the same tempdir.
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param min.score The minimum score for counting a match. Can be given as a
#' character string containing a percentage (default: "90%") of the highest
#' possible score or as a single number.
#' @param use.SavedPWM use local motif position information. This data is
#' download or generate by users. it must be a rds file and the information
#' saved as GRanges. Once this parameter is used, parameters "motifPWM" will be set to NULL.
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
#' # for a quick example(only CTCF will be processing)
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
#'        motifPWM = getMotifPWM(motif.file = system.file("extdata", "CTCF.txt", package="esATAC"), is.PWM = FALSE))
#'
#' # call pipeline
#' # for overall example(all human motif in JASPAR will be processed)
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
#' @export
#'
atacRepsPipe <- function(genome, fastqInput1,fastqInput2=NULL, refdir=NULL, tmpdir=NULL, threads=2, adapter1 = NULL, adapter2 = NULL,
                         interleave = FALSE,  createReport = TRUE, motifPWM = NULL, prefix = NULL,
                         chr = c(1:22, "X", "Y"), min.score = "90%", use.SavedPWM = NULL, ...){
    if(is.null(fastqInput2)&&!interleave&&is.null(adapter1)){
        stop("adapter1 should not be NULL for single end sequencing data")
    }
    if(!is.null(fastqInput1)){
        if(class(fastqInput1)=="list"){
            for(n in 1:length(fastqInput1)){
                fastqInput1[[n]] <- normalizePath(fastqInput1[[n]])
                checkFilePathExist(fastqInput1[[n]])
            }
            if(!is.null(adapter1)){
                if(length(adapter1)!=length(fastqInput1)){
                    stop("the number of adapters in adapter1 do not match the number of replicates in fastqInput1")
                }
            }
        }else{
            stop("fastqInput1 must be a list")
        }
        
    }else{
        stop("fastqInput1 is required!")
    }
    
    if(!is.null(fastqInput2)){
        if(class(fastqInput2)=="list"){
            if(length(fastqInput1)!=length(fastqInput2)){
                stop(sprintf("fastqInput1 replicates size: %d and fastqInput2 replicates size %d do not match!",length(fastqInput1),length(fastqInput2)))
            }
            for(n in 1:length(fastqInput2)){
                if(length(fastqInput1[[n]])!=length(fastqInput2[[n]])){
                    stop(sprintf("fastq file number does not match in fastqInput1 and fastqInput2 replicates %d",n))
                }
                fastqInput2[[n]] <- normalizePath(fastqInput2[[n]])
                checkFilePathExist(fastqInput2[[n]])
            }
            if(!is.null(adapter2)){
                if(length(adapter2)!=length(fastqInput2)){
                    stop("the number of adapters in adapter2 do not match the number of replicates in fastqInput2")
                }
            }
        }else{
            stop("fastqInput2 must be a list")
        }
    }
    
    if(is.null(motifPWM)){
        opts <- list()
        opts[["tax_group"]] <- "vertebrates"
        pwm <- getMatrixSet(JASPAR2016::JASPAR2016, opts)
        pwm <- TFBSTools::toPWM(pwm)
        names(pwm) <- TFBSTools::name(pwm)
        pwm <- lapply(X = pwm, FUN = TFBSTools::as.matrix)
        names(pwm) <- gsub(pattern = "[^a-zA-Z0-9]", replacement = "", x = names(pwm), perl = TRUE)
    }else{
        pwm <- motifPWM
    }
    
    param.tmp <- list(...)
    if(!(!is.null(param.tmp[["dontSet"]])&&param.tmp[["dontSet"]])){
        if(!is.null(refdir)){
            options(atacConf=setConfigure("refdir",refdir))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","refdir"))){
                dir.create(file.path("esATAC_pipeline","refdir"))
            }
            options(atacConf=setConfigure("refdir",file.path("esATAC_pipeline","refdir")))
        }
        if(!is.null(threads)){
            options(atacConf=setConfigure("threads",as.numeric(threads)))
            message(getConfigure("threads"))
        }else{
            options(atacConf=setConfigure("threads",2))
        }
        if(!is.null(genome)){
            options(atacConf=setConfigure("genome",genome))
        }else{
            stop("parameter genome is required")
        }
        if(!is.null(tmpdir)){
            options(atacConf=setConfigure("tmpdir",tmpdir))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","intermediate_results"))){
                dir.create(file.path("esATAC_pipeline","intermediate_results"))
            }else{
                warning(sprintf("path '%s' is exist",file.path("esATAC_pipeline","intermediate_results")))
            }
            options(atacConf=setConfigure("tmpdir",file.path("esATAC_pipeline","intermediate_results")))
            tmpdir <- "esATAC_pipeline"
        }
        if(is.null(param.tmp[["esATAC_result"]])){
            esATAC_result<-file.path(dirname(.obtainConfigure("tmpdir")),"esATAC_result")
            dir.create(esATAC_result)
        }else{
            esATAC_result<-param.tmp[["esATAC_result"]]
        }
        if(is.null(param.tmp[["esATAC_report"]])){
            esATAC_report<-file.path(dirname(.obtainConfigure("tmpdir")),"esATAC_report")
            dir.create(esATAC_report)
        }else{
            esATAC_report<-param.tmp[["esATAC_report"]]
        }
    }else{
        if(is.null(param.tmp[["esATAC_result"]])){
            esATAC_result<-NULL
        }else{
            esATAC_result<-param.tmp[["esATAC_result"]]
        }
        if(is.null(param.tmp[["esATAC_report"]])){
            esATAC_report<-NULL
        }else{
            esATAC_report<-param.tmp[["esATAC_report"]]
        }
    }
    tmpdir <- .obtainConfigure("tmpdir")
    
    
    conclusions <- list()
    
    for(n in 1:length(fastqInput1)){
        if(!dir.exists(file.path(tmpdir,sprintf("replicate%d",n)))){
            dir.create(file.path(tmpdir,sprintf("replicate%d",n)))
        }
        options(atacConf=setConfigure("tmpdir",file.path(tmpdir,sprintf("replicate%d",n))))
        if(!is.null(fastqInput2)){
            fastqInput2n <- fastqInput2[[n]]
        }else{
            fastqInput2n <- NULL
        }
        if(!is.null(adapter1)){
            adapter1n <- adapter1[n]
        }else{
            adapter1n <- NULL
        }
        if(!is.null(adapter2)){
            adapter2n <- adapter2[n]
        }else{
            adapter2n <- NULL
        }
        if(!is.null(esATAC_report)&&!is.null(esATAC_result)){
            dir.create(file.path(esATAC_result,sprintf("replicate%d",n)))
            dir.create(file.path(esATAC_report,sprintf("replicate%d",n)))
            if(is.null(use.SavedPWM)){
                conclusions[[n]]<-atacPipe(genome = genome, fastqInput1 = fastqInput1[[n]],fastqInput2 = fastqInput2n, refdir = refdir,
                                           tmpdir = tmpdir, threads = threads, adapter1 = adapter1n, adapter2 = adapter2n,
                                           interleave = interleave,  basicAnalysis = TRUE, createReport = createReport, motifPWM = motifPWM, prefix = prefix,
                                           chr = chr, min.score = min.score, dontSet=TRUE,
                                           esATAC_result=normalizePath(file.path(esATAC_result,sprintf("replicate%d",n))),
                                           esATAC_report=normalizePath(file.path(esATAC_report,sprintf("replicate%d",n))),...)
            }else{
                conclusions[[n]]<-atacPipe(genome = genome, fastqInput1 = fastqInput1[[n]],fastqInput2 = fastqInput2n, refdir = refdir,
                                           tmpdir = tmpdir, threads = threads, adapter1 = adapter1n, adapter2 = adapter2n,
                                           interleave = interleave,  basicAnalysis = TRUE, createReport = createReport, prefix = prefix,
                                           chr = chr, use.SavedPWM = use.SavedPWM, dontSet=TRUE,
                                           esATAC_result=normalizePath(file.path(esATAC_result,sprintf("replicate%d",n))),
                                           esATAC_report=normalizePath(file.path(esATAC_report,sprintf("replicate%d",n))),...)
            }
            
        }else{
            if(is.null(use.SavedPWM)){
                conclusions[[n]]<-atacPipe(genome = genome, fastqInput1 = fastqInput1[[n]],fastqInput2 = fastqInput2n, refdir = refdir,
                                           tmpdir = tmpdir, threads = threads, adapter1 = adapter1n, adapter2 = adapter2n,
                                           interleave = interleave,  basicAnalysis = TRUE, createReport = createReport, motifPWM = motifPWM, prefix = prefix,
                                           chr = chr, min.score = min.score, dontSet=TRUE,...)
            }else{
                conclusions[[n]]<-atacPipe(genome = genome, fastqInput1 = fastqInput1[[n]],fastqInput2 = fastqInput2n, refdir = refdir,
                                           tmpdir = tmpdir, threads = threads, adapter1 = adapter1n, adapter2 = adapter2n,
                                           interleave = interleave,  basicAnalysis = TRUE, createReport = createReport, prefix = prefix,
                                           chr = chr, use.SavedPWM = use.SavedPWM, dontSet=TRUE,...)
            }
        }
    }
    if(!dir.exists(file.path(tmpdir,"rep_concord_merge"))){
        dir.create(file.path(tmpdir,"rep_concord_merge"))
    }
    
    unionPeak <- GRanges()
    genome <- seqinfo(.obtainConfigure("bsgenome"))
    peaklist <- list()
    for(n in 1:length(fastqInput1)){
        peakn <-import.bed(getParam(conclusions[[n]][["atacProcs"]][["peakCalling"]],"bedOutput"),genome=genome)
        peaknumb <- length(peakn) / 10 #top 10%
        sortedscore <- sort(mcols(peakn)$score,decreasing = TRUE)
        lowscore <- sortedscore[peaknumb]
        peakn <- peakn[mcols(peakn)$score>=lowscore]
        peaklist[[n]] <- peakn
        unionPeak <- union(unionPeak,peakn)
    }
    mcols(unionPeak)<-list(number=1:length(unionPeak))
    peaknumset<-list()
    for(n in 1:length(peaklist)){
        peaknumset[[sprintf("replicate%d",n)]]<-subjectHits(findOverlaps(peaklist[[n]], unionPeak,ignore.strand = TRUE))
    }
    options(atacConf=setConfigure("tmpdir",file.path(tmpdir,"rep_concord_merge")))
    if(length(fastqInput1)>5){
        warning("'VennDiagram' will not be generated for more than 5 replicates")
    }else{
        venn.diagram(x = peaknumset,filename = file.path(.obtainConfigure("tmpdir"),"vennDiagram.tiff"))
    }
    
    binsList <- NULL
    for(n in 1:length(fastqInput1)){
        binsList<-cbind(binsList, getBinsReadsCount(bedInput = getParam(conclusions[[n]][["atacProcs"]][["sam2Bed"]],"bedOutput"),
                                                    bsgenome = .obtainConfigure("bsgenome"),binsize = 1000))
    }
    colnames(binsList) <- paste0("replicate_",1:length(fastqInput1))
    correlation <- cor(binsList)
    message("the correlation matrix:")
    print(correlation)
    
    pdf(filename="corrplot.pdf")
    corrplot(correlation, method = "color", addCoef.col = "grey")
    dev.off()
    
    mergedReadsBed <- file.path(.obtainConfigure("tmpdir"),"mergedReads.bed")
    file.create(mergedReadsBed)
    if(!is.null(fastqInput2)){
        for(n in 1:length(fastqInput1)){
            file.append(mergedReadsBed,
                        getParam(conclusions[[n]][["atacProcs"]][["sam2Bed"]],"bedOutput"))
        }
        sortedReadsBed <- bedUtils(bedInput = mergedReadsBed, bedOutput = paste0(mergedReadsBed,".sorted.bed"), chrFilterList = NULL,sortBed = TRUE)
        shortBed <- atacBedUtils(sortedReadsBed, maxFragLen = 100, chrFilterList = NULL)
        peakCalling <- atacPeakCalling(shortBed)
        
        
        bedToBigWig <- atacBedToBigWig(shortBed)
        tssqc100 <-atacTSSQC(sortedReadsBed,reportPrefix = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(mergedReadsBed),".tssqc100")),fragLenRange = c(0,100))
        tssqc180_247 <-atacTSSQC(sortedReadsBed,reportPrefix = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(mergedReadsBed),".tssqc180_247")),fragLenRange = c(180,247))
        fragLenDistr <- atacFragLenDistr(sortedReadsBed)
        
        DHSQC <- atacPeakQC(peakCalling,qcbedInput = "DHS",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(mergedReadsBed),".DHSQC")))
        blacklistQC <- atacPeakQC(peakCalling,qcbedInput = "blacklist",reportOutput = file.path(.obtainConfigure("tmpdir"),paste0(getSuffixlessFileName0(mergedReadsBed),".blacklistQC")))
        fripQC <- atacFripQC(atacProcReads = shortBed,atacProcPeak = peakCalling)
        
        Peakanno <- atacPeakAnno(atacProc = peakCalling)
        goAna <- atacGOAnalysis(atacProc = Peakanno, ont = "BP", pvalueCutoff = 0.01)
        if(is.null(use.SavedPWM)){
            output_motifscan <- atacMotifScan(atacProc = peakCalling, motifPWM = pwm, min.score = min.score, prefix = prefix)
        }else{
            output_motifscan <- atacMotifScan(atacProc = peakCalling, use.SavedPWM = use.SavedPWM, prefix = prefix)
        }
        cs_output <- atacExtractCutSite(atacProc = sortedReadsBed, prefix = prefix)
        footprint <- atacCutSiteCount(atacProcCutSite = cs_output, atacProcMotifScan = output_motifscan,
                                      strandLength = 100, prefix = prefix, chr = chr)
        
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
        
        
        
        atacProcs=list(
            bedToBigWig = bedToBigWig,
            tssqc100 = tssqc100,
            tssqc180_247 = tssqc180_247,
            fragLenDistr = fragLenDistr,
            peakCalling = peakCalling,
            DHSQC = DHSQC,
            blacklistQC = blacklistQC,
            fripQC = fripQC,
            sortedReadsBed = sortedReadsBed,
            shortBed = shortBed,
            Peakanno = Peakanno,
            goAna = goAna,
            output_motifscan = output_motifscan,
            cs_output = cs_output,
            footprint = footprint
        )
        conclusion <- list(
            conclusions = conclusions,
            peaknumset = peaknumset,
            correlation = correlation,
            binsList = binsList,
            atacProcs = atacProcs)
        
        if(createReport){
            message("Begin to generate report")
            
            
            rmdfile<-system.file(package="esATAC", "extdata", "Rep_Report.Rmd")
            rmdtext<-readChar(rmdfile,nchars=file.info(rmdfile)$size,useBytes = TRUE)
            #rmdtext<-sprintf(rmdtext,filename)
            
            workdir <- getwd()
            save(atacProcs,peaknumset,correlation,binsList,file = file.path(.obtainConfigure("tmpdir"),"Report.Rdata"))
            
            writeChar(rmdtext,con = file.path(.obtainConfigure("tmpdir"),"Rep_Report.Rmd"),useBytes = TRUE)
            render(file.path(.obtainConfigure("tmpdir"),"Rep_Report.Rmd"))
            #knit(file.path(.obtainConfigure("tmpdir"),"Report.Rmd"), file.path(.obtainConfigure("tmpdir"),"Report.md"))
            #markdownToHTML(file.path(.obtainConfigure("tmpdir"),"Report.md"), file.path(.obtainConfigure("tmpdir"),"Report.html"))
            #browseURL(paste0('file://', file.path(.obtainConfigure("tmpdir"),"Report.html")))
            message("Generate report done")
            
            
            
            
            if(!is.null(esATAC_report)&&!is.null(esATAC_result)){
                
                #generate html index file
                replicateNum <- paste0("replicate ",1:length(fastqInput1))
                urllink <- paste0("<a href='./replicate",1:length(fastqInput1),"/Report.html'>","replicate ",1:length(fastqInput1)," report link</a>")
                singleRep <- data.frame(Report_Name=replicateNum,Link=urllink)
                
                mergeConRep <- data.frame(Report_Name="concordance and merge",link="<a href='./rep_concord_merge/Rep_Report.html'>concordance and merge analysis report link</a>")
                save(singleRep,mergeConRep,file = file.path(.obtainConfigure("tmpdir"),"ReportIdx.Rdata"))
                rmdidxfile<-system.file(package="esATAC", "extdata", "Rep_Report_Index.Rmd")
                file.copy(rmdidxfile,.obtainConfigure("tmpdir"))
                render(file.path(.obtainConfigure("tmpdir"),"Rep_Report_Index.Rmd"))
                file.copy(file.path(.obtainConfigure("tmpdir"),"Rep_Report_Index.html"),file.path(esATAC_report,"Report.html"))
                
                ## copy other files
                esATAC_report <- file.path(esATAC_report,"rep_concord_merge")
                esATAC_result <- file.path(esATAC_result,"rep_concord_merge")
                dir.create(esATAC_report)
                dir.create(esATAC_result)
                file.copy(file.path(.obtainConfigure("tmpdir"),"Rep_Report.html"),esATAC_report, overwrite = TRUE)
                file.copy(getReportVal(goAna,"goOutput"),esATAC_report, overwrite = TRUE)
                dir.create(file.path(esATAC_result,"peak"))
                file.copy(getParam(peakCalling,"bedOutput"),file.path(esATAC_result,"peak"), overwrite = TRUE)
                file.copy(getReportVal(Peakanno,"annoOutput"),esATAC_result, overwrite = TRUE)
                file.copy(from = getReportVal(atacProcs$footprint,"pdf.dir"),
                          to = esATAC_result, overwrite = TRUE, recursive = TRUE)
                
                
                
                
                message(sprintf("type `browseURL(\"%s\")` to view Report in web browser",file.path(esATAC_report,"Report.html")))
            }else{
                message(sprintf("type `browseURL(\"%s\")` to view Report in web browser",file.path(.obtainConfigure("tmpdir"),"Report.html")))
            }
        }
    }
    
    
    invisible(conclusion)
    
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
#' @param motifPWM \code{List} scalar. Motif PWM list.
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param min.score The minimum score for counting a match. Can be given as a
#' character string containing a percentage (default: "90%") of the highest
#' possible score or as a single number.
#' @param use.SavedPWM use local motif position information. This data is
#' download or generate by users. it must be a rds file and the information
#' saved as GRanges. Once this parameter is used, parameters "motifPWM" will be set to NULL.
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
#'      motifPWM = getMotifPWM(motif.file = system.file("extdata", "CTCF.txt", package="esATAC"), is.PWM = FALSE))
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
                          refdir=NULL, tmpdir=NULL, threads=2, interleave = FALSE, createReport = TRUE, motifPWM = NULL,
                          chr = c(1:22, "X", "Y"), min.score = "90%", use.SavedPWM = NULL, ...){ #saveTmp = TRUE,
    
    stopifnot(is.list(caseFastqInput1))
    stopifnot(length(caseFastqInput1)>0)
    for(n in 1:length(caseFastqInput1)){
        caseFastqInput1[[n]] <- normalizePath(caseFastqInput1[[n]])
        checkFilePathExist(caseFastqInput1[[n]])
    }
    if(!is.null(caseAdapter1)){
        if(length(caseAdapter1)!=length(caseFastqInput1)){
            stop("the number of adapters in caseAdapter1 do not match the number of replicates in caseFastqInput1")
        }
    }
    
    stopifnot(is.list(caseFastqInput2))
    stopifnot(length(caseFastqInput2)>0)
    if(length(caseFastqInput1)!=length(caseFastqInput2)){
        stop(sprintf("caseFastqInput1 replicates size: %d and caseFastqInput2 replicates size %d do not match!",length(caseFastqInput1),length(caseFastqInput2)))
    }
    for(n in 1:length(caseFastqInput2)){
        if(length(caseFastqInput1[[n]])!=length(caseFastqInput2[[n]])){
            stop(sprintf("fastq file number does not match in caseFastqInput1 and caseFastqInput2 replicates %d",n))
        }
        caseFastqInput2[[n]] <- normalizePath(caseFastqInput2[[n]])
        checkFilePathExist(caseFastqInput2[[n]])
    }
    if(!is.null(caseAdapter1)){
        if(length(caseAdapter2)!=length(caseFastqInput2)){
            stop("the number of adapters in caseAdapter2 do not match the number of replicates in caseFastqInput2")
        }
    }
    
    stopifnot(is.list(ctrlFastqInput1))
    stopifnot(length(ctrlFastqInput1)>1)
    for(n in 1:length(ctrlFastqInput1)){
        ctrlFastqInput1[[n]] <- normalizePath(ctrlFastqInput1[[n]])
        checkFilePathExist(ctrlFastqInput1[[n]])
    }
    if(!is.null(ctrlAdapter1)){
        if(length(ctrlAdapter1)!=length(ctrlFastqInput1)){
            stop("the number of adapters in ctrlAdapter1 do not match the number of replicates in ctrlFastqInput1")
        }
    }
    
    stopifnot(is.list(ctrlFastqInput2))
    stopifnot(length(ctrlFastqInput2)>1)
    if(length(ctrlFastqInput1)!=length(ctrlFastqInput2)){
        stop(sprintf("ctrlFastqInput1 replicates size: %d and ctrlFastqInput2 replicates size %d do not match!",length(ctrlFastqInput1),length(ctrlFastqInput2)))
    }
    for(n in 1:length(ctrlFastqInput2)){
        if(length(ctrlFastqInput1[[n]])!=length(ctrlFastqInput2[[n]])){
            stop(sprintf("fastq file number does not match in ctrlFastqInput1 and ctrlFastqInput2 replicates %d",n))
        }
        ctrlFastqInput2[[n]] <- normalizePath(ctrlFastqInput2[[n]])
        checkFilePathExist(ctrlFastqInput2[[n]])
    }
    if(!is.null(ctrlAdapter1)){
        if(length(ctrlAdapter2)!=length(ctrlFastqInput2)){
            stop("the number of adapters in ctrlAdapter2 do not match the number of replicates in ctrlFastqInput2")
        }
    }
    
    
    
    param.tmp <- list(...)
    if(!(!is.null(param.tmp[["dontSet"]])&&param.tmp[["dontSet"]])){
        if(!is.null(refdir)){
            options(atacConf=setConfigure("refdir",refdir))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","refdir"))){
                dir.create(file.path("esATAC_pipeline","refdir"))
            }
            options(atacConf=setConfigure("refdir",file.path("esATAC_pipeline","refdir")))
        }
        if(!is.null(threads)){
            options(atacConf=setConfigure("threads",as.numeric(threads)))
            message(getConfigure("threads"))
        }else{
            options(atacConf=setConfigure("threads",2))
        }
        if(!is.null(genome)){
            options(atacConf=setConfigure("genome",genome))
        }else{
            stop("parameter genome is required")
        }
        if(!is.null(tmpdir)){
            options(atacConf=setConfigure("tmpdir",tmpdir))
        }else{
            if(!dir.exists("esATAC_pipeline")){
                dir.create("esATAC_pipeline")
            }
            if(!dir.exists(file.path("esATAC_pipeline","intermediate_results"))){
                dir.create(file.path("esATAC_pipeline","intermediate_results"))
            }else{
                warning(sprintf("path '%s' is exist",file.path("esATAC_pipeline","intermediate_results")))
            }
            options(atacConf=setConfigure("tmpdir",file.path("esATAC_pipeline","intermediate_results")))
            tmpdir <- "esATAC_pipeline"
        }
        if(is.null(param.tmp[["esATAC_result"]])){
            esATAC_result<-file.path(dirname(.obtainConfigure("tmpdir")),"esATAC_result")
            dir.create(esATAC_result)
        }else{
            esATAC_result <- param.tmp[["esATAC_result"]]
        }
        if(is.null(param.tmp[["esATAC_report"]])){
            esATAC_report<-file.path(dirname(.obtainConfigure("tmpdir")),"esATAC_report")
            dir.create(esATAC_report)
        }else{
            esATAC_report <- param.tmp[["esATAC_report"]]
        }
    }else{
        if(is.null(param.tmp[["esATAC_result"]])){
            esATAC_result<-NULL
        }else{
            esATAC_result <- param.tmp[["esATAC_result"]]
        }
        if(is.null(param.tmp[["esATAC_report"]])){
            esATAC_report<-NULL
        }else{
            esATAC_report <- param.tmp[["esATAC_report"]]
        }
    }
    
    if(is.null(use.SavedPWM)){
        if(is.null(motifPWM)){
            opts <- list()
            opts[["tax_group"]] <- "vertebrates"
            pwm <- getMatrixSet(JASPAR2016::JASPAR2016, opts)
            pwm <- TFBSTools::toPWM(pwm)
            names(pwm) <- TFBSTools::name(pwm)
            pwm <- lapply(X = pwm, FUN = TFBSTools::as.matrix)
            names(pwm) <- gsub(pattern = "[^a-zA-Z0-9]", replacement = "", x = names(pwm), perl = TRUE)
            pwm <- pwm[1:4]
        }else{
            pwm <- motifPWM
        }
    }
    
    message("Begin to process case sample...")
    tmpdir <- .obtainConfigure("tmpdir")
    dir.create(file.path(tmpdir,"case"))
    dir.create(file.path(esATAC_report,"case"))
    dir.create(file.path(esATAC_result,"case"))
    options(atacConf=setConfigure("tmpdir",file.path(tmpdir,"case")))
    
    
    if(is.null(use.SavedPWM)){
        caselist <- atacRepsPipe(genome = genome, fastqInput1 = caseFastqInput1,fastqInput2 = caseFastqInput2,
                                 adapter1 = caseAdapter1, adapter2 = caseAdapter2,interleave = interleave,
                                 createReport = TRUE, motifPWM =pwm, prefix = "Case_data", chr = chr, dontSet=TRUE, min.score = min.score,
                                 esATAC_report = file.path(esATAC_report,"case"), esATAC_result = file.path(esATAC_result,"case")) #saveTmp = TRUE,
    }else{
        caselist <- atacRepsPipe(genome = genome, fastqInput1 = caseFastqInput1,fastqInput2 = caseFastqInput2,
                                 adapter1 = caseAdapter1, adapter2 = caseAdapter2,interleave = interleave,
                                 createReport = TRUE, prefix = "Case_data", chr = chr, dontSet=TRUE, use.SavedPWM = use.SavedPWM,
                                 esATAC_report = file.path(esATAC_report,"case"), esATAC_result = file.path(esATAC_result,"case")) #saveTmp = TRUE,
    }
    
    message("Case sample process done")
    message(" ")
    message("Begin to process control sample")
    dir.create(file.path(tmpdir,"control"))
    dir.create(file.path(esATAC_report,"control"))
    dir.create(file.path(esATAC_result,"control"))
    options(atacConf=setConfigure("tmpdir",file.path(tmpdir,"control")))
    
    
    if(is.null(use.SavedPWM)){
        ctrllist <- atacRepsPipe(genome = genome, fastqInput1 = ctrlFastqInput1,fastqInput2 = ctrlFastqInput2,
                                 adapter1 = ctrlAdapter1, adapter2 = ctrlAdapter2,interleave = interleave,
                                 createReport = TRUE, motifPWM =pwm, prefix = "Control_data", chr = chr, dontSet=TRUE, min.score = min.score,
                                 esATAC_report = file.path(esATAC_report,"control"), esATAC_result = file.path(esATAC_result,"control")) #saveTmp = TRUE,
    }else{
        ctrllist <- atacRepsPipe(genome = genome, fastqInput1 = ctrlFastqInput1,fastqInput2 = ctrlFastqInput2,
                                 adapter1 = ctrlAdapter1, adapter2 = ctrlAdapter2,interleave = interleave,
                                 createReport = TRUE, prefix = "Control_data", chr = chr, dontSet=TRUE, use.SavedPWM = use.SavedPWM,
                                 esATAC_report = file.path(esATAC_report,"control"), esATAC_result = file.path(esATAC_result,"control")) #saveTmp = TRUE,
    }
    
    message("control sample process done")
    message(" ")
    message("Begin to generate case_control summary")
    
    dir.create(file.path(tmpdir,"case_control"))
    options(atacConf=setConfigure("tmpdir",file.path(tmpdir,"case_control")))
    
    bed.case <- getParam(caselist$atacProcs$sortedReadsBed, "bedOutput")
    bed.ctrl <- getParam(ctrllist$atacProcs$sortedReadsBed, "bedOutput")
    
    case.peak <- getParam(caselist$atacProcs$peakCalling, "bedOutput")
    ctrl.peak <- getParam(ctrllist$atacProcs$peakCalling, "bedOutput")
    
    CaseControlPeakCompPath <- paste(file.path(tmpdir,"case_control"), "/",
                                     c("Case_mergedReads_specific.bed", "Control_mergedReads_specific.bed", "overlap_mergedReads.bed"),
                                     sep = "")
    
    peakCom <- peakcomp(bedInput1 = case.peak, bedInput2 = ctrl.peak, bedOutput = CaseControlPeakCompPath)
    case_specific.peak <- getParam(peakCom, "bedOutput")[1]
    ctrl_specific.peak <- getParam(peakCom, "bedOutput")[2]
    overlap.peak <- getParam(peakCom, "bedOutput")[3]
    
    # for case
    Peakanno.case <- peakanno(peakInput = case_specific.peak)
    goAna.case <- atacGOAnalysis(atacProc = Peakanno.case, ont = "BP", pvalueCutoff = 0.01)
    
    # for ctrl
    Peakanno.ctrl <- peakanno(peakInput = ctrl_specific.peak)
    goAna.ctrl <- atacGOAnalysis(atacProc = Peakanno.ctrl, ont = "BP", pvalueCutoff = 0.01)
    
    if(is.null(use.SavedPWM)){
        print(case_specific.peak)
        print(ctrl_specific.peak)
        print(overlap.peak)
        mout <- atacMotifScanPair(atacProc = peakCom, motifPWM = pwm, min.score = min.score)
    }else{
        mout <- atacMotifScanPair(atacProc = peakCom, use.SavedPWM = use.SavedPWM)
    }
    
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
    
    
    
    conclusion <- list(caselist = caselist,
                       ctrllist = ctrllist,
                       comp_result = comp_result
    )
    
    message("Generate summary done")
    if(createReport){
        message("Begin to generate Report")
        #filename <- strsplit(case[["fastqInput1"]],".fastq|.FASTQ|.FQ|.fq")[[1]][1]
        #filename <- basename(filename)
        
        rmdfile<-system.file(package="esATAC", "extdata", "Rep_Report2.Rmd")
        rmdtext<-readChar(rmdfile,nchars=file.info(rmdfile)$size,useBytes = TRUE)
        #rmdtext<-sprintf(rmdtext,filename)
        
        workdir <- getwd()
        save(caselist,ctrllist,comp_result,workdir,file = file.path(.obtainConfigure("tmpdir"),"Report2.Rdata"))
        
        writeChar(rmdtext,con = file.path(.obtainConfigure("tmpdir"),"Rep_Report2.Rmd"),useBytes = TRUE)
        render(file.path(.obtainConfigure("tmpdir"),"Rep_Report2.Rmd"))
        #knit(file.path(.obtainConfigure("tmpdir"),"Report.Rmd"), file.path(.obtainConfigure("tmpdir"),"Report.md"))
        #markdownToHTML(file.path(.obtainConfigure("tmpdir"),"Report.md"), file.path(.obtainConfigure("tmpdir"),"Report.html"))
        #browseURL(paste0('file://', file.path(.obtainConfigure("tmpdir"),"Report.html")))
        message("Generate report done")
        if(!is.null(esATAC_report)&&!is.null(esATAC_result)){
            #generate report index
            replicateNum <- paste0("replicate ",1:length(caseFastqInput1))
            urllink <- paste0("<a href='./case/replicate",1:length(caseFastqInput1),"/Report.html'>","replicate ",1:length(caseFastqInput1)," report link</a>")
            caseSingleRep <- data.frame(Report_Name=replicateNum,Link=urllink)
            caseMergeConRep <- data.frame(Report_Name="concordance and merge",link="<a href='./case/rep_concord_merge/Rep_Report.html'>concordance and merge analysis report link</a>")
            
            replicateNum <- paste0("replicate ",1:length(ctrlFastqInput1))
            urllink <- paste0("<a href='./control/replicate",1:length(ctrlFastqInput1),"/Report.html'>","replicate ",1:length(ctrlFastqInput1)," report link</a>")
            ctrlSingleRep <- data.frame(Report_Name=replicateNum,Link=urllink)
            ctrlMergeConRep <- data.frame(Report_Name="concordance and merge",link="<a href='./control/rep_concord_merge/Rep_Report.html'>concordance and merge analysis report link</a>")
            
            case_control <- data.frame(Report_Name="case and control",link="<a href='./case_control/Rep_Report2.html'>case and control analysis report link</a>")
            save(caseSingleRep,caseMergeConRep,ctrlSingleRep,ctrlMergeConRep,case_control,file = file.path(.obtainConfigure("tmpdir"),"ReportIdx.Rdata"))
            
            rmdidxfile<-system.file(package="esATAC", "extdata", "Report_Index.Rmd")
            file.copy(rmdidxfile,.obtainConfigure("tmpdir"))
            render(file.path(.obtainConfigure("tmpdir"),"Report_Index.Rmd"))
            file.copy(file.path(.obtainConfigure("tmpdir"),"Report_Index.html"),file.path(esATAC_report,"Report.html"))
            
            # copy other files
            esATAC_report <- file.path(esATAC_report,"case_control")
            esATAC_result <- file.path(esATAC_result,"case_control")
            dir.create(esATAC_report)
            dir.create(esATAC_result)
            file.copy(file.path(.obtainConfigure("tmpdir"),"Rep_Report2.html"),esATAC_report, overwrite = TRUE)
            file.copy(getReportVal(caselist$atacProcs$goAna,"goOutput"),esATAC_report, overwrite = TRUE)
            file.copy(getReportVal(ctrllist$atacProcs$goAna,"goOutput"),esATAC_report, overwrite = TRUE)
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
            
            
            message(sprintf("type `browseURL(\"%s\")` to view Report in web browser",file.path(esATAC_report,"Rep_Report.html")))
        }else{
            message(sprintf("type `browseURL(\"%s\")` to view Report in web browser",file.path(.obtainConfigure("tmpdir"),"Rep_Report2.html")))
        }
    }
    
    invisible(conclusion)
    
}


