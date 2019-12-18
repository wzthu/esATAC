setClass(Class = "Bowtie2Mapping",
         contains = "ATACProc"
)



setMethod(
    f = "init",
    signature = "Bowtie2Mapping",
    definition = function(.Object, prevSteps = list(),...){
        allparam <- list(...)
        fastqInput1 <- allparam[["fastqInput1"]]
        fastqInput2 <- allparam[["fastqInput2"]]
        samOutput <- allparam[["samOutput"]]
        bt2Idx <- allparam[["bt2Idx"]]
        reportOutput <- allparam[["reportOutput"]]
        interleave <- allparam[["interleave"]]
        threads <- allparam[["threads"]]
        paramList <- allparam[["paramList"]]
        
        if(length(prevSteps) > 0){
            fastqStep <- prevSteps[[1]]
            input(.Object)[["fastqInput1"]] <- output(fastqStep)[["fastqOutput1"]]
            input(.Object)[["fastqInput2"]] <- output(fastqStep)[["fastqOutput2"]]
            param(.Object)[["interleave"]] <- property(fastqStep)[["interleave"]]
            param(.Object)[["singleEnd"]] <- property(fastqStep)[["singleEnd"]]
        }else{
            param(.Object)[["interleave"]] <- interleave
            property(.Object)[["interleave"]] <- interleave
            if(interleave){
                if(!is.null(fastqInput2)){
                    stop("interleave data should put in one fastq file")
                }else{
                    property(.Object)$singleEnd <- TRUE
                    param(.Object)$singleEnd <- TRUE
                }
            }else{
                property(.Object)$singleEnd <- is.null(fastqInput2)
                param(.Object)$singleEnd<-is.null(fastqInput2)
            }
        }

        if(!is.null(fastqInput1)){
            input(.Object)[["fastqInput1"]] <- fastqInput1;
        }
        if(!is.null(fastqInput2)){
            input(.Object)[["fastqInput2"]] <- fastqInput2;
        }


        if(is.null(samOutput)){
            if(!is.null(input(.Object)[["fastqInput1"]])){
                output(.Object)[["samOutput"]] <- 
                    getAutoPath(.Object,input(.Object)[["fastqInput1"]],"fastq|fq","sam")
            }
        }else{
            output(.Object)[["samOutput"]] <- samOutput;
        }

        if(is.null(bt2Idx)){
            param(.Object)[["bt2Idx"]] <- file.path(getRefDir(),getRefRc("bt2Idx"))
            input(.Object)[["bt2Idxs"]] <- getRefFiles("bt2Idx")
        }else{
            param(.Object)[["bt2Idx"]] <- bt2Idx
            input(.Object)[["bt2Idxs"]] <- paste0(bt2Idx,c(".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"))
        }


        if(!is.null(paramList)){
            paramList<-trimws(as.character(paramList))
            paramList<-paste(paramList,collapse = " ")
            paramList <- strsplit(paramList,"\\s+")[[1]]
            if(length(paramList)>0){
                rejectp<-"-p|--threads|-x|-1|-2|-U|-S|--interleaved"
                checkParam(paramList,rejectp)
                param(.Object)[["paramList"]]<-paramList
            }
        }

        if(is.null(reportOutput)){
            if(!is.null(input(.Object)[["fastqInput1"]])){
                output(.Object)[["reportOutput"]] <- getAutoPath(.Object, input(.Object)[["fastqInput1"]], "fq|fastq", "report")
            }
        }else{
            output(.Object)[["reportOutput"]] <- reportOutput;
        }
        
        stopifnot(is.numeric(threads))
        param(.Object)[["threads"]] <- threads
        
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "Bowtie2Mapping",
    definition = function(.Object,...){
        paramList <- NULL
        if(!is.null(param(.Object)[["paramList"]])){
            paramList <- paste(param(.Object)[["paramList"]],collapse = " ")
        }

        paramList<-paste("-p",as.character(param(.Object)[["threads"]]),paramList)
 

        writeLog(.Object,"start mapping with parameters:")
        writeLog(.Object,paste0("bowtie2 index:",param(.Object)[["bt2Idx"]]))
        writeLog(.Object,paste0("samOutput:",output(.Object)[["samOutput"]]))
        writeLog(.Object,paste0("report:",output(.Object)[["reportOutput"]]))
        writeLog(.Object,paste0("fastqInput1:",input(.Object)[["fastqInput1"]]))
        writeLog(.Object,paste0("fastqInput2:",input(.Object)[["fastqInput2"]]))
        writeLog(.Object,paste0("other parameters:",paramList))
        if(length(paramList)>0){
            rs<-bowtie2(bt2Index = param(.Object)[["bt2Idx"]],
                        samOutput = output(.Object)[["samOutput"]],
                        seq1 = input(.Object)[["fastqInput1"]],
                        paramList,
                        seq2 = input(.Object)[["fastqInput2"]],
                        interleaved = param(.Object)[["interleave"]],
                        overwrite=TRUE)
        }else{
            rs<-bowtie2(bt2Index = param(.Object)[["bt2Idx"]],
                        samOutput = output(.Object)[["samOutput"]],
                        seq1 = input(.Object)[["fastqInput1"]],
                        seq2 = input(.Object)[["fastqInput2"]],
                        interleaved = param(.Object)[["interleave"]],
                        overwrite=TRUE)
        }

        writeLines(text = rs,con = output(.Object)[["reportOutput"]])
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "Bowtie2Mapping",
    definition = function(.Object, ...){
        txt <- readLines(output(.Object)[["reportOutput"]])
        s<-strsplit(txt[1]," ")
        report(.Object)$total <- as.integer(s[[1]][1])
        
        s<-strsplit(txt[length(txt)],"% ")
        report(.Object)$maprate <- (as.numeric(s[[1]][1])/100)
        
        report(.Object)$detail <- txt
        .Object
    }
)




#' @name Bowtie2Mapping
#' @importFrom Rbowtie2 bowtie2
#' @title Use bowtie2 aligner to map reads to reference genome
#' @description
#' Use bowtie2 aligner to map reads to reference genome
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacRemoveAdapter}}
#' \code{\link{removeAdapter}}
#' @param reportOutput \code{Character} scalar.
#' The prefix of report files path.
#' @param bt2Idx \code{Character} scalar.
#' bowtie2 index files
#' prefix: 'dir/basename'
#' (minus trailing '.*.bt2' of 'dir/basename.*.bt2').
#' @param samOutput \code{Character} scalar.
#' A path to a SAM file
#' used for the alignment output.
#' @param fastqInput1 \code{Character} vector. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates
#' paired with file paths in fastqInput2.
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}
#' @param fastqInput2 \code{Character} vector. It contains file paths with
#' #2 mates paired with file paths in fastqInput1.
#' For single-end sequencing files and interleaved paired-end
#' sequencing files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' @param paramList Additional arguments to be passed on to the binaries.
#' See below for details.
#' @param interleave \code{Logical}. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param threads \code{Integer} scalar.
#' The threads will be created in this process. default: getThreads()
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{bowtie2Mapping} instead.
#’ All additional arguments in paramList are interpreted as
#' additional parameters to be passed on to
#' bowtie2. You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8 --no-mixed")
#' with white space splited just like command line,
#' or put them as \code{Character} vector
#' (e.g. c("--threads","8","--no-mixed")). Note that some
#' arguments("-x","--interleaved","-U","-1","-2","-S","threads") to the
#' bowtie2 are invalid if they are already handled as explicit
#' function arguments. See the output of
#' \code{bowtie2_usage()} for details about available parameters.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacRemoveAdapter}}
#' \code{\link{removeAdapter}}
#' \code{\link{bowtie2}}
#' \code{\link{bowtie2_build}}
#' \code{\link{bowtie2_usage}}
#' \code{\link{atacSam2Bam}}
#' \code{\link{atacSamToBed}}
#' \code{\link{atacLibComplexQC}}
#' @examples
#' td <- tempdir()
#' setTmpDir(td)
#'
#' ## Building a bowtie2 index
#' library("Rbowtie2")
#' refs <- dir(system.file(package="esATAC", "extdata", "bt2","refs"),
#' full=TRUE)
#' bowtie2_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads 4 --quiet",overwrite=TRUE)
#' ## Alignments
#' reads_1 <- system.file(package="esATAC", "extdata", "bt2", "reads",
#' "reads_1.fastq")
#' reads_2 <- system.file(package="esATAC", "extdata", "bt2", "reads",
#' "reads_2.fastq")
#' if(file.exists(file.path(td, "lambda_virus.1.bt2"))){
#'     (bowtie2Mapping(bt2Idx = file.path(td, "lambda_virus"),
#'        samOutput = file.path(td, "result.sam"),
#'        fastqInput1=reads_1,fastqInput2=reads_2,threads=3))
#'     head(readLines(file.path(td, "result.sam")))
#' }



setGeneric("atacBowtie2Mapping",function(atacProc,samOutput=NULL,
                                         reportOutput =NULL, bt2Idx=NULL,
                                         fastqInput1=NULL, fastqInput2=NULL, 
                                         interleave = FALSE, threads = getThreads(), 
                                         paramList="--no-discordant --no-unal --no-mixed -X 2000", ...) 
                                         standardGeneric("atacBowtie2Mapping"))



#' @rdname Bowtie2Mapping
#' @aliases atacBowtie2Mapping
#' @export
setMethod(
    f = "atacBowtie2Mapping",
    signature = "ATACProc",
    definition = function(atacProc,samOutput=NULL,
                          reportOutput =NULL, bt2Idx=NULL,
                          fastqInput1=NULL, fastqInput2=NULL, 
                          interleave = FALSE, threads = getThreads(), 
                          paramList="--no-discordant --no-unal --no-mixed -X 2000", ...){
        allpara <- c(list(Class = "Bowtie2Mapping", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname Bowtie2Mapping
#' @aliases bowtie2Mapping
#' @export
bowtie2Mapping <- function(fastqInput1, fastqInput2=NULL,
                           samOutput=NULL,reportOutput =NULL, 
                           bt2Idx=NULL, interleave = FALSE, 
                           threads = getThreads(), 
                           paramList="--no-discordant --no-unal --no-mixed -X 2000", ...){
    allpara <- c(list(Class = "Bowtie2Mapping", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
