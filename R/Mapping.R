Bowtie2Mapping <-R6Class(
    classname = "Bowtie2Mapping",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc,samOutput=NULL, bt2Idx=NULL, 
                              fastqInput1=NULL, fastqInput2=NULL, 
                              interleave = FALSE, 
                              paramList="--no-discordant --no-unal --no-mixed -X 2000",
                              reportOutput = NULL, threads = NULL, editable=FALSE){
            super$initialize("Bowtie2Mapping",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1");
                private$paramlist[["fastqInput2"]] <- atacProc$getParam("fastqOutput2");
                regexProcName<-sprintf("(fastq|fq|%s)",atacProc$getProcName())
                private$paramlist[["interleave"]] <- atacProc$getParam("interleave")
            }else{
                regexProcName<-"(fastq|fq)"
                private$paramlist[["interleave"]] <- interleave
                if(is.null(fastqInput2)){
                    private$singleEnd<-TRUE
                }else{
                    private$singleEnd<-FALSE
                }
            }

            if(!is.null(fastqInput1)){
                private$paramlist[["fastqInput1"]] <- fastqInput1;
            }
            if(!is.null(fastqInput2)){
                private$paramlist[["fastqInput2"]] <- fastqInput2;
            }


            if(is.null(samOutput)){
                if(!is.null(private$paramlist[["fastqInput1"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput1"]],regexProcName)
                    private$paramlist[["samOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".sam"))
                }
            }else{
                private$paramlist[["samOutput"]] <- samOutput;
            }

            if(is.null(bt2Idx)){
                private$paramlist[["bt2Idx"]]=.obtainConfigure("bt2Idx")
            }else{
                private$paramlist[["bt2Idx"]]<-bt2Idx
            }

            
            if(!is.null(paramList)){
                paramList<-trimws(as.character(paramList))
                paramList<-paste(paramList,collapse = " ")
                paramList <- strsplit(paramList,"\\s+")[[1]]
                if(length(paramList)>0){
                    rejectp<-"-p|--threads|-x|-1|-2|-U|-S|--interleaved"
                    private$checkParam(paramList,rejectp)
                    private$paramlist[["paramList"]]<-paramList
                }
            }

            if(is.null(reportOutput)){
                if(!is.null(private$paramlist[["fastqInput1"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput1"]],regexProcName)
                    private$paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
            }else{
                private$paramlist[["reportOutput"]] <- reportOutput;
            }
            if(!is.null(threads)){
                private$paramlist[["threads"]] <- as.integer(threads)
            }
            private$paramValidation()
        }


    ),
    private = list(
        processing = function(){
            if(!is.null(private$paramlist[["threads"]])){
                if(private$paramlist[["threads"]]>1){
                    paramList<-paste(c("-p",as.character(private$paramlist[["threads"]]),private$paramlist[["paramList"]]),collapse = " ")
                }
            }else if(.obtainConfigure("threads")>1){
                paramList<-paste(c("-p",as.character(.obtainConfigure("threads")),private$paramlist[["paramList"]]),collapse = " ")
            }
            
            private$writeLog("start mapping with parameters:")
            private$writeLog(paste0("bowtie2 index:",private$paramlist[["bt2Idx"]]))
            private$writeLog(paste0("samOutput:",private$paramlist[["samOutput"]]))
            private$writeLog(paste0("report:",private$paramlist[["reportOutput"]]))
            private$writeLog(paste0("fastqInput1:",private$paramlist[["fastqInput1"]]))
            private$writeLog(paste0("fastqInput2:",private$paramlist[["fastqInput2"]]))
            private$writeLog(paste0("other parameters:",paste(private$paramlist[["paramList"]],collapse = " ")))
            if(length(paramList>0)){
                rs<-bowtie2(bt2Index = private$paramlist[["bt2Idx"]],
                        samOutput = private$paramlist[["samOutput"]],
                        seq1 = private$paramlist[["fastqInput1"]],
                        paramList,
                        seq2 = private$paramlist[["fastqInput2"]],
                        interleaved = private$paramlist[["interleave"]],
                        overwrite=TRUE)
            }else{
                rs<-bowtie2(bt2Index = private$paramlist[["bt2Idx"]],
                        samOutput = private$paramlist[["samOutput"]],
                        seq1 = private$paramlist[["fastqInput1"]],
                        seq2 = private$paramlist[["fastqInput2"]],
                        interleaved = private$paramlist[["interleave"]],
                        overwrite=TRUE)
            }
            
            writeLines(text = rs,con = private$paramlist[["reportOutput"]])



        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["fastqInput1"]])){
                stop("fastqInput1 is required.")
            }
            if(is.null(private$paramlist[["bt2Idx"]])){
                stop("bt2Idx is required")
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["fastqInput1"]]);
            private$checkFileExist(private$paramlist[["fastqInput2"]]);
            private$checkFileCreatable(private$paramlist[["samOutput"]]);
            private$checkFileCreatable(private$paramlist[["reportOutput"]]);
        },
        getReportValImp = function(item){
            txt <- readLines(private$paramlist[["reportOutput"]])
            if(item == "total"){
                s<-strsplit(txt[1]," ")
                return(as.integer(s[[1]][1]))
            }
            if(item == "maprate"){
                s<-strsplit(txt[length(txt)],"% ") 
                return(as.numeric(s[[1]][1])/100)
            }
            if(item == "detail"){
                return(txt)
            }
            stop(paste0(item," is not an item of report value."))
        },
        getReportItemsImp = function(){
            return(c("total","maprate","detail"))
        }
    )


)

#' @name atacBowtie2Mapping
#' @aliases atacBowtie2Mapping
#' @aliases bowtie2Mapping
#' @importFrom Rbowtie2 bowtie2
#' @title Use bowtie2 aligner to map reads to reference genome 
#' @description 
#' Use bowtie2 aligner to map reads to reference genome 
#' @param atacProc \code{\link{ATACProc}} object scalar. 
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
#' The threads will be created in this process. default: 1
#' @details The parameter related to input and output file path
#' will be automatically 
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or 
#' generated based on known parameters 
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently, 
#' \code{atacProc} should be set \code{NULL} 
#' or you can use \code{fregLenDistr} instead.
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
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
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
#' setConfigure("tmpdir",td)
#' 
#' ## Building a bowtie2 index
#' library("Rbowtie2")
#' refs <- dir(system.file(package="ATACFlow", "extdata", "bt2","refs"),
#' full=TRUE)
#' bowtie2_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads 4 --quiet",overwrite=TRUE)
#' ## Alignments
#' reads_1 <- system.file(package="ATACFlow", "extdata", "bt2", "reads",
#' "reads_1.fastq")
#' reads_2 <- system.file(package="ATACFlow", "extdata", "bt2", "reads",
#' "reads_2.fastq")
#' if(file.exists(file.path(td, "lambda_virus.1.bt2"))){
#'     (atacBowtie2Mapping(NULL,bt2Idx = file.path(td, "lambda_virus"),
#'        samOutput = file.path(td, "result.sam"),
#'        fastqInput1=reads_1,fastqInput2=reads_2,threads=3))
#'     head(readLines(file.path(td, "result.sam")))
#' }

#' @rdname atacBowtie2Mapping
#' @export 
atacBowtie2Mapping <- function(atacProc,samOutput=NULL,reportOutput =NULL, bt2Idx=NULL,fastqInput1=NULL, fastqInput2=NULL, interleave = FALSE, threads = NULL, paramList="--no-discordant --no-unal --no-mixed -X 2000"){
    atacproc<-Bowtie2Mapping$new(atacProc=atacProc,bt2Idx=bt2Idx,samOutput=samOutput, fastqInput1=fastqInput1,
                                 fastqInput2=fastqInput2, interleave = interleave, paramList=paramList,reportOutput=reportOutput, threads= threads)
    atacproc$process()
    invisible(atacproc)
}

#' @rdname atacBowtie2Mapping
#' @export
bowtie2Mapping <- function(fastqInput1, fastqInput2=NULL,samOutput=NULL,reportOutput =NULL, bt2Idx=NULL, interleave = FALSE, threads = NULL, paramList="--no-discordant --no-unal --no-mixed -X 2000"){
    atacproc<-Bowtie2Mapping$new(atacProc=NULL,bt2Idx=bt2Idx,samOutput=samOutput, fastqInput1=fastqInput1,
                                 fastqInput2=fastqInput2, interleave = interleave, paramList=paramList,reportOutput=reportOutput, threads= threads)
    atacproc$process()
    invisible(atacproc)
}
