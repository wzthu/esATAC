#' @name identify_adapters
#' @title identify adapters for paired-end reads
#' @description This function can be use to call \code{AdapterRemoval}
#' that wrapped in shared library for adapters identifying.
#' @param file1 \code{Character} vector. It can be file paths with #1
#' mates paired with file paths in file2
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}
#' @param file2 \code{Character} vector. It contains file paths with
#' #2 mates paired with file paths in seq1.
#' For interleaved paired-end sequencing files(argument
#' interleaved=\code{TRUE}),it must to be setted to \code{NULL}.
#' @param ... Additional arguments to be passed on to the binaries.
#' See below for details.
#' @param basename \code{Character}. The outputfile path prefix.
#' Default: your_output
#' @param overwrite \code{Logical}. Force overwriting of existing
#' files if setting \code{TRUE}.
#' @details All additional arguments in ... are interpreted as
#' additional parameters to be passed on to
#' bowtie2_build. All of them should be \code{Character} or
#' \code{Numeric} scalar. You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8") with white
#'  space splited just like command line,
#' or put them in different \code{Character}(e.g. "--threads","8").
#' Note that some arguments("--identify-adapters",
#' "--file1","--file2","--basename") to the
#' identify_adapters are invalid if they are already handled as
#' explicit function arguments. See the output of
#' \code{adapterremoval_usage()} for details about available parameters.
#' @return An invisible \code{Character} vector of adapters for each mate.
#' @author Zheng Wei
#' @references    Schubert, Lindgreen, and Orlando (2016).
#' AdapterRemoval v2: rapid adapter trimming, identification, and read merging.
#' BMC Research Notes, 12;9(1):88.
#' @export identify_adapters
#' @examples
#' td <- tempdir()
#' reads_1 <- system.file(package="ATACFlow", "extdata", "adrm", "reads_1.fq")
#' reads_2 <- system.file(package="ATACFlow", "extdata", "adrm", "reads_2.fq")
#' adapters <- identify_adapters(file1=reads_1,file2=reads_2,
#' basename = file.path(td,"reads")
#' ,"--threads 2",overwrite=TRUE)
#' adapters
identify_adapters <- function(file1,file2,...,basename = NULL,
    overwrite = FALSE){
    file1<-trimws(as.character(file1))
    if(!is.null(file2)){
        file2<-trimws(as.character(file2))
    }
    if(!is.null(basename)){
        basename<-trimws(as.character(basename))
    }
    checkFileExist(file1,"file1")
    checkFileExist(file2,"file2")
    checkFileCreatable(paste0(basename,".adapter1"),"file1",overwrite)
    checkFileCreatable(paste0(basename,".adapter2"),"file2",overwrite)

    paramArray<-
        checkAddArgus("--identify-adapters|--file1|--file2|--basename",...)
    if(is.null(file2)){
        argvs<-c("--identify-adapters","--file1",
                 file1,"--interleaved",paramArray);
    }else{
        argvs<-c("--identify-adapters","--file1",
                 file1,"--file2",file2,paramArray);
    }
    if(!is.null(basename)){
        argvs<-c(argvs,"--basename",basename)
    }

    .callbinary("AdapterRemoval",paste(argvs,collapse = " "))
    if(is.null(basename)){
        basename<-"your_output"
    }
    adapter1tb<-readLines(paste0(basename,".adapter1"));
    adapter2tb<-readLines(paste0(basename,".adapter2"));
    return(c(adapter1tb,adapter2tb))
}


#' @name remove_adapters
#' @title Interface to bowtie2 of adapterremoval-2.2.1a
#' @description This function can be use to call \code{AdapterRemoval}
#' that wrapped in shared library.
#' @param file1 \code{Character} vector. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates paired
#' with file paths in file2
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}
#' @param adapter1 \code{Character}. It is an adapter sequence for file1.
#' Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
#' @param file2 \code{Character} vector. It contains file paths with #2
#' mates paired with file paths in file1.
#' For single-end sequencing files and interleaved paired-end sequencing
#' files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' @param output1 \code{Character}. The trimmed mate1 reads output file
#' path for file1. Defualt:
#' basename.pair1.truncated (paired-end),
#' basename.truncated (single-end), or
#' basename.paired.truncated (interleaved)
#' @param adapter2 \code{Character}. It is an adapter sequence for file2.
#' Defualt: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#' @param output2 \code{Character}. The trimmed mate2 reads output file
#' path for file2. Default:
#' BASENAME.pair2.truncated (only used in PE mode, but not if
#' --interleaved-output is enabled)
#' @param ... Additional arguments to be passed on to the binaries.
#' See below for details.
#' @param basename \code{Character}. The outputfile path prefix.
#' Default: your_output
#' @param interleaved \code{Logical}. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param overwrite \code{Logical}. Force overwriting of existing files
#' if setting \code{TRUE}.
#' @details All additional arguments in ... are interpreted as additional
#' parameters to be passed on to
#' bowtie2_build. All of them should be \code{Character} or \code{Numeric}
#' scalar. You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8") with white space
#' splited just like command line,
#' or put them in different \code{Character}(e.g. "--threads","8").
#' Note that some arguments(
#' "--file1","--file2","--adapter1","--adapter2","--output1","--output2",
#' "--basename","--interleaved") to the
#' identify_adapters are invalid if they are already handled as explicit
#' function arguments. See the output of
#' \code{adapterremoval_usage()} for details about available parameters.
#' @author Zheng Wei
#' @return An invisible \code{Integer} of the shared library call status.
#' The value is 0 when there is not any mistake.
#' Otherwise the value is non-zero.
#' @references    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval
#' v2: rapid
#' adapter trimming, identification, and read merging.
#' BMC Research Notes, 12;9(1):88.
#' @export remove_adapters
#' @examples
#' td <- tempdir()
#'
#' # Identify adapters
#' reads_1 <- system.file(package="ATACFlow", "extdata", "adrm", "reads_1.fq")
#' reads_2 <- system.file(package="ATACFlow", "extdata", "adrm", "reads_2.fq")
#' adapters <- identify_adapters(file1=reads_1,file2=reads_2,
#' basename=file.path(td,"reads"), "--threads 3",overwrite=TRUE)
#'
#' # Remove adapters
#' cmdout<-remove_adapters(file1=reads_1,file2=reads_2,adapter1 = adapters[1],
#' adapter2 = adapters[2],
#' output1=file.path(td,"reads_1.trimmed.fq"),
#' output2=file.path(td,"reads_2.trimmed.fq"),
#' basename=file.path(td,"reads.base"),overwrite=TRUE,"--threads 3");cmdout
#'
remove_adapters <- function(file1,...,adapter1 = NULL,output1 = NULL,
                            file2 = NULL,adapter2 = NULL,output2 = NULL,
                            basename = NULL,interleaved = FALSE,
                            overwrite = FALSE){
    file1<-trimws(as.character(file1))
    if(!is.null(adapter1)){
        adapter1<-trimws(as.character(adapter1))
    }
    if(!is.null(output1)){
        output1<-trimws(as.character(output1))
    }
    if(!is.null(file2)){
        if(interleaved){
            stop("Argumnet `seq2` has to be NULL when interleaved=TRUE")
        }else{
            file2<-trimws(as.character(file2))
            if(length(file1)!=length(file2)){
                stop(paste0("The lengths of arguments `file1` ",
                            "and `file2` should be the same length"))
            }
        }
    }
    if(!is.null(adapter2)){
        adapter2<-trimws(as.character(adapter2))
    }
    if(!is.null(output2)){
        output2<-trimws(as.character(output2))
    }
    if(!is.null(basename)){
        basename<-trimws(as.character(basename))
    }

    paramArray<-
        checkAddArgus(paste0("--file1|--file2|--adapter1|--adapter2|",
                             "--output1|--output2|--basename|--interleaved"),...)
    checkFileExist(file1,"file1")
    checkFileExist(file2,"file2")
    checkFileCreatable(output1,"output1",overwrite)
    checkFileCreatable(output2,"output2",overwrite)
    checkPathExist(basename,"basename")

    argvs<-c("--file1",file1)

    if(!is.null(adapter1)){
        argvs<-c(argvs,"--adapter1",adapter1)
    }
    if(!is.null(output1)){
        argvs<-c(argvs,"--output1",output1)
    }
    if(!is.null(file2)){
        argvs<-c(argvs,"--file2",file2)
    }
    if(!is.null(adapter2)){
        argvs<-c(argvs,"--adapter2",adapter2)
    }
    if(!is.null(output2)){
        argvs<-c(argvs,"--output2",output2)
    }
    if(!is.null(basename)){
        argvs<-c(argvs,"--basename",basename)
    }
    if(interleaved){
        argvs<-c(argvs,"--interleaved")
    }
    argvs<-c(argvs,paramArray)
    invisible(.callbinary("AdapterRemoval",paste(argvs,collapse = " ")))


}


#' @name adapterremoval_usage
#' @title Print available arguments for adapterremoval
#' @description Print available arguments for adapterremoval.
#' Note that some arguments to the
#' adapterremoval are invalid if they are
#' already handled as explicit function arguments.
#' @return An invisible \code{Integer} of the shared library call status.
#' The value is 0 when there is not any mistakes.
#' Otherwise the value is non-zero.
#' @author Zheng Wei
#' @references Schubert, Lindgreen, and Orlando (2016).
#' AdapterRemoval v2: rapid adapter trimming, identification, and read merging.
#' BMC Research Notes, 12;9(1):88.
#' @export adapterremoval_usage
#' @examples
#' adapterremoval_usage()
adapterremoval_usage<- function(){
    invisible(.callbinary("AdapterRemoval","-h"))
}


#' @name adapterremoval_version
#' @title Print version information of adapterremoval
#' @description Print version information of adapterremoval
#' @return An invisible \code{Integer} of the shared library call status.
#' The value is 0 when there is not any mistakes
#' @author  Zheng Wei
#' @references Schubert, Lindgreen, and Orlando (2016).
#' AdapterRemoval v2: rapid adapter trimming, identification, and read merging.
#' BMC Research Notes, 12;9(1):88.
#' @export adapterremoval_version
#' @examples
#' adapterremoval_version()
adapterremoval_version<- function(){
    invisible(.callbinary("AdapterRemoval","--version"))
}

.callbinary<- function(bin, args)
{
    args <- gsub("^ *| *$", "", args)
    call <- paste(shQuote(file.path(system.file(package="ATACFlow"), bin)), args)
    output <- system(call, intern=TRUE,show.output.on.console=TRUE)
    return(output)
}
