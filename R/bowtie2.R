#' @name bowtie2
#' @title Interface to bowtie2 of bowtie2-2.2.3
#' @description This function can be use to call \code{bowtie2}
#' that wrapped in shared library.
#' @param bt2Index \code{Character} scalar. bowtie2 index files
#' prefix: 'dir/basename'
#' (minus trailing '.*.bt2' of 'dir/basename.*.bt2').
#' @param samOutput \code{Character} scalar. A path to a SAM file
#' used for the alignment output.
#' @param seq1 \code{Character} vector. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates
#' paired with file paths in seq2.
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}
#' @param seq2 \code{Character} vector. It contains file paths with
#' #2 mates paired with file paths in seq1.
#' For single-end sequencing files and interleaved paired-end
#' sequencing files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' @param ... Additional arguments to be passed on to the binaries.
#' See below for details.
#' @param interleaved \code{Logical}. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param overwrite \code{Logical}. Force overwriting of existing
#' files if setting \code{TRUE}.
#' @details All additional arguments in ... are interpreted as
#' additional parameters to be passed on to
#' bowtie2_build. All of them should be \code{Character} or
#' \code{Numeric} scalar. You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8 --no-mixed")
#' with white space splited just like command line,
#' or put them in different \code{Character}
#' (e.g. "--threads","8","--no-mixed"). Note that some
#' arguments("-x","--interleaved","-U","-1","-2","-S") to the
#' bowtie2 are invalid if they are already handled as explicit
#' function arguments. See the output of
#' \code{bowtie2_usage()} for details about available parameters.
#' @author Zheng Wei
#' @return An invisible \code{Integer} of the shared library call
#' status. The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012).
#' Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export bowtie2
#' @examples
#' td <- tempdir()
#' ## Building a bowtie2 index
#' refs <- dir(system.file(package="Rbowtie2", "extdata", "bt2","refs"),
#' full=TRUE)
#' bowtie2_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads 4 --quiet",overwrite=TRUE)
#' ## Alignments
#' reads_1 <- system.file(package="Rbowtie2", "extdata", "bt2", "reads",
#' "reads_1.fastq")
#' reads_2 <- system.file(package="Rbowtie2", "extdata", "bt2", "reads",
#' "reads_2.fastq")
#' if(file.exists(file.path(td, "lambda_virus.1.bt2"))){
#'     cmdout<-bowtie2(bt2Index = file.path(td, "lambda_virus"),
#'        samOutput = file.path(td, "result.sam"),
#'        seq1=reads_1,seq2=reads_2,overwrite=TRUE,"--threads 3");cmdout
#'     head(readLines(file.path(td, "result.sam")))
#' }
#'

bowtie2 <- function(bt2Index,samOutput,seq1,...,seq2=NULL,interleaved=FALSE,
                    overwrite=FALSE){
    if(R.Version()$arch=="i386"){
        return("bowtie2 is not available for 32bit, please use 64bit R instead")
    }
    bt2Index <-trimws(as.character(bt2Index))
    samOutput<-trimws(as.character(samOutput))

    seq1<-trimws(as.character(seq1))


    if(!is.null(seq2)){
        seq2<-trimws(as.character(seq2))
        if(length(seq1)!=length(seq2)){
            stop(paste0("The lengths of arguments ",
                        "`seq1` and `seq2` should be the same length"))
        }
    }
    paramArray<-checkAddArgus("-x|--interleaved|-U|-1|-2|-S",...)


    if(interleaved){
        if(length(seq1)>1){
            stop(paste0("Argumnet `seq1` has to be a SINGLE file",
                        " path rather than a vector of paths"))
        }
        if(!is.null(seq2)){
            stop("Argumnet `seq2` has to be NULL when interleaved=TRUE")
        }
    }


    checkFileExist(seq1,"seq1")
    checkFileExist(seq2,"seq2")
    checkPathExist(bt2Index,"bt2Index")
    checkFileExist(paste0(bt2Index,".1.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".2.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".3.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".4.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".rev.1.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".rev.2.bt2"),"bt2Index")
    checkFileCreatable(samOutput,"samOutput",overwrite)

    argvs = c("-x",bt2Index)
    seq1<-paste0(seq1,collapse = ",")
    if(is.null(seq2)){
        if(interleaved){
            argvs <- c(argvs,"--interleaved",seq1)
        }else{
            argvs <- c(argvs,"-U",seq1)
        }
    }else{
        seq2<-paste0(seq2,collapse = ",")
        argvs <- c(argvs,"-1",seq1,"-2",seq2)
    }

    argvs <- c(paramArray,argvs,"-S",samOutput)

    invisible(.callbinary("bowtie2-align-s",paste(argvs,collapse = " ")))

}


#' @name bowtie2-build
#' @title Interface to bowtie2-build of bowtie2-2.2.3
#' @description This function can be use to call \code{bowtie2-build}
#' that wrapped in shared library.
#' @param references \code{Character} vector. The path to the files containing
#' the references for which to
#' build a bowtie index.
#' @param bt2Index \code{Character} scalar. Write bowtie2 index data to files
#' with this prefix: 'dir/basename'.
#' If the files with path like 'dir/basename.*.bt2' already exists,
#' the function function will cast an error,
#' unless argument overwrite is \code{TRUE}.
#' @param ... Additional arguments to be passed on to the binaries.
#' See below for details.
#' @param overwrite \code{Logical}. Force overwriting of existing files
#' if setting \code{TRUE}.
#' @details All additional arguments in ... are interpreted as additional
#' parameters to be passed on to
#' bowtie2_build. All of them should be \code{Character} or
#' \code{Numeric} scalar. You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8 --quiet") with white
#' space splited just like command line,
#' or put them in different \code{Character}(e.g. "--threads","8","--quiet").
#' See the output of
#' \code{bowtie2_build_usage()} for details about available parameters.
#' @author Zheng Wei
#' @return An invisible \code{Integer} of the shared library call status.
#' The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read
#' alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export bowtie2_build
#' @examples
#' td <- tempdir()
#' ## Building a bowtie2 index
#' refs <- dir(system.file(package="Rbowtie2", "extdata", "bt2","refs"),
#' full=TRUE)
#' cmdout<-bowtie2_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads 4 --quiet",overwrite=TRUE);cmdout
#' ## Use additional arguments in another way
#' cmdout<-bowtie2_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' "--threads",4,"--quiet",overwrite=TRUE);cmdout
#' ## The function will print the output
#' ## during the process without "--quiet" argument.
#' cmdout<-bowtie2_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
#' overwrite=TRUE);cmdout

bowtie2_build <- function(references,bt2Index,...,overwrite=FALSE){
    if(R.Version()$arch=="i386"){
        return("bowtie2 is not available for 32bit, please use 64bit R instead")
    }
    references<- trimws(as.character(references))
    bt2Index <- trimws(as.character(bt2Index))

    paramArray<-checkAddArgus("noinvalid",...)

    checkFileExist(references,"references")
    checkPathExist(bt2Index,"bt2Index")
    checkFileCreatable(paste0(bt2Index,".1.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".2.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".3.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".4.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".rev.1.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".rev.2.bt2"),"bt2Index",overwrite)

    references<-paste0(references,collapse = ",")
    argvs <- c(paramArray,references,bt2Index)


    invisible(.callbinary("bowtie2-build-s",paste(argvs,collapse = " ")))

}

#' @name bowtie2_version
#' @title Print version information of bowtie2-2.2.3
#' @description Print version information of bowtie2-2.2.3
#' @author Zheng Wei
#' @return An invisible \code{Integer} of the shared library call status.
#' The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read
#' alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export bowtie2_version
#' @examples
#' cmdout<-bowtie2_version();cmdout
bowtie2_version <- function(){
    if(R.Version()$arch=="i386"){
        return("bowtie2 is not available for 32bit, please use 64bit R instead")
    }
    invisible(.callbinary("bowtie2-align-s","--version"))
}

#' @name bowtie2_usage
#' @title Print available arguments for bowtie2
#' @description Note that some arguments to the
#' bowtie2 are invalid if they are
#' already handled as explicit function arguments.
#' @author Zheng Wei
#' @return An invisible \code{Integer} of the shared library call status.
#' The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read
#' alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export bowtie2_usage
#' @examples
#' cmdout<-bowtie2_usage();cmdout
bowtie2_usage <- function(){
    if(R.Version()$arch=="i386"){
        return("bowtie2 is not available for 32bit, please use 64bit R instead")
    }
    invisible(.callbinary("bowtie2-align-s","-h"))
}

#' @name bowtie2_build_usage
#' @title Print available arguments for bowtie2_build_usage
#' @description Note that some arguments to the
#' bowtie2_build_usage are invalid if they are
#' already handled as explicit function arguments.
#' @author Zheng Wei
#' @return An invisible \code{Integer} of the shared library call status.
#' The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead B, Salzberg S.
#' Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
#' @export bowtie2_build_usage
#' @examples
#' cmdout<-bowtie2_build_usage();cmdout
bowtie2_build_usage <- function() {
    if(R.Version()$arch=="i386"){
        return("bowtie2 is not available for 32bit, please use 64bit R instead")
    }
    invisible(.callbinary("bowtie2-build-s","-h"))
}


