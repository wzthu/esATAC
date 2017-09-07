.renamer_call<-function(inputFile,outputFile,fileType="fq",interleave = FALSE){
  argv<-list(inputFile=inputFile,outputFile=outputFile,fileType=fileType,interleave=interleave);
  return(renamer(argv));
}

.sam2bed_merge_call <- function(samfile, bedfile,posOffset,negOffset,sortBed,uniqueBed,filterList,minFregLen,maxFregLen,saveExtLen,downSample = 2e9)
{
    argv <- list(samfile = samfile, bedfile = bedfile ,posOffset = posOffset,negOffset = negOffset,
                 sort = sortBed,unique = uniqueBed, minFregLen = minFregLen, maxFregLen = maxFregLen, saveExtLen = saveExtLen, memSize = 8 ,downSample = downSample,removeXS = TRUE)
    print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    return(R_sam2bed_merge_wrapper(argv,filterList))
}

.sam2bed_call <- function(samfile, bedfile,posOffset,negOffset,sortBed,uniqueBed,filterList,downSample = 2e9){
    argv <- list(samfile = samfile, bedfile = bedfile ,posOffset = posOffset,negOffset = negOffset,
                 sort = sortBed,unique = uniqueBed, memSize = 8 ,downSample = downSample,removeXS = TRUE)
    print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    return(R_sam2bed_wrapper(argv,filterList))
}

.lib_complex_qc_call <- function(bedfile, sortedBed, max_reads){

    argv <- list(bedfile = bedfile ,sortedBed = sortedBed,max_reads = max_reads)
    print(argv)
    rs<-lib_complex_qc(argv)
    if(rs["PBC2"]<0){
        rs["PBC2"]=NA
    }
    print(rs)
    return(rs)
}

# only chr1-chrY will be saved, chrM and others will be removed.
.chr_separate_call <- function(ReadsIfile, ReadsOpath, Name){
  argv <- list(readsIfile = ReadsIfile, readsOpath = ReadsOpath, name = Name)
  ChrDivi_wrapper(argv)
  return(TRUE)
}

.CutSite_call <- function(InputFile, OutputFile){
  argv <- list(readsIfile = InputFile, readsOpath = OutputFile)
  print(argv)
  return(CutCountPre_wrapper(argv))
}


.CutSiteCount <- function(readsfile, motiffile, matrixfile, motif_len, strand_len){
  argv <- list(readsfile = readsfile, motiffile = motiffile, matrixfile = matrixfile,
               motif_len = motif_len, strand_len = strand_len)
  return(CutSiteCount_wrapper(argv))
}




.bedOprUtils_call<-function(ibedfile, obedfile,reportPrefix,
                            mergePair,downSample ,posOffset,negOffset,
                            sortBed,uniqueBed,minFregLen,maxFregLen,filterList,select){
    argv <- list(ibedfile=ibedfile, obedfile=obedfile,reportPrefix=reportPrefix,
                 memSize=8,mergePair=mergePair,downSample = downSample,
                 posOffset=posOffset,negOffset=negOffset,
                 sortBed=sortBed,uniqueBed=uniqueBed,
                 minFregLen=minFregLen,maxFregLen=maxFregLen,
                 filterList=filterList,select=select)
    print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    return(bedOprUtils(argv,filterList))

}


checkFileExist <- function(filePath,argname){
    if(!is.null(filePath)){
        for(i in 1:length(filePath)){
            if(!file.exists(filePath[i])){
                stop(sprintf("For argument `%s`, file does not exist: `%s`",
                             argname,filePath[i]))
            }
        }
    }
}
checkPathExist <- function(filePath,argname){
    if(!is.null(filePath)){
        if(!dir.exists(dirname(filePath))){
            stop(sprintf("For argument `%s`,path does not exist: `%s`",
                         argname,filePath))
        }
    }
}
checkFileCreatable <- function(filePath,argname,overwrite){
    if(!is.null(filePath)){
        if(file.exists(filePath)){
            if(overwrite){
                warning(sprintf(paste0("For argument `%s`, file exist:%s. ",
                                       "It will be overwrited"),
                                argname,filePath));
            }else{
                stop(sprintf(paste0("For argument `%s`,file exist: %s. ",
                                    "Use 'overwrite=TRUE' to overwrite"),
                             argname,filePath));
            }
        }else if(!file.create(filePath)){
            stop(sprintf(paste0("For argument `%s`, cannot create file `%s`.",
                                "\nNo such directory or permission denied."),
                         argname,filePath));
            stop("")
        }else{
            unlink(filePath)
        }
    }
}

checkAddArgus<- function(pattern,...){
    paramlist<-trimws(as.character(list(...)))
    paramArray<-c()
    if(length(paramlist)>0){
        for(i in 1:length(paramlist)){
            paramArray<-c(paramArray,strsplit(paramlist[i],"\\s+")[[1]])
        }
    }
    fixed<-grepl(pattern,paramArray)
    if(sum(fixed)>0){
        invalidp<-paste0(paramArray[fixed],collapse = " ")
        stop(sprintf(paste0("Argument(s) `%s` are invalid for additional ",
                            "argument. Please set them as fixed arguments."),
                     invalidp))
    }
    return(paramArray)
}


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
#' if(file.exists(file.path(td, "lambda_virus"))){
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
#' reads_1 <- system.file(package="Rbowtie2", "extdata", "adrm", "reads_1.fq")
#' reads_2 <- system.file(package="Rbowtie2", "extdata", "adrm", "reads_2.fq")
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
#' reads_1 <- system.file(package="Rbowtie2", "extdata", "adrm", "reads_1.fq")
#' reads_2 <- system.file(package="Rbowtie2", "extdata", "adrm", "reads_2.fq")
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
    call <- paste(shQuote(file.path(system.file(package="Rbowtie2"), bin)), args)
    output <- system(call, intern=TRUE,show.output.on.console=TRUE)
    return(output)
}




