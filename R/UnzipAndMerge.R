setClass(Class = "UnzipAndMerge",
         contains = "ATACProc"
         )

setMethod(
    f = "init",
    signature = "UnzipAndMerge",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        fastqInput1 <- allparam[["fastqInput1"]]
        fastqInput2 <- allparam[["fastqInput2"]]
        fastqOutput1 <- allparam[["fastqOutput1"]]
        fastqOutput2 <- allparam[["fastqOutput2"]]
        interleave <- allparam[["interleave"]]
        param(.Object)$interleave <- interleave
        property(.Object)$interleave <- interleave
        if(interleave){
            property(.Object)$singleEnd <- FALSE
            input(.Object)$fastqInput1 <- fastqInput1
            if(!is.null(fastqOutput1)){
                output(.Object)$fastqOutput1 <- fastqOutput1
            }else{
                output(.Object)$fastqOutput1 <- getAutoPath(.Object,.Object$inputList[["fastqInput1"]][1], "(fastq|fq|bz2|gz)","fq")
            }
        }else{
            if(is.null(fastqInput2)){
                property(.Object)$singleEnd <- TRUE
                input(.Object)$fastqInput1 <-fastqInput1
                if(!is.null(fastqOutput1)){
                    output(.Object)$fastqOutput1 <-fastqOutput1
                }else{
                    output(.Object)$fastqOutput1 <-getAutoPath(.Object,.Object$inputList[["fastqInput1"]][1], "(fastq|fq|bz2|gz)","fq")
                }
            }else{
                property(.Object)$singleEnd <-FALSE
                input(.Object)$fastqInput1 <-fastqInput1
                input(.Object)$fastqInput2 <-fastqInput2
                if(length(.Object$inputList[["fastqInput1"]])!=length(.Object$inputList[["fastqInput2"]])){
                    stop("The number of pair-end fastq files should be equal.")
                }
                if(!is.null(fastqOutput1)){
                    #add check of private$paramlist[["fastqInput1"]][i]!=fastqOutput1
                    output(.Object)$fastqOutput1 <- fastqOutput1
                }else{
                    output(.Object)$fastqOutput1 <- getAutoPath(.Object,.Object$inputList[["fastqInput1"]][1], "(fastq|fq|bz2|gz)","fq")
                }
                if(!is.null(fastqOutput2)){
                    #add check of private$paramlist[["fastqInput2"]][i]!=fastqOutput2
                    output(.Object)$fastqOutput2 <- fastqOutput2
                }else{
                    output(.Object)$fastqOutput2 <- getAutoPath(.Object,.Object$inputList[["fastqInput2"]][1], "(fastq|fq|bz2|gz)","fq")
                }
            }
        }
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "UnzipAndMerge",
    definition = function(.Object,...){
        if(property(.Object)[["singleEnd"]]||(!property(.Object)[["singleEnd"]]&&.Object$paramList[["interleave"]])){
            fileNumber<-length(.Object$inputList[["fastqInput1"]])
            decompress(.Object,.Object$inputList[["fastqInput1"]][1],.Object$outputList[["fastqOutput1"]])
            if(fileNumber>1){
                for(i in 2:fileNumber){
                    tempfastqfile<-decompressFastq(.Object,.Object$inputList[["fastqInput1"]][i],dirname(.Object$outputList[["fastqOutput1"]]));
                    file.append(.Object$outputList[["fastqOutput1"]],tempfastqfile)
                    if(tempfastqfile!=.Object$inputList[["fastqInput1"]][i]){
                        unlink(tempfastqfile)
                    }
                }
            }
        }else{
            fileNumber<-length(.Object$inputList[["fastqInput1"]])
            decompress(.Object,.Object$inputList[["fastqInput1"]][1],.Object$outputList[["fastqOutput1"]])
            decompress(.Object,.Object$inputList[["fastqInput2"]][1],.Object$outputList[["fastqOutput2"]])
            if(fileNumber>1){
                for(i in 2:fileNumber){
                    tempfastqfile<-decompressFastq(.Object,.Object$inputList[["fastqInput1"]][i],dirname(.Object$outputList[["fastqOutput1"]]));
                    file.append(.Object$outputList[["fastqOutput1"]],tempfastqfile)
                    if(tempfastqfile!=.Object$inputList[["fastqInput1"]][i]){
                        unlink(tempfastqfile)
                    }
                    tempfastqfile<-decompressFastq(.Object,.Object$inputList[["fastqInput2"]][i],dirname(.Object$outputList[["fastqOutput2"]]));
                    file.append(.Object$outputList[["fastqOutput2"]],tempfastqfile)
                    if(tempfastqfile!=.Object$inputList[["fastqInput2"]][i]){
                        unlink(tempfastqfile)
                    }
                }
            }
        }
        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "UnzipAndMerge",
    definition = function(.Object,...){
        if(is.null(.Object$inputList[["fastqInput1"]])){
            stop("fastqInput1 is required.")
        }
        if(.Object$paramList[["interleave"]]&&.Object$propList[["singleEnd"]]){
            stop("Single end data should not be interleave")
        }
    }
)




setGeneric(
    name = "decompressFastq",
    def = function(.Object,filename,destpath,...){
        standardGeneric("decompressFastq")
    }
)
setMethod(
    f = "decompressFastq",
    signature = "UnzipAndMerge",
    definition = function(.Object,filename,destpath,...){
        destname<-file.path(destpath,basename(filename))
        writeLog(.Object,paste0("processing file:"))
        writeLog(.Object,sprintf("source:%s",filename))
        writeLog(.Object,sprintf("destination:%s",destname))
        if(isBzipped(filename)){
            destname<-gsub(sprintf("[.]%s$", "bz2"), "", destname, ignore.case=TRUE)
            return(bunzip2(filename,destname=destname,overwrite=TRUE,remove=FALSE))
        }else if(isGzipped(filename)){
            destname<-gsub(sprintf("[.]%s$", "gz"), "", destname, ignore.case=TRUE)
            return(gunzip(filename,destname=destname,overwrite=TRUE,remove=FALSE))
        }else{
            return(filename)
        }
    }
)

setGeneric(
    name = "decompress",
    def = function(.Object,filename,destname,...){
        standardGeneric("decompress")
    }
)
setMethod(
    f = "decompress",
    signature = "UnzipAndMerge",
    definition = function(.Object,filename,destname,...){
        writeLog(.Object,paste0("processing file:"))
        writeLog(.Object,sprintf("source:%s",filename))
        writeLog(.Object,sprintf("destination:%s",destname))
        if(isBzipped(filename)){
            return(bunzip2(filename,destname=destname,overwrite=TRUE,remove=FALSE))
        }else if(isGzipped(filename)){
            return(gunzip(filename,destname=destname,overwrite=TRUE,remove=FALSE))
        }else if(normalizePath(dirname(filename))!=normalizePath(dirname(destname))||
                 basename(filename)!=basename(destname)){
            file.copy(filename,destname,overwrite = TRUE)
        }

        return(destname)
    }
)


setGeneric(
    name = "removeCompressSuffix",
    def = function(.Object,filename,...){
        standardGeneric("removeCompressSuffix")
    }
)
setMethod(
    f = "removeCompressSuffix",
    signature = "UnzipAndMerge",
    definition = function(.Object,filename,...){
        filename<-gsub(sprintf("[.]%s$", "bz2"), "", filename, ignore.case=TRUE)
        filename<-gsub(sprintf("[.]%s$", "gz"), "", filename, ignore.case=TRUE)
        filename<-gsub(sprintf("[.]%s$", "fastq"), "", filename, ignore.case=TRUE)
        filename<-gsub(sprintf("[.]%s$", "fq"), "", filename, ignore.case=TRUE)
        filename<-paste0(filename,".",getProcName(.Object),".fq")
        return(filename)
    }
)



#' @name UnzipAndMerge
#' @aliases atacUnzipAndMerge
#' @aliases removeAdapter
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
#' @param fastqOutput1 \code{Character}. The trimmed mate1 reads output file
#' path for fastqInput2.
#' @param fastqOutput2 \code{Character}. The trimmed mate2 reads output file
#' path for fastqInput2.
#' @param interleave \code{Logical}. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param ... Additional arguments, currently unused.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacRenamer}}
#' \code{\link{atacQCReport}}
#' @examples
#' td<-tempdir()
#' setTmpDir(td)
#' 
#' # Identify adapters
#' prefix<-system.file(package="esATAC", "extdata", "uzmg")
#' (reads_1 <-file.path(prefix,"m1",dir(file.path(prefix,"m1"))))
#' (reads_2 <-file.path(prefix,"m2",dir(file.path(prefix,"m2"))))
#' 
#' reads_merged_1 <- file.path(td,"reads_1.fq")
#' reads_merged_2 <- file.path(td,"reads_2.fq")
#' atacproc <- atacUnzipAndMerge(fastqInput1 = reads_1,fastqInput2 = reads_2)
#' dir(td)
#' 
#' @rdname UnzipAndMerge
#' @export
atacUnzipAndMerge<- function(fastqInput1, fastqInput2=NULL,
                             fastqOutput1=NULL,fastqOutput2=NULL,
                             interleave = FALSE, ...){
    allpara <- c(list(Class = "UnzipAndMerge", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}

#' @rdname UnzipAndMerge
#' @export
unzipAndMerge<- function(fastqInput1, fastqInput2=NULL,
                             fastqOutput1=NULL,fastqOutput2=NULL,
                             interleave = FALSE, ...){
    allpara <- c(list(Class = "UnzipAndMerge", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
    
}

