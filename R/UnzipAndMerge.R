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
                output(.Object)$fastqOutput1 <- getAutoPath(.Object,input(.Object)[["fastqInput1"]][1], "(fastq|fq|bz2|gz)","fq")
            }
        }else{
            if(is.null(fastqInput2)){
                property(.Object)$singleEnd <- TRUE
                input(.Object)$fastqInput1 <-fastqInput1
                if(!is.null(fastqOutput1)){
                    output(.Object)$fastqOutput1 <-fastqOutput1
                }else{
                    output(.Object)$fastqOutput1 <-getAutoPath(.Object,input(.Object)[["fastqInput1"]][1], "(fastq|fq|bz2|gz)","fq")
                }
            }else{
                property(.Object)$singleEnd <-FALSE
                input(.Object)$fastqInput1 <-fastqInput1
                input(.Object)$fastqInput2 <-fastqInput2
                if(length(input(.Object)[["fastqInput1"]])!=length(input(.Object)[["fastqInput2"]])){
                    stop("The number of pair-end fastq files should be equal.")
                }
                if(!is.null(fastqOutput1)){
                    #add check of private$paramlist[["fastqInput1"]][i]!=fastqOutput1
                    output(.Object)$fastqOutput1 <- fastqOutput1
                }else{
                    output(.Object)$fastqOutput1 <- getAutoPath(.Object,input(.Object)[["fastqInput1"]][1], "(fastq|fq|bz2|gz)","fq")
                }
                if(!is.null(fastqOutput2)){
                    #add check of private$paramlist[["fastqInput2"]][i]!=fastqOutput2
                    output(.Object)$fastqOutput2 <- fastqOutput2
                }else{
                    output(.Object)$fastqOutput2 <- getAutoPath(.Object,input(.Object)[["fastqInput2"]][1], "(fastq|fq|bz2|gz)","fq")
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
        if(property(.Object)[["singleEnd"]]||(!property(.Object)[["singleEnd"]]&&param(.Object)[["interleave"]])){
            fileNumber<-length(input(.Object)[["fastqInput1"]])
            decompress(.Object,input(.Object)[["fastqInput1"]][1],output(.Object)[["fastqOutput1"]])
            if(fileNumber>1){
                for(i in 2:fileNumber){
                    tempfastqfile<-decompressFastq(.Object,input(.Object)[["fastqInput1"]][i],dirname(output(.Object)[["fastqOutput1"]]));
                    file.append(output(.Object)[["fastqOutput1"]],tempfastqfile)
                    if(tempfastqfile!=input(.Object)[["fastqInput1"]][i]){
                        unlink(tempfastqfile)
                    }
                }
            }
        }else{
            fileNumber<-length(input(.Object)[["fastqInput1"]])
            decompress(.Object,input(.Object)[["fastqInput1"]][1],output(.Object)[["fastqOutput1"]])
            decompress(.Object,input(.Object)[["fastqInput2"]][1],output(.Object)[["fastqOutput2"]])
            if(fileNumber>1){
                for(i in 2:fileNumber){
                    tempfastqfile<-decompressFastq(.Object,input(.Object)[["fastqInput1"]][i],dirname(output(.Object)[["fastqOutput1"]]));
                    file.append(output(.Object)[["fastqOutput1"]],tempfastqfile)
                    if(tempfastqfile!=input(.Object)[["fastqInput1"]][i]){
                        unlink(tempfastqfile)
                    }
                    tempfastqfile<-decompressFastq(.Object,input(.Object)[["fastqInput2"]][i],dirname(output(.Object)[["fastqOutput2"]]));
                    file.append(output(.Object)[["fastqOutput2"]],tempfastqfile)
                    if(tempfastqfile!=input(.Object)[["fastqInput2"]][i]){
                        unlink(tempfastqfile)
                    }
                }
            }
        }
        
        .Object
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




setMethod(
    f = "genReport",
    signature = "UnzipAndMerge",
    definition = function(.Object, ...){
        if(param(.Object)$interleave){
            report(.Object)$seqtype <- "paired end (PE,interleave)"
            report(.Object)$frag <- 2
        }else if(is.null(input(.Object)$fastqInput2)){
            report(.Object)$seqtype  <- "single end (SE)"
            report(.Object)$frag <- 1
        }else{
            report(.Object)$seqtype  <- "paired end (PE)"
            report(.Object)$frag <- 2
        }
        if(is.null(input(.Object)$fastqInput2)){
            report(.Object)$filelist <- data.frame(`File(s)`= input(.Object)$fastqInput1)
        }else{
            report(.Object)$filelist <- data.frame(`Mate1 files`=input(.Object)$fastqInput1,
                                                   `Mate2 files`=input(.Object)$fastqInput2)
        }
        .Object
    }
)


#' @name UnzipAndMerge
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
#' 
#' ignoreCheck() # warnning: run this for fast test only
#' 
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

