setClass(Class = "SCUnzipAndMerge",
         contains = "ATACProc"
         )

setMethod(
    f = "init",
    signature = "SCUnzipAndMerge",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        fastqInput1 <- allparam[["fastqInput1"]]
        fastqInput2 <- allparam[["fastqInput2"]]
        fastqBarcodeInput <- allparam[["fastqBarcodeInput"]] 
        fastqOutput1 <- allparam[["fastqOutput1"]]
        fastqOutput2 <- allparam[["fastqOutput2"]]
        fastqBarcodeOutput <- allparam[["fastqBarcodeOutput"]]
        interleave <- FALSE
        input(.Object)$fastqInput1 <-fastqInput1
        input(.Object)$fastqInput2 <-fastqInput2
        if(length(input(.Object)[["fastqInput1"]])!=length(input(.Object)[["fastqInput2"]])){
            stop("The number of pair-end fastq files should be equal.")
        }
        if(length(input(.Object)[["fastqInput1"]])!=length(input(.Object)[["fastqBarcodeInput"]])){
            stop("The number of Barcode fastq files should match fastq.")
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
        if(!is.null(fastqBarcodeOutput)){
            output(.Object)$fastqBarcodeOutput <- fastqBarcodeOutput
        }else{
            output(.Object)$fastqBarcodeOutput <- getAutoPath(.Object,input(.Object)[["fastqBarcodeInput"]][1], "(fastq|fq|bz2|gz)","fq")
        }
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "SCUnzipAndMerge",
    definition = function(.Object,...){
        fileNumber<-length(input(.Object)[["fastqInput1"]])
        decompress(.Object,input(.Object)[["fastqInput1"]][1],output(.Object)[["fastqOutput1"]])
        decompress(.Object,input(.Object)[["fastqBarcodeInput"]][1],output(.Object)[["fastqBarcodeOutput"]])
        decompress(.Object,input(.Object)[["fastqInput2"]][1],output(.Object)[["fastqOutput2"]])

        if(fileNumber>1){
            for(i in 2:fileNumber){
                tempfastqfile<-decompressFastq(.Object,input(.Object)[["fastqInput1"]][i],dirname(output(.Object)[["fastqOutput1"]]));
                file.append(output(.Object)[["fastqOutput1"]],tempfastqfile)
                if(tempfastqfile!=input(.Object)[["fastqInput1"]][i]){
                    unlink(tempfastqfile)
                }
                tempfastqfile<-decompressFastq(.Object,input(.Object)[["fastqBarcodeInput"]][i],dirname(output(.Object)[["fastqBarcodeOutput"]]));
                file.append(output(.Object)[["fastqBarcodeOutput"]],tempfastqfile)
                if(tempfastqfile!=input(.Object)[["fastqBarcodeInput"]][i]){
                    unlink(tempfastqfile)
                }
                tempfastqfile<-decompressFastq(.Object,input(.Object)[["fastqInput2"]][i],dirname(output(.Object)[["fastqOutput2"]]));
                file.append(output(.Object)[["fastqOutput2"]],tempfastqfile)
                if(tempfastqfile!=input(.Object)[["fastqInput2"]][i]){
                    unlink(tempfastqfile)
                }
            }
        }
            
        .Object
    }
)



setMethod(
    f = "decompressFastq",
    signature = "SCUnzipAndMerge",
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

setMethod(
    f = "decompress",
    signature = "SCUnzipAndMerge",
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
    signature = "SCUnzipAndMerge",
    definition = function(.Object, ...){
        report(.Object)$seqtype  <- "paired end (PE)"
        report(.Object)$frag <- 2
        report(.Object)$filelist <- data.frame(`Mate1 files`=input(.Object)$fastqInput1,
                                                   `Mate2 files`=input(.Object)$fastqInput2,
`Barcode files`=input(.Object)$fastqBarcodeInput,
)
        
        .Object
    }
)


#' @name SCUnzipAndMerge
#' @title Unzip and merge fastq files
#' @description
#' Unzip and merge fastq files that are in format of bzip, gzip or fastq
#' @param fastqInput1 \code{Character} vector. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates paired
#' with file paths in fastqInput2
#' @param fastqBarcodeInput \code{Character} vector. Fastq file that save barcode.
#' @param fastqInput2 \code{Character} vector. It contains file paths with #2
#' mates paired with file paths in fastqInput1
#' @param fastqOutput1 \code{Character}. The trimmed mate1 reads output file
#' path for fastqInput1.
#' @param fastqBarcodeOutput \code{Character}. The trimmed mate1 reads output file
#' path for fastqInput1.
#' @param fastqOutput2 \code{Character}. The trimmed mate2 reads output file
#' path for fastqInput2.
#' @param ... Additional arguments, currently unused.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSCRenamer}}
#' @examples
#' 
#' ignoreCheck() # warnning: run this for fast test only
#' 
#' td<-tempdir()
#' setTmpDir(td)
#' 
# # Identify adapters
# prefix<-system.file(package="esATAC", "extdata", "uzmg")
# (reads_1 <-file.path(prefix,"m1",dir(file.path(prefix,"m1"))))
# (reads_2 <-file.path(prefix,"m2",dir(file.path(prefix,"m2"))))
# 
# reads_merged_1 <- file.path(td,"reads_1.fq")
# reads_merged_2 <- file.path(td,"reads_2.fq")
# atacproc <- atacUnzipAndMerge(fastqInput1 = reads_1,fastqInput2 = reads_2)
# dir(td)
# 
#' @rdname SCUnzipAndMerge
#' @export
atacSCUnzipAndMerge<- function(fastqInput1,fastqBarcodeInput, fastqInput2,
                             fastqOutput1=NULL,fastqBarcodeOutput=NULL,fastqOutput2=NULL, ...){
    allpara <- c(list(Class = "SCUnzipAndMerge", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}

#' @rdname SCUnzipAndMerge
#' @export
scUnzipAndMerge<- function(fastqInput1,fastqBarcodeInput, fastqInput2,
                             fastqOutput1=NULL,fastqBarcodeOutput=NULL,fastqOutput2=NULL, ...){
    allpara <- c(list(Class = "SCUnzipAndMerge", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
    
}

