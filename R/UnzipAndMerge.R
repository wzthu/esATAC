setClass(Class = "UnzipAndMerge",
         contains = "ATACProc"
         )

setMethod(
    f = "initialize",
    signature = "UnzipAndMerge",
    definition = function(.Object,fastqInput1,..., fastqInput2=NULL,fastqOutput1=NULL,fastqOutput2=NULL,interleave = FALSE,editable = FALSE){
        .Object <- init(.Object,"UnzipAndMerge",editable,list())
        .Object@paramlist[["interleave"]]<-interleave
        if(interleave){
            .Object@singleEnd<-FALSE
            .Object@paramlist[["fastqInput1"]]<-fastqInput1
            for(i in 1:length(.Object@paramlist[["fastqInput1"]])){
                checkFileExist(.Object,.Object@paramlist[["fastqInput1"]][i]);
            }
            if(!is.null(fastqOutput1)){
                .Object@paramlist[["fastqOutput1"]]<-fastqOutput1
                checkFileCreatable(.Object,.Object@paramlist[["fastqOutput1"]])
            }else{
                .Object@paramlist[["fastqOutput1"]]<-file.path(.obtainConfigure("tmpdir"),basename(.Object@paramlist[["fastqInput1"]][1]))
                .Object@paramlist[["fastqOutput1"]]<-removeCompressSuffix(.Object,.Object@paramlist[["fastqOutput1"]])
            }
        }else{
            if(is.null(fastqInput2)){
                .Object@singleEnd<-TRUE
                .Object@paramlist[["fastqInput1"]]<-fastqInput1
                for(i in 1:length(.Object@paramlist[["fastqInput1"]])){
                    checkFileExist(.Object,.Object@paramlist[["fastqInput1"]][i]);
                }
                if(!is.null(fastqOutput1)){
                    .Object@paramlist[["fastqOutput1"]]<-fastqOutput1
                    checkFileCreatable(.Object,.Object@paramlist[["fastqOutput1"]])
                }else{
                    .Object@paramlist[["fastqOutput1"]]<-file.path(.obtainConfigure("tmpdir"),basename(.Object@paramlist[["fastqInput1"]][1]))
                    .Object@paramlist[["fastqOutput1"]]<-removeCompressSuffix(.Object,.Object@paramlist[["fastqOutput1"]])
                }
                
            }else{
                .Object@singleEnd<-FALSE
                .Object@paramlist[["fastqInput1"]]<-fastqInput1
                for(i in 1:length(.Object@paramlist[["fastqInput1"]])){
                    checkFileExist(.Object,.Object@paramlist[["fastqInput1"]][i]);
                }
                .Object@paramlist[["fastqInput2"]]<-fastqInput2
                if(length(.Object@paramlist[["fastqInput1"]])!=length(.Object@paramlist[["fastqInput2"]])){
                    stop("The number of pair-end fastq files should be equal.")
                }
                for(i in 1:length(.Object@paramlist[["fastqInput2"]])){
                    checkFileExist(.Object,.Object@paramlist[["fastqInput2"]][i]);
                }
                if(!is.null(fastqOutput1)){
                    #add check of private$paramlist[["fastqInput1"]][i]!=fastqOutput1
                    .Object@paramlist[["fastqOutput1"]]<-fastqOutput1
                    checkFileCreatable(.Object,.Object@paramlist[["fastqOutput1"]])
                }else{
                    .Object@paramlist[["fastqOutput1"]]<-file.path(.obtainConfigure("tmpdir"),basename(.Object@paramlist[["fastqInput1"]][1]))
                    .Object@paramlist[["fastqOutput1"]]<-removeCompressSuffix(.Object,.Object@paramlist[["fastqOutput1"]])
                }
                if(!is.null(fastqOutput2)){
                    #add check of private$paramlist[["fastqInput2"]][i]!=fastqOutput2
                    .Object@paramlist[["fastqOutput2"]]<-fastqOutput2
                    checkFileCreatable(.Object,.Object@paramlist[["fastqOutput2"]])
                }else{
                    .Object@paramlist[["fastqOutput2"]]<-file.path(.obtainConfigure("tmpdir"),basename(.Object@paramlist[["fastqInput2"]][1]))
                    .Object@paramlist[["fastqOutput2"]]<-removeCompressSuffix(.Object,.Object@paramlist[["fastqOutput2"]])
                }
            }
        }
        
        paramValidation(.Object)
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "UnzipAndMerge",
    definition = function(.Object,...){
        if(.Object@singleEnd||(!.Object@singleEnd&&.Object@paramlist[["interleave"]])){
            fileNumber<-length(.Object@paramlist[["fastqInput1"]])
            decompress(.Object,.Object@paramlist[["fastqInput1"]][1],.Object@paramlist[["fastqOutput1"]])
            if(fileNumber>1){
                for(i in 2:fileNumber){
                    tempfastqfile<-decompressFastq(.Object,.Object@paramlist[["fastqInput1"]][i],dirname(.Object@paramlist[["fastqOutput1"]]));
                    file.append(.Object@paramlist[["fastqOutput1"]],tempfastqfile)
                    if(tempfastqfile!=.Object@paramlist[["fastqInput1"]][i]){
                        unlink(tempfastqfile)
                    }
                }
            }
        }else{
            fileNumber<-length(.Object@paramlist[["fastqInput1"]])
            decompress(.Object,.Object@paramlist[["fastqInput1"]][1],.Object@paramlist[["fastqOutput1"]])
            decompress(.Object,.Object@paramlist[["fastqInput2"]][1],.Object@paramlist[["fastqOutput2"]])
            if(fileNumber>1){
                for(i in 2:fileNumber){
                    tempfastqfile<-decompressFastq(.Object,.Object@paramlist[["fastqInput1"]][i],dirname(.Object@paramlist[["fastqOutput1"]]));
                    file.append(.Object@paramlist[["fastqOutput1"]],tempfastqfile)
                    if(tempfastqfile!=.Object@paramlist[["fastqInput1"]][i]){
                        unlink(tempfastqfile)
                    }
                    tempfastqfile<-decompressFastq(.Object,.Object@paramlist[["fastqInput2"]][i],dirname(.Object@paramlist[["fastqOutput2"]]));
                    file.append(.Object@paramlist[["fastqOutput2"]],tempfastqfile)
                    if(tempfastqfile!=.Object@paramlist[["fastqInput2"]][i]){
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
        if(is.null(.Object@paramlist[["fastqInput1"]])){
            stop("fastqInput1 is required.")
        }
        if(.Object@paramlist[["interleave"]]&&.Object@singleEnd){
            stop("Single end data should not be interleave")
        }
    }
)



setMethod(
    f = "checkAllPath",
    signature = "UnzipAndMerge",
    definition = function(.Object,...){
        checkFileCreatable(.Object,.Object@paramlist[["fastqOutput1"]])
        checkFileCreatable(.Object,.Object@paramlist[["fastqOutput2"]])
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
        .Object<-writeLog(.Object,paste0("processing file:"))
        .Object<-writeLog(.Object,sprintf("source:%s",filename))
        .Object<-writeLog(.Object,sprintf("destination:%s",destname))
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
        .Object<-writeLog(.Object,paste0("processing file:"))
        .Object<-writeLog(.Object,sprintf("source:%s",filename))
        .Object<-writeLog(.Object,sprintf("destination:%s",destname))
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






UnzipAndMerge <-R6Class(
    classname = "UnzipAndMerge",
    inherit = ATACProc,
    public = list(
        initialize = function(fastqInput1, fastqInput2=NULL,fastqOutput1=NULL,fastqOutput2=NULL,interleave = FALSE,editable = FALSE){
            super$initialize("UnzipAndMerge",editable,list())
            private$paramlist[["interleave"]]<-interleave
            if(interleave){
                private$singleEnd<-FALSE
                private$paramlist[["fastqInput1"]]<-fastqInput1
                for(i in 1:length(private$paramlist[["fastqInput1"]])){
                    private$checkFileExist(private$paramlist[["fastqInput1"]][i]);
                }
                if(!is.null(fastqOutput1)){
                    private$paramlist[["fastqOutput1"]]<-fastqOutput1
                    private$checkFileCreatable(private$paramlist[["fastqOutput1"]])
                }else{
                    private$paramlist[["fastqOutput1"]]<-file.path(.obtainConfigure("tmpdir"),basename(private$paramlist[["fastqInput1"]][1]))
                    private$paramlist[["fastqOutput1"]]<-private$removeCompressSuffix(private$paramlist[["fastqOutput1"]])
                }
            }else{
                if(is.null(fastqInput2)){
                    private$singleEnd<-TRUE
                    private$paramlist[["fastqInput1"]]<-fastqInput1
                    for(i in 1:length(private$paramlist[["fastqInput1"]])){
                        private$checkFileExist(private$paramlist[["fastqInput1"]][i]);
                    }
                    if(!is.null(fastqOutput1)){
                        private$paramlist[["fastqOutput1"]]<-fastqOutput1
                        private$checkFileCreatable(private$paramlist[["fastqOutput1"]])
                    }else{
                        private$paramlist[["fastqOutput1"]]<-file.path(.obtainConfigure("tmpdir"),basename(private$paramlist[["fastqInput1"]][1]))
                        private$paramlist[["fastqOutput1"]]<-private$removeCompressSuffix(private$paramlist[["fastqOutput1"]])
                    }
                    
                }else{
                    private$singleEnd<-FALSE
                    private$paramlist[["fastqInput1"]]<-fastqInput1
                    for(i in 1:length(private$paramlist[["fastqInput1"]])){
                        private$checkFileExist(private$paramlist[["fastqInput1"]][i]);
                    }
                    private$paramlist[["fastqInput2"]]<-fastqInput2
                    if(length(private$paramlist[["fastqInput1"]])!=length(private$paramlist[["fastqInput2"]])){
                        stop("The number of pair-end fastq files should be equal.")
                    }
                    for(i in 1:length(private$paramlist[["fastqInput2"]])){
                        private$checkFileExist(private$paramlist[["fastqInput2"]][i]);
                    }
                    if(!is.null(fastqOutput1)){
                        #add check of private$paramlist[["fastqInput1"]][i]!=fastqOutput1
                        private$paramlist[["fastqOutput1"]]<-fastqOutput1
                        private$checkFileCreatable(private$paramlist[["fastqOutput1"]])
                    }else{
                        private$paramlist[["fastqOutput1"]]<-file.path(.obtainConfigure("tmpdir"),basename(private$paramlist[["fastqInput1"]][1]))
                        private$paramlist[["fastqOutput1"]]<-private$removeCompressSuffix(private$paramlist[["fastqOutput1"]])
                    }
                    if(!is.null(fastqOutput2)){
                        #add check of private$paramlist[["fastqInput2"]][i]!=fastqOutput2
                        private$paramlist[["fastqOutput2"]]<-fastqOutput2
                        private$checkFileCreatable(private$paramlist[["fastqOutput2"]])
                    }else{
                        private$paramlist[["fastqOutput2"]]<-file.path(.obtainConfigure("tmpdir"),basename(private$paramlist[["fastqInput2"]][1]))
                        private$paramlist[["fastqOutput2"]]<-private$removeCompressSuffix(private$paramlist[["fastqOutput2"]])
                    }
                }
            }
            
            private$paramValidation()


        }
    ),
    private = list(
        processing = function(){
            if(private$singleEnd||(!private$singleEnd&&private$paramlist[["interleave"]])){
                fileNumber<-length(private$paramlist[["fastqInput1"]])
                private$decompress(private$paramlist[["fastqInput1"]][1],private$paramlist[["fastqOutput1"]])
                if(fileNumber>1){
                    for(i in 2:fileNumber){
                        tempfastqfile<-private$decompressFastq(private$paramlist[["fastqInput1"]][i],dirname(private$paramlist[["fastqOutput1"]]));
                        file.append(private$paramlist[["fastqOutput1"]],tempfastqfile)
                        if(tempfastqfile!=private$paramlist[["fastqInput1"]][i]){
                            unlink(tempfastqfile)
                        }
                    }
                }
            }else{
                fileNumber<-length(private$paramlist[["fastqInput1"]])
                private$decompress(private$paramlist[["fastqInput1"]][1],private$paramlist[["fastqOutput1"]])
                private$decompress(private$paramlist[["fastqInput2"]][1],private$paramlist[["fastqOutput2"]])
                if(fileNumber>1){
                    for(i in 2:fileNumber){
                        tempfastqfile<-private$decompressFastq(private$paramlist[["fastqInput1"]][i],dirname(private$paramlist[["fastqOutput1"]]));
                        file.append(private$paramlist[["fastqOutput1"]],tempfastqfile)
                        if(tempfastqfile!=private$paramlist[["fastqInput1"]][i]){
                            unlink(tempfastqfile)
                        }
                        tempfastqfile<-private$decompressFastq(private$paramlist[["fastqInput2"]][i],dirname(private$paramlist[["fastqOutput2"]]));
                        file.append(private$paramlist[["fastqOutput2"]],tempfastqfile)
                        if(tempfastqfile!=private$paramlist[["fastqInput2"]][i]){
                            unlink(tempfastqfile)
                        }
                    }
                }
            }
        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["fastqInput1"]])){
                stop("fastqInput1 is required.")
            }
            if(private$paramlist[["interleave"]]&&private$singleEnd){
                stop("Single end data should not be interleave")
            }
        },
        checkAllPath = function(){
            private$checkFileCreatable(private$paramlist[["fastqOutput1"]])
            private$checkFileCreatable(private$paramlist[["fastqOutput2"]])
        },
        decompressFastq = function(filename,destpath){
            destname<-file.path(destpath,basename(filename))
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("source:%s",filename))
            private$writeLog(sprintf("destination:%s",destname))
            if(isBzipped(filename)){
                destname<-gsub(sprintf("[.]%s$", "bz2"), "", destname, ignore.case=TRUE)
                return(bunzip2(filename,destname=destname,overwrite=TRUE,remove=FALSE))
            }else if(isGzipped(filename)){
                destname<-gsub(sprintf("[.]%s$", "gz"), "", destname, ignore.case=TRUE)
                return(gunzip(filename,destname=destname,overwrite=TRUE,remove=FALSE))
            }else{
                return(filename)
            }


        },
        decompress = function(filename,destname){
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("source:%s",filename))
            private$writeLog(sprintf("destination:%s",destname))
            if(isBzipped(filename)){
                return(bunzip2(filename,destname=destname,overwrite=TRUE,remove=FALSE))
            }else if(isGzipped(filename)){
                return(gunzip(filename,destname=destname,overwrite=TRUE,remove=FALSE))
            }else if(normalizePath(dirname(filename))!=normalizePath(dirname(destname))||
                     basename(filename)!=basename(destname)){
                file.copy(filename,destname,overwrite = TRUE)
            }

            return(destname)
        },
        removeCompressSuffix= function(filename){
            filename<-gsub(sprintf("[.]%s$", "bz2"), "", filename, ignore.case=TRUE)
            filename<-gsub(sprintf("[.]%s$", "gz"), "", filename, ignore.case=TRUE)
            filename<-gsub(sprintf("[.]%s$", "fastq"), "", filename, ignore.case=TRUE)
            filename<-gsub(sprintf("[.]%s$", "fq"), "", filename, ignore.case=TRUE)
            filename<-paste0(filename,".",self$getProcName(),".fq")
            return(filename)
        }
    )

)




#' @name atacUnzipAndMerge
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
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacRenamer}} 
#' \code{\link{atacQCReport}} 
#' @examples 
#' td <- tempdir()
#' options(atacConf=setConfigure("tmpdir",td))
#' 
#' # Identify adapters
#' prefix<-system.file(package="ATACpipe", "extdata", "uzmg")
#' (reads_1 <-file.path(prefix,"m1",dir(file.path(prefix,"m1"))))
#' (reads_2 <-file.path(prefix,"m2",dir(file.path(prefix,"m2"))))
#' 
#' reads_merged_1 <- file.path(td,"reads1.fastq")
#' reads_merged_2 <- file.path(td,"reads2.fastq")
#' atacproc <- atacUnzipAndMerge(fastqInput1 = reads_1,fastqInput2 = reads_2) 
#' dir(td)
#' 
#' @rdname atacUnzipAndMerge
#' @export 
atacUnzipAndMerge<- function(fastqInput1, fastqInput2=NULL,
                             fastqOutput1=NULL,fastqOutput2=NULL,
                             interleave = FALSE){
    atacproc <- new(
        "UnzipAndMerge",
        fastqInput1 = fastqInput1,
        fastqInput2 = fastqInput2,
        fastqOutput1 = fastqOutput1,
        fastqOutput2 = fastqOutput2,
        interleave = interleave);
    atacproc<-process(atacproc)
    invisible(atacproc)
}

#' @rdname atacUnzipAndMerge
#' @export 
unzipAndMerge<- function(fastqInput1, fastqInput2=NULL,
                             fastqOutput1=NULL,fastqOutput2=NULL,
                             interleave = FALSE){
    atacproc <- new(
        "UnzipAndMerge",
        fastqInput1 = fastqInput1,
        fastqInput2 = fastqInput2,
        fastqOutput1 = fastqOutput1,
        fastqOutput2 = fastqOutput2,
        interleave = interleave);
    atacproc<-process(atacproc)
    invisible(atacproc)
}

