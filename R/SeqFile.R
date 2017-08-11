UnzipAndMerge <-R6Class(
    classname = "UnzipAndMerge",
    inherit = BaseProc,
    public = list(
        initialize = function(fastqInput1, fastqInput2=NULL,fastqOutput1=NULL,fastqOutput2=NULL,editable=FALSE){
            super$initialize("UnzipAndMerge",editable,list())
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
                    private$paramlist[["fastqOutput1"]]<-fastqOutput1
                    private$checkFileCreatable(private$paramlist[["fastqOutput1"]])
                }else{
                    private$paramlist[["fastqOutput1"]]<-file.path(.obtainConfigure("tmpdir"),basename(private$paramlist[["fastqInput1"]][1]))
                    private$paramlist[["fastqOutput1"]]<-private$removeCompressSuffix(private$paramlist[["fastqOutput1"]])
                }
                if(!is.null(fastqOutput2)){
                    private$paramlist[["fastqOutput2"]]<-fastqOutput2
                    private$checkFileCreatable(private$paramlist[["fastqOutput2"]])
                }else{
                    private$paramlist[["fastqOutput2"]]<-file.path(.obtainConfigure("tmpdir"),basename(private$paramlist[["fastqInput2"]][1]))
                    private$paramlist[["fastqOutput2"]]<-private$removeCompressSuffix(private$paramlist[["fastqOutput2"]])
                }
            }

            private$checkRequireParam();

        },
        processing = function(){
            if(!super$processing()){
                return()
            }
            if(private$singleEnd){
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
            private$setFinish()
        },
        setResultParam = function(fastqInput1, fastqInput2=NULL){
            super$setResultParam();
            private$paramlist[["fastqInput1"]] <- fastqInput1
            private$paramlist[["fastqInput2"]] <- fastqInput2
            private$checkFileExist(private$paramlist[["fastqInput1"]])
            private$checkFileExist(private$paramlist[["fastqInput2"]])
            private$paramlist[["fastqOutput1"]] <- fastqInput1
            private$paramlist[["fastqOutput2"]] <- fastqInput2
        }
    ),
    private = list(
        checkRequireParam = function(){
            if(private$editable){
                return();
            }
            if(is.null(private$paramlist[["fastqInput1"]])){
                stop("fastqInput1 is required.")
            }
        },
        decompressFastq = function(filename,destpath){
            destname<-file.path(destpath,basename(filename))

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
            if(isBzipped(filename)){
                return(bunzip2(filename,destname=destname,overwrite=TRUE,remove=FALSE))
            }else if(isGzipped(filename)){
                return(gunzip(filename,destname=destname,overwrite=TRUE,remove=FALSE))
            }else if(normalizePath(dirname(filename))!=normalizePath(dirname(destname))||
                     basename(filename)!=basename(destname)){
                file.copy(filename,destname)
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





atacUnzipAndMerge<- function(fastqInput1,fastqInput2=NULL){
    atacproc <- UnzipAndMerge$new(fastqInput1,fastqInput2);
    atacproc$processing();
    return(atacproc);
}
