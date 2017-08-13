Bowtie2Mapping <-R6Class(
    classname = "Bowtie2Mapping",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc,bowtie2Index=NULL,samOutput=NULL,threads=NULL,
                              fastqInput1=NULL, fastqInput2=NULL,editable=FALSE){
            super$initialize("Bowtie2Mapping",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1");
                private$paramlist[["fastqInput2"]] <- atacProc$getParam("fastqOutput2");
            }

            if(!is.null(fastqInput1)){
                private$paramlist[["fastqInput1"]] <- fastqInput1;
            }
            if(!is.null(fastqInput2)){
                private$paramlist[["fastqInput2"]] <- fastqInput2;
            }

            if(is.null(samOutput)){
                private$paramlist[["samOutput"]] <- paste(private$paramlist[["fastqInput1"]],".sam",sep="");
            }else{
                private$paramlist[["samOutput"]] <- samOutput;
            }

            if(is.null(bowtie2Index)){
                datadir<-.obtainConfigure("datadir");
                genome<-.obtainConfigure("genome");
                idxprefix<-paste0(datadir,"/",genome)
                private$paramlist[["bowtie2Index"]] <- idxprefix;
            }else{
                private$paramlist[["bowtie2Index"]] <- bowtie2Index;
            }

            if(is.null(threads)){
                private$paramlist[["threads"]] <- .obtainConfigure("threads");
            }else{
                private$paramlist[["threads"]] <- as.numeric(threads)
            }

            private$checkFileExist(private$paramlist[["fastqInput1"]]);
            private$checkFileExist(private$paramlist[["fastqInput2"]]);
            private$checkFileCreatable(private$paramlist[["samOutput"]]);
            private$checkFileExist(paste0(private$paramlist[["bowtie2Index"]],".1.bt2"));
            private$checkFileExist(paste0(private$paramlist[["bowtie2Index"]],".2.bt2"));
            private$checkFileExist(paste0(private$paramlist[["bowtie2Index"]],".3.bt2"));
            private$checkFileExist(paste0(private$paramlist[["bowtie2Index"]],".4.bt2"));
            private$checkFileExist(paste0(private$paramlist[["bowtie2Index"]],".rev.1.bt2"));
            private$checkFileExist(paste0(private$paramlist[["bowtie2Index"]],".rev.2.bt2"));
            private$checkRequireParam();
            print(private$paramlist[["samOutput"]])
        },
        processing = function(){
            super$processing()

            if(is.null(private$paramlist[["fastqInput2"]])){
                stop("Single-end-mapping does not support yet!")
            }else{
                .bowtie2_paired_call(bowtie2Index=private$paramlist[["bowtie2Index"]],
                                     samOutput=private$paramlist[["samOutput"]],
                                     fastqInput1=private$paramlist[["fastqInput1"]],
                                     fastqInput2=private$paramlist[["fastqInput2"]],
                                     threads=private$paramlist[["threads"]])
            }

            private$finish <- TRUE
        },
        setResultParam = function(samOutput){
            super$setResultParam();
            private$paramlist[["samOutput"]] <- samOutput
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
            if(is.null(private$paramlist[["bowtie2Index"]])){
                stop("bowtieIdx")
            }
        }
    )


)


