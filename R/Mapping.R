Bowtie2Mapping <-R6Class(
    classname = "Bowtie2Mapping",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc,samOutput=NULL, bt2Idx=NULL, fastqInput1=NULL, fastqInput2=NULL, interleave = FALSE, paramList="default",reportPrefix =NULL,editable=FALSE){
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

            if(.obtainConfigure("threads")>1){
                private$paramlist[["paramList"]]<-c("-p",as.character(.obtainConfigure("threads")))
            }
            if(is.null(paramList)){

            }else if(paramList=="default"){
                private$paramlist[["paramList"]]<-c("--no-discordant","--no-unal","--no-mixed","-X","2000",
                                                    private$paramlist[["paramList"]])
            }else{
                private$paramlist[["paramList"]]<-c(paramList,private$paramlist[["paramList"]])
                rejectp<-"-p|--threads|-x|-1|-2|-U|-S|--interleaved"
                private$checkParam(paramlist,rejectp)
            }

            if(is.null(reportPrefix)){
                if(!is.null(private$paramlist[["fastqInput1"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput1"]],regexProcName)
                    private$paramlist[["reportPrefix"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
            }else{
                private$paramlist[["reportPrefix"]] <- reportPrefix;
            }
            private$paramValidation()
        }


    ),
    private = list(
        processing = function(){
            private$writeLog("start mapping with parameters:")
            private$writeLog(paste0("bowtie2 index:",private$paramlist[["bt2Idx"]]))
            private$writeLog(paste0("samOutput:",private$paramlist[["samOutput"]]))
            private$writeLog(paste0("report:",private$paramlist[["reportPrefix"]]))
            private$writeLog(paste0("fastqInput1:",private$paramlist[["fastqInput1"]]))
            private$writeLog(paste0("fastqInput2:",private$paramlist[["fastqInput2"]]))
            private$writeLog(paste0("other parameters:",paste(private$paramlist[["paramList"]],collapse = " ")))
            sink(private$paramlist[["reportPrefix"]])####------------------------------------------------------------std::cout>>Rcpp::Rout
#            .bowtie2_call(bowtie2Index=private$paramlist[["bt2Idx"]],samOutput=private$paramlist[["samOutput"]],
#                          fastqInput1=private$paramlist[["fastqInput1"]],fastqInput2=private$paramlist[["fastqInput2"]],
#                          paramlist=private$paramlist[["paramList"]])

            bowtie2(bt2Index = private$paramlist[["bt2Idx"]],
                    samOutput = private$paramlist[["samOutput"]],
                    seq1 = private$paramlist[["fastqInput1"]],
                    paste(private$paramlist[["paramList"]],collapse = " "),
                    seq2 = private$paramlist[["fastqInput2"]],
                    interleaved = private$paramlist[["interleave"]],
                    overwrite=TRUE)
            sink()



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
        }
    )


)


atacBowtie2Mapping <- function(atacProc,samOutput=NULL,reportPrefix =NULL, bt2Idx=NULL,fastqInput1=NULL, fastqInput2=NULL, interleave = FALSE, paramList="default"){
    atacproc<-Bowtie2Mapping$new(atacProc=atacProc,bt2Idx=bt2Idx,samOutput=samOutput, fastqInput1=fastqInput1,
                                 fastqInput2=fastqInput2, interleave = interleave, paramList=paramList,reportPrefix=reportPrefix)
    atacproc$process()
    return(atacproc)
}
