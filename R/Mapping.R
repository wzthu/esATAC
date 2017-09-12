Bowtie2Mapping <-R6Class(
    classname = "Bowtie2Mapping",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc,samOutput=NULL, bt2Idx=NULL, fastqInput1=NULL, fastqInput2=NULL, interleave = FALSE, paramList="default",reportOutput =NULL,editable=FALSE){
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

            
            if(is.null(paramList)){

            }else if(length(paramList)==1 && paramList=="default"){
                private$paramlist[["paramList"]]<-c("--no-discordant","--no-unal","--no-mixed","-X","2000",
                                                    private$paramlist[["paramList"]])
            }else{
                private$paramlist[["paramList"]]<-c(paramList,private$paramlist[["paramList"]])
                rejectp<-"-p|--threads|-x|-1|-2|-U|-S|--interleaved"
                private$checkParam(paramList,rejectp)
            }

            if(is.null(reportOutput)){
                if(!is.null(private$paramlist[["fastqInput1"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput1"]],regexProcName)
                    private$paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
            }else{
                private$paramlist[["reportOutput"]] <- reportOutput;
            }
            private$paramValidation()
        }


    ),
    private = list(
        processing = function(){
            if(.obtainConfigure("threads")>1){
                paramList<-paste(c(private$paramlist[["paramList"]],"-p",as.character(.obtainConfigure("threads"))),collapse = " ")
            }
            
            private$writeLog("start mapping with parameters:")
            private$writeLog(paste0("bowtie2 index:",private$paramlist[["bt2Idx"]]))
            private$writeLog(paste0("samOutput:",private$paramlist[["samOutput"]]))
            private$writeLog(paste0("report:",private$paramlist[["reportOutput"]]))
            private$writeLog(paste0("fastqInput1:",private$paramlist[["fastqInput1"]]))
            private$writeLog(paste0("fastqInput2:",private$paramlist[["fastqInput2"]]))
            private$writeLog(paste0("other parameters:",paste(private$paramlist[["paramList"]],collapse = " ")))
            if(length(paramList>0)){
                rs<-bowtie2(bt2Index = private$paramlist[["bt2Idx"]],
                        samOutput = private$paramlist[["samOutput"]],
                        seq1 = private$paramlist[["fastqInput1"]],
                        paramList,
                        seq2 = private$paramlist[["fastqInput2"]],
                        interleaved = private$paramlist[["interleave"]],
                        overwrite=TRUE)
            }else{
                rs<-bowtie2(bt2Index = private$paramlist[["bt2Idx"]],
                        samOutput = private$paramlist[["samOutput"]],
                        seq1 = private$paramlist[["fastqInput1"]],
                        seq2 = private$paramlist[["fastqInput2"]],
                        interleaved = private$paramlist[["interleave"]],
                        overwrite=TRUE)
            }
            
            writeLines(text = rs,con = private$paramlist[["reportOutput"]])



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
            private$checkFileCreatable(private$paramlist[["reportOutput"]]);
        },
        getReportValImp = function(item){
            if(sum(item == c("adapter1","adapter2"))>0){
                adapter<-readLines(paste0(private$paramlist[["reportOutput"]],item))
                return(adapter[1])
            }
            txt <- readLines(private$paramlist[["reportOutput"]])
            if(item == "total"){
                s<-strsplit(txt[1]," ")
                return(as.integer(s[[1]][1]))
            }
            if(item == "maprate"){
                s<-strsplit(txt[length(txt)],"% ") 
                return(as.numeric(s[[1]][1])/100)
            }
            if(item == "detail"){
                return(txt)
            }
            stop(paste0(item," is not an item of report value."))
        },
        getReportItemsImp = function(){
            return(c("total","maprate","detail"))
        }
    )


)


atacBowtie2Mapping <- function(atacProc,samOutput=NULL,reportOutput =NULL, bt2Idx=NULL,fastqInput1=NULL, fastqInput2=NULL, interleave = FALSE, paramList="default"){
    atacproc<-Bowtie2Mapping$new(atacProc=atacProc,bt2Idx=bt2Idx,samOutput=samOutput, fastqInput1=fastqInput1,
                                 fastqInput2=fastqInput2, interleave = interleave, paramList=paramList,reportOutput=reportOutput)
    atacproc$process()
    return(atacproc)
}
