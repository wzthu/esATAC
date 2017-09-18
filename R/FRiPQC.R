FRiPQC <-R6Class(
    classname = "FRiPQC",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProcReads,atacProcPeak,reportOutput=NULL,readsBedInput=NULL,peakBedInput=NULL,editable=FALSE){
            super$initialize("FRiPQC",editable,list(arg1=atacProcReads,arg2=atacProcPeak))
            if(!is.null(atacProcReads)){
                private$paramlist[["readsBedInput"]] <- atacProcReads$getParam("bedOutput");
                #private$paramlist[["readsCount"]] <- atacProcReads$getParam("readsCount")
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProcReads$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            if(!is.null(atacProcPeak)){
                private$paramlist[["peakBedInput"]] <- atacProcPeak$getParam("bedOutput");
            }

            if(!is.null(readsBedInput)){
                private$paramlist[["readsBedInput"]] <- readsBedInput
            }
            if(!is.null(peakBedInput)){
                private$paramlist[["peakBedInput"]] <- peakBedInput
            }

            if(is.null(reportOutput)){
                if(!is.null(private$paramlist[["peakBedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["peakBedInput"]],regexProcName)
                    private$paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report.txt"))
                }
            }else{
                private$paramlist[["reportOutput"]] <- reportOutput;
            }

            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){
            qcval=list();
            genome <- Seqinfo(genome = .obtainConfigure("genome"))
            gr_a <- import(private$paramlist[["readsBedInput"]], genome = genome)
            gr_b <- import(private$paramlist[["peakBedInput"]], genome = genome)
            qcval[["peakReads"]]<-length(subsetByOverlaps(gr_a, gr_b, ignore.strand = TRUE))
            qcval[["totalReads"]]<-R.utils::countLines(private$paramlist[["readsBedInput"]])
            qcval[["totalPeaks"]]<-length(gr_b)
            qcval[["FRiP"]]<-qcval[["peakReads"]]/qcval[["totalReads"]]
            #unlink(paste0(private$paramlist[["reportPrefix"]],".tmp"))
            ####private$paramlist[["qcval"]]<-qcval
            write.table(as.data.frame(qcval),file = private$paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["readsBedInput"]])){
                stop("readsBedInput is required.")
            }
            if(is.null(private$paramlist[["peakBedInput"]])){
                stop("peakBedInput is required.")
            }

        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["readsBedInput"]]);
            private$checkFileExist(private$paramlist[["peakBedInput"]]);
            private$checkFileCreatable(private$paramlist[["reportOutput"]]);
        },
        getReportValImp = function(item){
            qcval <- as.list(read.table(file= private$paramlist[["reportOutput"]],header=TRUE))
            if(item == "report"){
                return(data.frame(Item=names(qcval),Value=as.character(qcval)))
            }else{
                return(qcval[[item]])
            }
        },
        getReportItemsImp = function(){
            return(c("report","peakReads","totalReads","totalPeaks","FRiP"))
        }
    )


)


atacFripQC<-function(atacProcReads,atacProcPeak,reportOutput=NULL,readsBedInput=NULL,peakBedInput=NULL){
    fripQC<-FRiPQC$new(atacProcReads=atacProcReads,atacProcPeak=atacProcPeak,reportOutput=reportOutput,
                       readsBedInput=readsBedInput,peakBedInput=peakBedInput,editable=FALSE)
    fripQC$process()
    invisible(fripQC)
}
