FRiPQC <-R6Class(
    classname = "FRiPQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProcReads,atacProcPeak,reportPrefix=NULL,readsBedInput=NULL,peakBedInput=NULL,editable=FALSE){
            super$initialize("FRiPQC",editable,list(arg1=atacProcReads,arg2=atacProcPeak))
            if(!is.null(atacProcReads)){
                private$paramlist[["readsBedInput"]] <- atacProcReads$getParam("bedOutput");
                #private$paramlist[["readsCount"]] <- atacProcReads$getParam("readsCount")
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
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

            if(is.null(reportPrefix)){
                if(!is.null(private$paramlist[["peakBedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["peakBedInput"]],regexProcName)
                    private$paramlist[["reportPrefix"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
                #private$paramlist[["reportPrefix"]] <- paste0(private$paramlist[["peakBedInput"]],".FripQCreport");
            }else{
                private$paramlist[["reportPrefix"]] <- reportPrefix;
            }

            private$checkRequireParam();
        }
    ),
    private = list(
        processing = function(){
            qcval=list();
            genome <- Seqinfo(genome = NA_character_)
            gr_a <- import(private$paramlist[["readsBedInput"]], genome = genome)
            gr_b <- import(private$paramlist[["peakBedInput"]], genome = genome)
            qcval[["peakReads"]]<-length(subsetByOverlaps(gr_a, gr_b, ignore.strand = TRUE))
            if(is.null(private$paramlist[["readsCount"]])){
                qcval[["totalReads"]]<-R.utils::countLines(private$paramlist[["readsBedInput"]])
            }else{
                qcval[["totalReads"]]<-private$paramlist[["readsCount"]]
            }
            qcval[["FRiP"]]<-qcval[["peakReads"]]/qcval[["totalReads"]]
            #unlink(paste0(private$paramlist[["reportPrefix"]],".tmp"))
            ####private$paramlist[["qcval"]]<-qcval
            qcval<-as.matrix(qcval)
            write.table(qcval,file = paste0(private$paramlist[["reportPrefix"]],".txt"),sep="\t",quote = FALSE,col.names = FALSE)

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
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
        }
    )


)


atacFripQC<-function(atacProcReads,atacProcPeak,reportPrefix=NULL,readsBedInput=NULL,peakBedInput){
    fripQC<-FRiPQC$new(atacProcReads=atacProcReads,atacProcPeak=atacProcPeak,reportPrefix=reportPrefix,
                       readsBedInput=readsBedInput,peakBedInput=peakBedInput,editable=FALSE)
    fripQC$process()
    return(fripQC)
}
