PeakQC <-R6Class(
    classname = "PeakQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, reportPrefix=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"),bedInput = NULL,editable=FALSE){
            super$initialize("PeakQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            qcbedInput <- qcbedInput[1]
            if(qcbedInput == "DHS"){
                private$paramlist[["qcbedInput"]]<-.obtainConfigure("DHS");
            }else if(qcbedInput == "blacklist"){
                private$paramlist[["qcbedInput"]]<-.obtainConfigure("blacklist");
            }else{
                private$paramlist[["qcbedInput"]]<-qcbedInput;
            }

            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
            }

            if(is.null(reportPrefix)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["reportPrefix"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
                #private$paramlist[["reportPrefix"]] <- paste0(private$paramlist[["bedInput"]],".BlacklistQCreport");
            }else{
                private$paramlist[["reportPrefix"]] <- reportPrefix;
            }
            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){
            genome <- Seqinfo(genome = .obtainConfigure("genome"))

            inputbed <- import(private$paramlist[["bedInput"]], genome = genome)


            qcbedInput<-import(private$paramlist[["qcbedInput"]], genome = genome)



            qcval=list();

            qcval[["totalInput"]]<-length(inputbed)
            qcval[["qcbedInput"]]<-length(subsetByOverlaps(inputbed, qcbedInput,ignore.strand = TRUE))
            qcval[["qcbedRate"]]<-qcval[["qcbedInput"]]/qcval[["totalInput"]]

            qcval<-as.matrix(qcval)
            write.table(qcval,file = paste0(private$paramlist[["reportPrefix"]],".txt"),sep="\t",quote = FALSE,col.names = FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }
            if(is.null(private$paramlist[["qcbedInput"]])){
                stop("qcbedInput is required.")
            }

        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileExist(private$paramlist[["qcbedInput"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
        }
    )


)


atacPeakQC<-function(atacProc, reportPrefix=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"), bedInput = NULL){
    atacproc<-PeakQC$new(atacProc, reportPrefix=reportPrefix,qcbedInput = qcbedInput,bedInput = bedInput,editable=FALSE)
    atacproc$process()
    return(atacproc)
}
