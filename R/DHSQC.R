DHSQC <-R6Class(
    classname = "DHSQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, reportPrefix=NULL,bedDHS = NULL,bedInput = NULL,editable=FALSE){
            super$initialize("DHSQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            if(!is.null(bedDHS)){
                private$paramlist[["bedDHS"]] <- bedDHS;
            }else{
                private$paramlist[["bedDHS"]]<-.obtainConfigure("DHS");
            }

            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
            }

            if(is.null(reportPrefix)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["reportPrefix"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
                #private$paramlist[["reportPrefix"]] <- paste0(private$paramlist[["bedInput"]],".DHSQCreport");
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


            dhsbed<-import(private$paramlist[["bedDHS"]], genome = genome)


            qcval=list();

            qcval[["totalInput"]]<-length(inputbed)
            qcval[["dhsInput"]]<-length(subsetByOverlaps(inputbed, dhsbed,ignore.strand = TRUE))
            qcval[["dhsRate"]]<-qcval[["dhsInput"]]/qcval[["totalInput"]]

            qcval<-as.matrix(qcval)
            write.table(qcval,file = paste0(private$paramlist[["reportPrefix"]],".txt"),sep="\t",quote = FALSE,col.names = FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }
            if(is.null(private$paramlist[["bedDHS"]])){
                stop("bedDHS is required.")
            }

        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileExist(private$paramlist[["bedDHS"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
        }
    )


)

atacDHSQC<-function(atacProc, reportPrefix=NULL,bedDHS = NULL,bedInput = NULL){
    atacproc<-DHSQC$new(atacProc, reportPrefix=reportPrefix,bedDHS = bedDHS,bedInput = bedInput)
    atacproc$process()
    return(atacproc)
}
