BlacklistQC <-R6Class(
    classname = "BlacklistQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, reportPrefix=NULL,bedBlacklist = NULL,bedInput = NULL,editable=FALSE){
            super$initialize("BlacklistQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            if(!is.null(bedBlacklist)){
                private$paramlist[["bedBlacklist"]] <- bedBlacklist;
            }else{
                private$paramlist[["bedBlacklist"]]<-.obtainConfigure("blacklist");
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


            blacklistbed<-import(private$paramlist[["bedBlacklist"]], genome = genome)



            qcval=list();

            qcval[["totalInput"]]<-length(inputbed)
            qcval[["blacklistInput"]]<-length(subsetByOverlaps(inputbed, blacklistbed,ignore.strand = TRUE))
            qcval[["blacklistRate"]]<-qcval[["blacklistInput"]]/qcval[["totalInput"]]

            qcval<-as.matrix(qcval)
            write.table(qcval,file = paste0(private$paramlist[["reportPrefix"]],".txt"),sep="\t",quote = FALSE,col.names = FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }
            if(is.null(private$paramlist[["bedBlacklist"]])){
                stop("bedBlacklist is required.")
            }

        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileExist(private$paramlist[["bedBlacklist"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
        }
    )


)


atacBlacklistQC<-function(atacProc, reportPrefix=NULL,bedBlacklist = NULL,bedInput = NULL){
    atacproc<-BlacklistQC$new(atacProc, reportPrefix=reportPrefix,bedBlacklist = bedBlacklist,bedInput = bedInput,editable=FALSE)
    atacproc$process()
    return(atacproc)
}
