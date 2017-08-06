DHSQC <-R6Class(
    classname = "DHSQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, reportPrefix=NULL,bedDHS = NULL,bedInput = NULL,editable=FALSE){
            super$initialize("DHSQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
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
                private$paramlist[["reportPrefix"]] <- paste(private$paramlist[["bedInput"]],".DHSQCreport");
            }else{
                private$paramlist[["reportPrefix"]] <- reportPrefix;
            }



            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileExist(private$paramlist[["bedDHS"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
            private$checkRequireParam();
        },
        processing = function(){
            super$processing()
            genome <- Seqinfo(genome = NA_character_)

            inputbed <- unique(import(private$paramlist[["bedInput"]], genome = genome))


            dhsbed<-import(private$paramlist[["bedDHS"]], genome = genome)

            pairs<-findOverlapPairs(inputbed, dhsbed,ignore.strand = TRUE)
            ibed<-unique(ranges(first(pairs)))



            qcval=list();

            qcval[["totalUniqInput"]]<-length(inputbed)
            qcval[["dhsInput"]]<-length(ibed)
            qcval[["dhsRate"]]<-qcval[["dhsInput"]]/qcval[["totalUniqInput"]]

            qcval<-as.matrix(qcval)
            write.table(qcval,file = paste0(private$paramlist[["reportPrefix"]],".txt"),sep="\t",quote = FALSE,col.names = FALSE)
            private$finish <- TRUE
        },
        setResultParam = function(fastqOutput1, fastqOutput2=NULL){
            super$setResultParam();
            private$paramlist[["fastqOutput1"]] <- fastqOutput1
            private$paramlist[["fastqOutput2"]] <- fastqOutput2
        }
    ),
    private = list(
        checkRequireParam = function(){
            if(private$editable){
                return();
            }
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }
            if(is.null(private$paramlist[["bedDHS"]])){
                stop("bedDHS is required.")
            }

        }
    )


)
