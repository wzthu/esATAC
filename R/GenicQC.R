GenicQC <-R6Class(
    classname = "GenicQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, txdb.knownGene = NULL,reportPrefix=NULL,bedInput = NULL,promoterRange=c(-2000,2000), editable=FALSE){
            super$initialize("GenicQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            if(!is.null(txdb.knownGene)){
                private$paramlist[["knownGene"]] <- txdb.knownGene;
            }else{
                private$paramlist[["knownGene"]]<-.obtainConfigure("knownGene");
            }

            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
            }

            if(is.null(reportPrefix)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["reportPrefix"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
                #private$paramlist[["reportPrefix"]] <- paste(private$paramlist[["bedInput"]],".GenicQCreport");
            }else{
                private$paramlist[["reportPrefix"]] <- reportPrefix;
            }

            private$paramlist[["promoterRange"]] <- promoterRange



            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){
            genome <- Seqinfo(genome = .obtainConfigure("genome"))
            readsbed <- unique(import(private$paramlist[["bedInput"]], genome = genome))



            txdb<-private$paramlist[["knownGene"]]



            qcval=list();

            qcval[["totalUniqReadsOrPeak"]]<-length(readsbed)

            genedatalist<-exons(txdb)
            qcval[["exon"]]<-length(subsetByOverlaps(readsbed, genedatalist,ignore.strand = TRUE))
            qcval[["exonRate"]]<-qcval[["exon"]]/qcval[["totalUniqReadsOrPeak"]]

            genedatalist<-genes(txdb)
            qcval[["gene"]]<-length(subBed<-subsetByOverlaps(readsbed, genedatalist,ignore.strand = TRUE))
            qcval[["geneRate"]]<-qcval[["gene"]]/qcval[["totalUniqReadsOrPeak"]]

            genedatalist<-promoters(txdb,upstream = abs(private$paramlist[["promoterRange"]][1]),downstream = abs(private$paramlist[["promoterRange"]][2]))
            qcval[["promoter"]]<-length(subsetByOverlaps(readsbed, genedatalist,ignore.strand = TRUE))
            qcval[["promoterRate"]]<-qcval[["promoter"]]/qcval[["totalUniqReadsOrPeak"]]

            qcval[["intron"]]<-qcval[["gene"]]-qcval[["exon"]]
            qcval[["intronRate"]]<-qcval[["intron"]]/qcval[["totalUniqReadsOrPeak"]]

            qcval[["intergenic"]]<-qcval[["totalUniqReadsOrPeak"]]-qcval[["gene"]]-qcval[["promoter"]]
            qcval[["intergenicRate"]]<-qcval[["intergenic"]]/qcval[["totalUniqReadsOrPeak"]]

            qcval<-as.matrix(qcval)
            write.table(qcval,file = paste0(private$paramlist[["reportPrefix"]],".txt"),sep="\t",quote = FALSE,col.names = FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["knownGene"]])){
                stop("txdb.knownGene is required.")
            }
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }
            if(is.null(private$paramlist[["promoterRange"]])){
                stop("promoterRange is required.")
            }

        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
        }
    )


)

atacGenicQC<-function(atacProc, txdb.knownGene = NULL,reportPrefix=NULL,bedInput = NULL,promoterRange=c(-2000,2000)){
    atacproc<-GenicQC(atacProc, txdb.knownGene = txdb.knownGene,reportPrefix=reportPrefix,bedInput = bedInput,promoterRange=promoterRange)
    atacproc$process()
    return(atacproc)
}
