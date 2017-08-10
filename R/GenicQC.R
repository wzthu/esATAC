GenicQC <-R6Class(
    classname = "GenicQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, txdb.knownGene = NULL,reportPrefix=NULL,bedInput = NULL,promoterRange=c(-2000,2000), editable=FALSE){
            super$initialize("GenicQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
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
                private$paramlist[["reportPrefix"]] <- paste(private$paramlist[["bedInput"]],".GenicQCreport");
            }else{
                private$paramlist[["reportPrefix"]] <- reportPrefix;
            }

            private$paramlist[["promoterRange"]] <- promoterRange


            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
            private$checkRequireParam();
        },
        processing = function(){
            super$processing()
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
            if(is.null(private$paramlist[["knownGene"]])){
                stop("txdb.knownGene is required.")
            }
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }

            if(is.null(private$paramlist[["promoterRange"]])){
                stop("promoterRange is required.")
            }


        }
    )


)
