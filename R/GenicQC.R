GenicQC <-R6Class(
    classname = "GenicQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, txdb.knownGene = NULL,reportPrefix=NULL,bedInput = NULL,promoterRange=c(-2000,2000), editable=FALSE){
            super$initialize("TSSQC",editable,list(arg1=atacProcReads,arg2=atacProcPeak))
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
                private$paramlist[["reportPrefix"]] <- paste(private$paramlist[["bedInput"]],".TSSQCreport");
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
            genome <- Seqinfo(genome = NA_character_)
            readsbed <- unique(import(private$paramlist[["bedInput"]], genome = genome))



            txdb<-private$paramlist[["knownGene"]]
            #trans<-GenomicFeatures::genes(txdb)#check gene tss or transcripts tss
            trans<-GenomicFeatures::transcripts(txdb)


            qcval=list();

            qcval[["totalUniqReadsOrPeak"]]<-length(readsbed)

            genedatalist<-exons(txdb)
            pairs<-findOverlapPairs(readsbed, exonlst,ignore.strand = TRUE)
            qcval[["exon"]]<-length(unique(first(pairs)))
            qcval[["exonRate"]]<-qcval[["exonReads"]]/qcval[["totalUniqReadsOrPeak"]]

            genedatalist<-genes(txdb)
            pairs<-findOverlapPairs(readsbed, exonlst,ignore.strand = TRUE)
            qcval[["gene"]]<-length(unique(first(pairs)))
            qcval[["geneRate"]]<-qcval[["exonReads"]]/qcval[["totalUniqReadsOrPeak"]]

            genedatalist<-promoters(txdb,upstream = abs(private$paramlist[["promoterRange"]][1]),downstream = abs(private$paramlist[["promoterRange"]][2]))
            pairs<-findOverlapPairs(readsbed, genedatalist,ignore.strand = TRUE)
            qcval[["promoter"]]<-length(unique(first(pairs)))
            qcval[["promoterRate"]]<-qcval[["exonReads"]]/qcval[["totalUniqReadsOrPeak"]]


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
            if(is.null(private$paramlist[["txdb.knownGene"]])){
                stop("txdb.knownGene is required.")
            }
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }

            if(is.null(private$paramlist[["fregLenRange"]])){
                stop("fregLenRange is required.")
            }


        }
    )


)
