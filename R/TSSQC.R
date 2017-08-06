TSSQC <-R6Class(
    classname = "TSSQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, txdb.knownGene = NULL,reportPrefix=NULL,bedInput = NULL,fregLenRange=c(0,2000),tssUpdownstream=1000,editable=FALSE){
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

            private$paramlist[["updownstream"]] <- tssUpdownstream
            private$paramlist[["fregLenRange"]] <- fregLenRange

            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
            private$checkRequireParam();
        },
        processing = function(){
            super$processing()
            genome <- Seqinfo(genome = NA_character_)
            readsbed <- unique(import(private$paramlist[["bedInput"]], genome = genome))

            readsbed<-readsbed[(width(readsbed)>=private$paramlist[["fregLenRange"]][1])&
                                               (width(readsbed)<=private$paramlist[["fregLenRange"]][2])]

            txdb<-private$paramlist[["knownGene"]]
            #trans<-GenomicFeatures::genes(txdb)#check gene tss or transcripts tss
            trans<-GenomicFeatures::transcripts(txdb)
            end(trans)<-start(trans)+1
            trans<-unique(trans)


            end(trans)<-start(trans)+private$paramlist[["updownstream"]]
            start(trans)<-start(trans)-private$paramlist[["updownstream"]]

            pairs<-findOverlapPairs(readsbed, trans,ignore.strand = TRUE)
            reads<-ranges(first(pairs))
            transspan<-second(pairs)
            #gp<-paste(as.character(start(transspan)),as.character(end(transspan)))
            #readsgps<-split(reads,gp)
            rvlist<-as.vector(strand(transspan)=="-")


            rs<-start(reads)-start(transspan)
            #rs[rs<0] <- 0
            rre<-2*private$paramlist[["updownstream"]]+1-rs[rvlist]
            re<-end(reads)-start(transspan)
            #re[re>2*private$paramlist[["updownstream"]]+1] <- 2*private$paramlist[["updownstream"]]+1
            rrs<-2*private$paramlist[["updownstream"]]+1-re[rvlist]
            rs[rvlist]<-rrs
            re[rvlist]<-rre
            re[re>2*private$paramlist[["updownstream"]]+1] <- 2*private$paramlist[["updownstream"]]+1
            rs[rs<0] <- 0
            start(reads)<-0
            end(reads)<-0
            end(reads)<-re
            start(reads)<- rs



            totaldistr<-as.numeric(coverage(reads,width = 2*private$paramlist[["updownstream"]]+1))

            df<-data.frame(val=totaldistr,x=-private$paramlist[["updownstream"]]:private$paramlist[["updownstream"]])
            ggplot(df,aes(x,val))+geom_line()+xlab("upstream<-TSS->downstream")+ylab("reads count")
            ggsave(paste(private$paramlist[["reportPrefix"]],".pdf"))

            pm<-GenomicFeatures::promoters(txdb)

            ############# drawing heatmap
            # gp<-mcols(transspan)$gene_id#tx_name
            # readsgps<-split(reads[width(reads)>100],gp)
            # heatmapdistr<-as.matrix(coverage(readsgps,width = 2*private$paramlist[["updownstream"]]+1))
            # idx<-order(sapply(1:nrow(heatmapdistr),function(i) max(heatmapdistr[i,])),decreasing = TRUE)
            # heatmapdistr<-heatmapdistr[idx,]
            #
            #
            #
            # dtst<-melt(t(heatmapdistr[100:20,]))
            # dtst$Var1<-dtst$Var1-private$paramlist[["updownstream"]]-1
            # plt<-ggplot(dtst,aes(Var1,Var2, fill=value))+geom_raster()+coord_cartesian(expand = FALSE)+
            #     theme(axis.ticks.y = element_blank(),panel.grid=element_blank(),axis.text.y = element_blank())+xlab("upstream<-TSS->downstream")+ylab("gene")
            #
            # ggsave("test.pdf")
            #################

            qcval=list();

            qcval[["totalUniqReads"]]<-length(readsbed)
            qcval[["TSSReads"]]<-length(unique(reads))
            qcval[["TSSRate"]]<-qcval[["TSSReads"]]/qcval[["totalUniqReads"]]

            exonlst<-exons(txdb)

            pairs<-findOverlapPairs(readsbed, exonlst,ignore.strand = TRUE)


            qcval[["exonReads"]]<-length(unique(first(pairs)))
            qcval[["exonRate"]]<-qcval[["exonReads"]]/qcval[["totalUniqReads"]]
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



        }
    )


)
