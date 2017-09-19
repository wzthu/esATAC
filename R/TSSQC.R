TSSQC <-R6Class(
    classname = "TSSQC",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, txdbKnownGene = NULL,reportPrefix=NULL,bedInput = NULL,fregLenRange=c(0,2000),tssUpdownstream=1000,editable=FALSE){
            super$initialize("TSSQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            if(!is.null(txdbKnownGene)){
                private$paramlist[["knownGene"]] <- txdbKnownGene;
            }else{
                private$paramlist[["knownGene"]]<-.obtainConfigure("knownGene");
            }

            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
            }

            if(is.null(reportPrefix)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["tsspdfOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".pdf"))
                    private$paramlist[["tsstxtOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".txt"))
                    private$paramlist[["tssreportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report.txt"))
                }
                #private$paramlist[["reportPrefix"]] <- paste0(private$paramlist[["bedInput"]],".TSSQCreport");
            }else{
                private$paramlist[["tsspdfOutput"]] <- paste0(reportPrefix,".pdf")
                private$paramlist[["tsstxtOutput"]] <- paste0(reportPrefix,".txt")
                private$paramlist[["tssreportOutput"]] <- paste0(reportPrefix,".report.txt")
            }

            private$paramlist[["updownstream"]] <- tssUpdownstream
            private$paramlist[["fregLenRange"]] <- fregLenRange


            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){

            genome <- Seqinfo(genome = .obtainConfigure("genome"))
            readsbed <- unique(import(private$paramlist[["bedInput"]], genome = genome))

            readsbed<-readsbed[(width(readsbed)>=private$paramlist[["fregLenRange"]][1])&
                                   (width(readsbed)<=private$paramlist[["fregLenRange"]][2])]

            txdb<-private$paramlist[["knownGene"]]
            #trans<-GenomicFeatures::genes(txdb)#check gene tss or transcripts tss
            TSS <- promoters(txdb, upstream=private$paramlist[["updownstream"]], downstream=1+private$paramlist[["updownstream"]])
            #end(trans)<-start(trans)+1
            TSS<-unique(TSS)


            #end(trans)<-start(trans)+private$paramlist[["updownstream"]]
            #start(trans)<-start(trans)-private$paramlist[["updownstream"]]

            pairs<-findOverlapPairs(readsbed, TSS,ignore.strand = TRUE)
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

            df<-data.frame(counts=totaldistr,pos=-private$paramlist[["updownstream"]]:private$paramlist[["updownstream"]])
            ggplot(df,aes(pos,counts))+geom_line()+xlab("upstream<-TSS->downstream")+ylab("reads count")
            ggsave(private$paramlist[["tsspdfOutput"]])
            
            write.table(df,file = private$paramlist[["tsstxtOutput"]],sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)



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
            print(qcval[["totalUniqReads"]])
            qcval[["TSSReads"]]<-length(subsetByOverlaps(readsbed, TSS,ignore.strand = TRUE))
            print(qcval[["TSSReads"]])
            qcval[["TSSRate"]]<-qcval[["TSSReads"]]/qcval[["totalUniqReads"]]
            print(qcval[["TSSRate"]])
            qcval<-as.matrix(qcval)
            
            write.table(qcval,file = private$paramlist[["tssreportOutput"]],sep="\t",quote = FALSE,col.names = FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["knownGene"]])){
                stop("txdbKnownGene is required.")
            }
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileCreatable(private$paramlist[["tsspdfOutput"]]);
            private$checkFileCreatable(private$paramlist[["tsstxtOutput"]]);
            private$checkFileCreatable(private$paramlist[["tssreportOutput"]]);
        },
        getReportValImp = function(item){
            tss <- read.table(file= private$paramlist[["tsstxtOutput"]],header=TRUE)
            if(item == "tss"){
                return(tss)
            }
        },
        getReportItemsImp = function(){
            return(c("tss"))
        }
    )


)

#' @name atacTSSQC
#' @aliases atacTSSQC
#' @aliases tssQC
#' @title Quality control for transcription start site(TSS) reads enrichment
#' @description 
#' These functions are used to generate the reads coverage plot around TSS.  
#' @param atacProc \code{\link{ATACProc}} object scalar. 
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}, 
#' \code{\link{atacBedUtils}}.
#' @param txdbKnownGene \code{TxDb} object scalar. 
#' TxDb object for specific species. 
#' @param reportPrefix \code{Character} scalar. 
#' The prefix of report files path. 
#' @param bedInput \code{Character} scalar. 
#' BED file input path.
#' @param fregLenRange \code{Interger} vector of 2 element. 
#' The fregment length ranges.
#' @param tssUpdownstream \code{Interger} scalar.
#' The upstream and downstrem from TSS locations.
#' @details The parameter related to input and output file path
#' will be automatically 
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or 
#' generated based on known parameters 
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently, 
#' \code{atacProc} should be set \code{NULL} 
#' or you can use \code{tssQC} instead.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacSamToBed}} 
#' \code{\link{atacBedUtils}}

#' @rdname atacTSSQC
#' @export 
atacTSSQC<-function(atacProc, txdbKnownGene = NULL,reportPrefix=NULL,bedInput = NULL,fregLenRange=c(0,2000),tssUpdownstream=1000){
    tssQC<-TSSQC$new(atacProc=atacProc, txdbKnownGene=txdbKnownGene,reportPrefix=reportPrefix,bedInput=bedInput,fregLenRange=fregLenRange,tssUpdownstream=tssUpdownstream,editable=FALSE)
    tssQC$process()
    invisible(tssQC)
}
#' @rdname atacTSSQC
#' @export 
tssQC<-function(bedInput, txdbKnownGene = NULL,reportPrefix=NULL,fregLenRange=c(0,2000),tssUpdownstream=1000){
    tssQC<-TSSQC$new(atacProc=NULL, txdbKnownGene=txdbKnownGene,reportPrefix=reportPrefix,bedInput=bedInput,fregLenRange=fregLenRange,tssUpdownstream=tssUpdownstream,editable=FALSE)
    tssQC$process()
    invisible(tssQC)
}




