setClass(Class = "TSSQC",
         contains = "ATACProc"
)

setMethod(
    f = "initialize",
    signature = "TSSQC",
    definition = function(.Object,atacProc, ..., txdbKnownGene = NULL,
                          bsgenome = NULL,reportPrefix=NULL,bedInput = NULL,
                          fregLenRange=c(0,2000),tssUpdownstream=1000,editable=FALSE){
        .Object <- init(.Object,"TSSQC",editable,list(arg1=atacProc))
        if(!is.null(atacProc)){
            .Object@paramlist[["bedInput"]] <- getParam(atacProc, "bedOutput");
            regexProcName<-sprintf("(BED|bed|Bed|%s)", getProcName(atacProc))
        }else{
            regexProcName<-"(BED|bed|Bed)"
        }
        if(!is.null(txdbKnownGene)){
            .Object@paramlist[["knownGene"]] <- txdbKnownGene;
        }else{
            .Object@paramlist[["knownGene"]]<-.obtainConfigure("knownGene");
        }

        if(!is.null(bedInput)){
            .Object@paramlist[["bedInput"]] <- bedInput;
        }

        if(is.null(reportPrefix)){
            if(!is.null(.Object@paramlist[["bedInput"]])){
                prefix<-getBasenamePrefix(.Object,.Object@paramlist[["bedInput"]],regexProcName)
                .Object@paramlist[["tsspdfOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".pdf"))
                .Object@paramlist[["tsstxtOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".txt"))
                .Object@paramlist[["tssreportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".report.txt"))
            }
            #.Object@paramlist[["reportPrefix"]] <- paste0(.Object@paramlist[["bedInput"]],".TSSQCreport");
        }else{
            .Object@paramlist[["tsspdfOutput"]] <- paste0(reportPrefix,".pdf")
            .Object@paramlist[["tsstxtOutput"]] <- paste0(reportPrefix,".txt")
            .Object@paramlist[["tssreportOutput"]] <- paste0(reportPrefix,".report.txt")
        }

        .Object@paramlist[["updownstream"]] <- tssUpdownstream
        .Object@paramlist[["fregLenRange"]] <- fregLenRange

        .Object@paramlist[["bsgenome"]] <- bsgenome

        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "TSSQC",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bsgenome"]])){
            genome <- seqinfo(.obtainConfigure("bsgenome"))
        }else{
            genome <- seqinfo(.Object@paramlist[["bsgenome"]])
        }
        #unique confilict with rJava, if solved, uncommented:
        #readsbed <- unique(import(.Object@paramlist[["bedInput"]], genome = genome))
#        readsbed <- import(.Object@paramlist[["bedInput"]], genome = genome)
       readsbed <- import(.Object@paramlist[["bedInput"]])
        readsbed<-readsbed[(width(readsbed)>=.Object@paramlist[["fregLenRange"]][1])&
                               (width(readsbed)<=.Object@paramlist[["fregLenRange"]][2])]

        txdb<-.Object@paramlist[["knownGene"]]
        #trans<-GenomicFeatures::genes(txdb)#check gene tss or transcripts tss

        TSS <- promoters(txdb, upstream=.Object@paramlist[["updownstream"]], downstream=1+.Object@paramlist[["updownstream"]])
        #end(trans)<-start(trans)+1

        #unique confilict with rJava, if solved, uncommented:
        #TSS<-unique(TSS)


        #end(trans)<-start(trans)+.Object@paramlist[["updownstream"]]
        #start(trans)<-start(trans)-.Object@paramlist[["updownstream"]]

        pairs<-findOverlapPairs(readsbed, TSS,ignore.strand = TRUE)
        reads<-ranges(first(pairs))
        transspan<-second(pairs)
        #gp<-paste(as.character(start(transspan)),as.character(end(transspan)))
        #readsgps<-split(reads,gp)
        rvlist<-as.vector(strand(transspan)=="-")


        rs<-start(reads)-start(transspan)
        #rs[rs<0] <- 0
        rre<-2*.Object@paramlist[["updownstream"]]+1-rs[rvlist]
        re<-end(reads)-start(transspan)
        #re[re>2*.Object@paramlist[["updownstream"]]+1] <- 2*.Object@paramlist[["updownstream"]]+1
        rrs<-2*.Object@paramlist[["updownstream"]]+1-re[rvlist]
        rs[rvlist]<-rrs
        re[rvlist]<-rre
        re[re>2*.Object@paramlist[["updownstream"]]+1] <- 2*.Object@paramlist[["updownstream"]]+1
        rs[rs<0] <- 0
        start(reads)<-0
        end(reads)<-0
        end(reads)<-re
        start(reads)<- rs



        totaldistr<-as.numeric(coverage(reads,width = 2*.Object@paramlist[["updownstream"]]+1))

        df<-data.frame(counts=totaldistr,pos=-.Object@paramlist[["updownstream"]]:.Object@paramlist[["updownstream"]])
        ggplot(df,aes(pos,counts))+geom_line()+xlab("upstream<-TSS->downstream")+ylab("reads count")
        ggsave(.Object@paramlist[["tsspdfOutput"]])

        write.table(df,file = .Object@paramlist[["tsstxtOutput"]],sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)



        ############# drawing heatmap
        # gp<-mcols(transspan)$gene_id#tx_name
        # readsgps<-split(reads[width(reads)>100],gp)
        # heatmapdistr<-as.matrix(coverage(readsgps,width = 2*.Object@paramlist[["updownstream"]]+1))
        # idx<-order(sapply(1:nrow(heatmapdistr),function(i) max(heatmapdistr[i,])),decreasing = TRUE)
        # heatmapdistr<-heatmapdistr[idx,]
        #
        #
        #
        # dtst<-melt(t(heatmapdistr[100:20,]))
        # dtst$Var1<-dtst$Var1-.Object@paramlist[["updownstream"]]-1
        # plt<-ggplot(dtst,aes(Var1,Var2, fill=value))+geom_raster()+coord_cartesian(expand = FALSE)+
        #     theme(axis.ticks.y = element_blank(),panel.grid=element_blank(),axis.text.y = element_blank())+xlab("upstream<-TSS->downstream")+ylab("gene")
        #
        # ggsave("test.pdf")
        #################

        qcval=list();

        qcval[["totalUniqReads"]]<-length(readsbed)
        writeLog(.Object,sprintf("Total Unique Reads: %.0f",qcval[["totalUniqReads"]]))
        qcval[["TSSReads"]]<-length(subsetByOverlaps(readsbed, TSS,ignore.strand = TRUE))
        writeLog(.Object,sprintf("TSS Reads: %.0f",qcval[["TSSReads"]]))
        qcval[["TSSRate"]]<-qcval[["TSSReads"]]/qcval[["totalUniqReads"]]
        writeLog(.Object,sprintf("TSS Rate: %f",qcval[["TSSRate"]]))
        qcval<-as.matrix(qcval)

        write.table(qcval,file = .Object@paramlist[["tssreportOutput"]],sep="\t",quote = FALSE,col.names = FALSE)
        .Object
    }
)



setMethod(
    f = "checkRequireParam",
    signature = "TSSQC",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["knownGene"]])){
            stop("txdbKnownGene is required.")
        }
        if(is.null(.Object@paramlist[["bedInput"]])){
            stop("bedInput is required.")
        }
    }
)

setMethod(
    f = "checkAllPath",
    signature = "TSSQC",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["bedInput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["tsspdfOutput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["tsstxtOutput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["tssreportOutput"]]);
    }
)


setMethod(
    f = "getReportValImp",
    signature = "TSSQC",
    definition = function(.Object, item,...){
        tss <- read.table(file= .Object@paramlist[["tsstxtOutput"]],header=TRUE)
        if(item == "tss"){
            return(tss)
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "TSSQC",
    definition = function(.Object){
        return(c("tss"))
    }
)



#' @name TSSQC
#' @title Quality control for transcription start site(TSS) reads enrichment
#' @description
#' These functions are used to generate the reads coverage plot around TSS.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}},
#' \code{\link{atacBedUtils}}.
#' @param txdbKnownGene \code{TxDb} object scalar.
#' TxDb object for specific species.
#' @param bsgenome \code{BSGenome} object scalar.
#' BSGenome object for specific species.
#' @param reportPrefix \code{Character} scalar.
#' The prefix of report files path.
#' @param bedInput \code{Character} scalar.
#' BED file input path.
#' @param fragLenRange \code{Interger} vector of 2 element.
#' The fragment length ranges.
#' @param tssUpdownstream \code{Interger} scalar.
#' The upstream and downstrem from TSS locations.
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{tssQC} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' \code{\link{atacBedUtils}}
#' \code{\link{bedUtils}}
#'
#' @examples
#' library(R.utils)
#' td <- tempdir()
#' options(atacConf=setConfigure("tmpdir",td))
#'
#' bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' tssQC(bedfile,TxDb.Hsapiens.UCSC.hg19.knownGene,BSgenome.Hsapiens.UCSC.hg19,fragLenRange=c(180,247))
#'
#' dir(td)
#' @importFrom BiocGenerics strand


setGeneric("atacTSSQC",function(atacProc, txdbKnownGene = NULL,bsgenome = NULL,
                                  reportPrefix=NULL,bedInput = NULL,
                                  fragLenRange=c(0,2000),tssUpdownstream=1000, ...) standardGeneric("atacTSSQC"))

#' @rdname TSSQC
#' @aliases atacTSSQC
#' @export
setMethod(
    f = "atacTSSQC",
    signature = "ATACProc",
    definition = function(atacProc, txdbKnownGene = NULL,bsgenome = NULL,
                          reportPrefix=NULL,bedInput = NULL,
                          fragLenRange=c(0,2000),tssUpdownstream=1000, ...){
        atacproc <- new(
            "TSSQC",
            atacProc = atacProc,
            txdbKnownGene = txdbKnownGene,
            bsgenome = bsgenome,
            reportPrefix = reportPrefix,
            bedInput = bedInput,
            fregLenRange = fragLenRange,
            tssUpdownstream = tssUpdownstream)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)




#' @rdname TSSQC
#' @aliases tssQC
#' @export
tssQC<-function(bedInput, txdbKnownGene = NULL,bsgenome = NULL,reportPrefix=NULL,fragLenRange=c(0,2000),tssUpdownstream=1000, ...){
    atacproc <- new(
        "TSSQC",
        atacProc = NULL,
        txdbKnownGene = txdbKnownGene,
        bsgenome = bsgenome,
        reportPrefix = reportPrefix,
        bedInput = bedInput,
        fregLenRange = fragLenRange,
        tssUpdownstream = tssUpdownstream)
    atacproc <- process(atacproc)
    invisible(atacproc)
}




