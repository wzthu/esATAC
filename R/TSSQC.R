setClass(Class = "TSSQC",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "TSSQC",
    definition = function(.Object,prevSteps, ...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        txdbKnownGene <- allparam[["txdbKnownGene"]]
        bsgenome <- allparam[["bsgenome"]]
        reportPrefix <- allparam[["reportPrefix"]]
        fragLenRange <- allparam[["fragLenRange"]]
        tssUpdownstream <- allparam[["tssUpdownstream"]]
        
        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                atacProc <- prevSteps[[1]]
                atacProc<-c(unlist(atacProc),list())
                atacProc <- atacProc[[length(atacProc)]]
                input(.Object)[["bedInput"]] <- output(atacProc)[["bedOutput"]]
            }
        }
        
        if(!is.null(txdbKnownGene)){
            param(.Object)[["knownGene"]] <- txdbKnownGene;
        }else{
            param(.Object)[["knownGene"]]<- getRefRc("knownGene")
        }

        if(!is.null(bedInput)){
            input(.Object)[["bedInput"]] <- bedInput;
        }

        if(is.null(reportPrefix)){
            if(!is.null(input(.Object)[["bedInput"]])){
                output(.Object)[["tsspdfOutput"]] <- getAutoPath(.Object,input(.Object)[["bedInput"]],"BED|bed|Bed","pdf")
                output(.Object)[["tsstxtOutput"]] <- getAutoPath(.Object,input(.Object)[["bedInput"]],"BED|bed|Bed","txt")
                output(.Object)[["tssreportOutput"]] <- getAutoPath(.Object,input(.Object)[["bedInput"]],"BED|bed|Bed","report.txt")
            }
            #.Object@paramlist[["reportPrefix"]] <- paste0(.Object@paramlist[["bedInput"]],".TSSQCreport");
        }else{
            output(.Object)[["tsspdfOutput"]] <- paste0(reportPrefix,".pdf")
            output(.Object)[["tsstxtOutput"]] <- paste0(reportPrefix,".txt")
            output(.Object)[["tssreportOutput"]] <- paste0(reportPrefix,".report.txt")
        }

        param(.Object)[["updownstream"]] <- tssUpdownstream
        param(.Object)[["fragLenRange"]] <- fragLenRange

        param(.Object)[["bsgenome"]] <- bsgenome

        .Object
    }
)


setMethod(
    f = "processing",
    signature = "TSSQC",
    definition = function(.Object,...){
        if(is.null(param(.Object)[["bsgenome"]])){
            genome <- seqinfo(getRefRc("bsgenome"))
        }else{
            genome <- seqinfo(param(.Object)[["bsgenome"]])
        }
        #unique confilict with rJava, if solved, uncommented:
        #readsbed <- unique(import(.Object@paramlist[["bedInput"]], genome = genome))
#        readsbed <- import(.Object@paramlist[["bedInput"]], genome = genome)
       readsbed <- import(input(.Object)[["bedInput"]])
        readsbed<-readsbed[(width(readsbed)>=param(.Object)[["fragLenRange"]][1])&
                               (width(readsbed)<=param(.Object)[["fragLenRange"]][2])]

        txdb<-param(.Object)[["knownGene"]]
        if(is.character(txdb)){
            library(txdb,character.only = TRUE)
            txdb <- get0(txdb)
        }
        #trans<-GenomicFeatures::genes(txdb)#check gene tss or transcripts tss

        TSS <- promoters(txdb, upstream=param(.Object)[["updownstream"]], downstream=1+param(.Object)[["updownstream"]])
        #end(trans)<-start(trans)+1

        #unique confilict with rJava, if solved, uncommented:
        #TSS<-unique(TSS)


        #end(trans)<-start(trans)+param(.Object)[["updownstream"]]
        #start(trans)<-start(trans)-param(.Object)[["updownstream"]]

        pairs<-findOverlapPairs(readsbed, TSS,ignore.strand = TRUE)
        reads<-ranges(first(pairs))
        transspan<-second(pairs)
        #gp<-paste(as.character(start(transspan)),as.character(end(transspan)))
        #readsgps<-split(reads,gp)
        rvlist<-as.vector(strand(transspan)=="-")


        rs<-start(reads)-start(transspan)
        #rs[rs<0] <- 0
        rre<-2*param(.Object)[["updownstream"]]+1-rs[rvlist]
        re<-end(reads)-start(transspan)
        #re[re>2*param(.Object)[["updownstream"]]+1] <- 2*param(.Object)[["updownstream"]]+1
        rrs<-2*param(.Object)[["updownstream"]]+1-re[rvlist]
        rs[rvlist]<-rrs
        re[rvlist]<-rre
        re[re>2*param(.Object)[["updownstream"]]+1] <- 2*param(.Object)[["updownstream"]]+1
        rs[rs<0] <- 0
        start(reads)<-0
        end(reads)<-0
        end(reads)<-re
        start(reads)<- rs



        totaldistr<-as.numeric(coverage(reads,width = 2*param(.Object)[["updownstream"]]+1))

        df<-data.frame(counts=totaldistr,pos=-param(.Object)[["updownstream"]]:param(.Object)[["updownstream"]])
        ggplot(df,aes(pos,counts))+geom_line()+xlab("upstream<-TSS->downstream")+ylab("reads count")
        ggsave(output(.Object)[["tsspdfOutput"]])

        write.table(df,file = output(.Object)[["tsstxtOutput"]],sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)



        ############# drawing heatmap
        # gp<-mcols(transspan)$gene_id#tx_name
        # readsgps<-split(reads[width(reads)>100],gp)
        # heatmapdistr<-as.matrix(coverage(readsgps,width = 2*param(.Object)[["updownstream"]]+1))
        # idx<-order(sapply(1:nrow(heatmapdistr),function(i) max(heatmapdistr[i,])),decreasing = TRUE)
        # heatmapdistr<-heatmapdistr[idx,]
        #
        #
        #
        # dtst<-melt(t(heatmapdistr[100:20,]))
        # dtst$Var1<-dtst$Var1-param(.Object)[["updownstream"]]-1
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
        

        write.table(data.frame(qcval),file = output(.Object)[["tssreportOutput"]],sep="\t",quote = FALSE,col.names = FALSE)
        
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "TSSQC",
    definition = function(.Object, ...){
        tss <- read.table(file= output(.Object)[["tsstxtOutput"]],header=TRUE)
        report(.Object)$tss <- tss
        qcval <- as.list(read.table(output(.Object)[["tssreportOutput"]],header = TRUE,sep = "\t"))

        for(n in names(qcval)){
            report(.Object)[[n]] <- qcval[[n]]
        }  
        .Object
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
#' @param newStepType \code{Character} scalar.
#' New class name
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
#' setTmpDir(td)
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
                                fragLenRange=c(0,2000),tssUpdownstream=1000, 
                                newStepType = "TSSQC", ...) standardGeneric("atacTSSQC"))

#' @rdname TSSQC
#' @aliases atacTSSQC
#' @export
setMethod(
    f = "atacTSSQC",
    signature = "ATACProc",
    definition = function(atacProc, txdbKnownGene = NULL,bsgenome = NULL,
                          reportPrefix=NULL,bedInput = NULL,
                          fragLenRange=c(0,2000),tssUpdownstream=1000, newStepType = "TSSQC", ...){
        allpara <- c(list(Class = regAttachedStep(newStepType,"TSSQC"), prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)




#' @rdname TSSQC
#' @aliases tssQC
#' @export
tssQC<-function(bedInput, txdbKnownGene = NULL,
                bsgenome = NULL,reportPrefix=NULL,
                fragLenRange=c(0,2000),tssUpdownstream=1000, newStepType = "TSSQC", ...){
    allpara <- c(list(Class = regAttachedStep(newStepType,"TSSQC"), prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}




