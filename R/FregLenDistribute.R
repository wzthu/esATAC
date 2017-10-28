setClass(Class = "FregLenDistr",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "FregLenDistr",
    definition = function(.Object,atacProc,...,reportPrefix=NULL,bedInput=NULL,editable=FALSE){
        .Object <- init(.Object,"FregLenDistr",editable,list(arg1=atacProc))
        if(.Object@singleEnd){
            .Object <- writeLog(.Object,"This process is for pair-end sequencing data.",isWarnning=TRUE)
        }
        if(!is.null(atacProc)){
            .Object@paramlist[["bedInput"]] <- getParam(atacProc, "bedOutput");
            regexProcName<-sprintf("(BED|bed|Bed|%s)",getProcName(atacProc))
        }else{
            regexProcName<-"(BED|bed|Bed)"
        }
        if(!is.null(bedInput)){
            .Object@paramlist[["bedInput"]] <- bedInput;
        }
        if(is.null(reportPrefix)){
            if(!is.null(.Object@paramlist[["bedInput"]])){
                prefix<-getBasenamePrefix(.Object,.Object@paramlist[["bedInput"]],regexProcName)
                .Object@paramlist[["lendistrpdfOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".lendistr.pdf"))
                .Object@paramlist[["lendistrtxtOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".lendistr.txt"))
                .Object@paramlist[["dnagroovepdfOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".dnagroove.pdf"))
                .Object@paramlist[["histonepdfOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".histone.pdf"))
            }
        }else{
            .Object@paramlist[["lendistrpdfOutput"]] <- paste0(reportPrefix,".lendistr.pdf")
            .Object@paramlist[["lendistrtxtOutput"]] <- paste0(reportPrefix,".lendistr.txt")
            .Object@paramlist[["dnagroovepdfOutput"]] <- paste0(reportPrefix,".dnagroove.pdf")
            .Object@paramlist[["histonepdfOutput"]] <- paste0(reportPrefix,".histone.pdf")
        }
        paramValidation(.Object)
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "FregLenDistr",
    definition = function(.Object,...){
        readslist<-read.table(file = .Object@paramlist[["bedInput"]],nrows = 1)
        bedcol=length(colnames(readslist))
        if(bedcol>3){
            readslist<-read.table(file = .Object@paramlist[["bedInput"]],colClasses = c("NULL","integer","integer",rep("NULL",bedcol-3)))
        }else{
            readslist<-read.table(file = .Object@paramlist[["bedInput"]],colClasses = c("NULL","integer","integer") )
        }



        readlens<-readslist[[2]]-readslist[[1]]

        #read length distribution
        #          allreadslen <-as.data.frame(readlens)
        #          colnames(allreadslen)<-"length"
        #          ggplot(allreadslen)+geom_histogram(bins = max(allreadslen),aes(x=length))
        #          ggsave(paste0(.Object@paramlist[["reportPrefix"]],".lendistr.pdf"))

        #readslen<-names(readscounts)
        #readslen<-cbind(readlen,readscounts)
        #ggplot(allreadslen)+geom_density(aes(x="length",fill="clarity"))

        #period distribution

        readscounts<-table(readlens)

        readscounts<-data.frame(readscounts)
        colnames(readscounts)<-c("length","counts")



        readscounts$length=as.integer(as.character(readscounts$length))
        mg<-data.frame(length=1:max(readscounts$length))
        readscounts <- merge(readscounts,mg,by="length",all = TRUE)
        readscounts$counts[is.na(readscounts$counts)]<-0

        
        
        write.table(x=readscounts,file = .Object@paramlist[["lendistrtxtOutput"]],quote = FALSE,row.names = FALSE,sep="\t")
        ggplot(readscounts[1:1000,], aes(length,counts))+geom_path(color="Red")+xlab("Fragment length (bp)")+ylab("Read counts") + theme_bw() + theme(panel.grid =element_blank()) 
        ggsave(.Object@paramlist[["lendistrpdfOutput"]])
        
        strength<-Mod(fft(readscounts$counts))/length(readscounts$counts)
        periodx<-length(readscounts$counts)/(1:(length(readscounts$counts)-1))
        strength<-strength[2:length(strength)]
        
        rs1<-as.data.frame(cbind(periodx[periodx<20&periodx>2],strength[periodx<20&periodx>2],0))
        rs2<-as.data.frame(cbind(periodx[periodx<400&periodx>2],strength[periodx<400&periodx>2],1))
        rs<-rbind(rs1,rs2)
        colnames(rs)<-c("period","strength","check")
        period<-"period"
        strength<-"strength"
        g1<-ggplot(rs[rs["check"]==0,]) + 
            geom_vline(xintercept = 10.4, linetype=2)+ 
            geom_line(aes(x=period,y=strength),color="Red")+ 
            theme_bw() + theme(panel.grid =element_blank()) + 
            annotate("text", x = 10.4, y = max(rs[rs["check"]==0,2]), 
                     label = "10.4bp") +xlab("period") + ylab("strength")
        ggsave(.Object@paramlist[["dnagroovepdfOutput"]])
        g2<-ggplot(rs[rs["check"]==1,]) + 
            geom_vline(xintercept = 186, linetype=2)+ 
            geom_line(aes(x=period,y=strength),color="Red")+ 
            theme_bw() + theme(panel.grid =element_blank()) + 
            annotate("text", x = 186, y = max(rs[rs["check"]==1,2]), 
                     label = "186bp") +xlab("period") + ylab("strength")  
        ggsave(.Object@paramlist[["histonepdfOutput"]])

        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "FregLenDistr",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bedInput"]])){
            stop("bedInput is required.")
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "FregLenDistr",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["bedInput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["lendistrpdfOutput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["lendistrtxtOutput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["dnagroovepdfOutput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["histonepdfOutput"]]);
    }
)


setMethod(
    f = "getReportValImp",
    signature = "FregLenDistr",
    definition = function(.Object, item){
        readscounts <- read.table(file= .Object@paramlist[["lendistrtxtOutput"]],header=TRUE)
        if(item == "readsCounts"){
            return(readscounts)
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "FregLenDistr",
    definition = function(.Object){
        return(c("readsCounts"))
    }
)

#' @name atacFregLenDistr
#' @title Quality control for fregment length distribution
#' @description
#' These functions are used to generate fregment distribution plot.
#' The fourier transform of fregment distribution will be calculated.
#' Strength distribution around period at 10.4bp and 180bp
#' will be shown in another two plots.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' \code{\link{atacBedUtils}}
#' \code{\link{bedUtils}}
#' @param reportPrefix \code{Character} scalar.
#' The prefix of report files path.
#' @param bedInput \code{Character} scalar.
#' BED file input path.
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{fregLenDistr} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' \code{\link{atacBedUtils}}
#' \code{\link{bedUtils}}
#'
#' @examples
#' 
#' library(R.utils)
#' td <- tempdir()
#' options(atacConf=setConfigure("tmpdir",td))
#'
#' bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#' fregLenDistr(bedfile)
#'
#' dir(td)
#' 
#' @importFrom BiocGenerics counts 
#' @importFrom ggplot2 geom_path ggplot geom_vline geom_line theme_bw theme annotate xlab ggsave element_blank
#' @name atacFregLenDistr
#' @export
#' @docType methods
#' @rdname atacFregLenDistr-methods

setGeneric("atacFregLenDistr",function(atacProc,reportPrefix=NULL,bedInput=NULL, ...) standardGeneric("atacFregLenDistr"))


#' @rdname atacFregLenDistr-methods
#' @aliases atacFregLenDistr
setMethod(
    f = "atacFregLenDistr",
    signature = "ATACProc",
    definition = function(atacProc,reportPrefix=NULL,bedInput=NULL, ...){
        atacproc <- new(
            "FregLenDistr",
            atacProc = atacProc,
            reportPrefix = reportPrefix,
            bedInput = bedInput)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)



#' @rdname atacFregLenDistr-methods
#' @export

fregLenDistr<-function(bedInput, reportPrefix=NULL, ...){
    atacproc <- new(
        "FregLenDistr",
        atacProc = NULL,
        reportPrefix = reportPrefix,
        bedInput = bedInput)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
