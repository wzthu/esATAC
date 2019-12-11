setClass(Class = "FragLenDistr",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "FragLenDistr",
    definition = function(.Object,prevSteps = list(), ...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        reportPrefix <- allparam[["reportPrefix"]]
        

#        if(property(.Object)[["singleEnd"]]){
#            writeLog(.Object,"This process is for pair-end sequencing data.",isWarnning=TRUE)
#        }
        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                atacProc <- prevSteps[[1]]
                atacProc<-c(unlist(atacProc),list())
                atacProc <- atacProc[[length(atacProc)]]
                input(.Object)[["bedInput"]] <- output(atacProc)[["bedOutput"]]
            }
        }
        print('--------------------')
        if(!is.null(bedInput)){
            input(.Object)[["bedInput"]] <- bedInput;
        }
        print('--------------------')
        if(is.null(reportPrefix)){
            if(!is.null(input(.Object)[["bedInput"]])){
                output(.Object)[["lendistrpdfOutput"]] <- getAutoPath(.Object, input(.Object)[["bedInput"]], "BED|Bed|bed","lendistr.pdf")
                output(.Object)[["lendistrtxtOutput"]] <- getAutoPath(.Object, input(.Object)[["bedInput"]], "BED|Bed|bed","lendistr.txt")
                output(.Object)[["dnagroovepdfOutput"]] <- getAutoPath(.Object, input(.Object)[["bedInput"]], "BED|Bed|bed","dnagroove.pdf")
                output(.Object)[["histonepdfOutput"]] <- getAutoPath(.Object, input(.Object)[["bedInput"]], "BED|Bed|bed","histone.pdf")
            }
        }else{
            output(.Object)[["lendistrpdfOutput"]] <- paste0(reportPrefix,".lendistr.pdf")
            output(.Object)[["lendistrtxtOutput"]] <- paste0(reportPrefix,".lendistr.txt")
            output(.Object)[["dnagroovepdfOutput"]] <- paste0(reportPrefix,".dnagroove.pdf")
            output(.Object)[["histonepdfOutput"]] <- paste0(reportPrefix,".histone.pdf")
        }
        print('--------------------')
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "FragLenDistr",
    definition = function(.Object,...){
        readslist<-read.table(file = input(.Object)[["bedInput"]],nrows = 1)
        bedcol=length(colnames(readslist))
        if(bedcol>3){
            readslist<-read.table(file = input(.Object)[["bedInput"]],colClasses = c("NULL","integer","integer",rep("NULL",bedcol-3)))
        }else{
            readslist<-read.table(file = input(.Object)[["bedInput"]],colClasses = c("NULL","integer","integer") )
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

        
        
        write.table(x=readscounts,file = output(.Object)[["lendistrtxtOutput"]],quote = FALSE,row.names = FALSE,sep="\t")
        ggplot(readscounts[1:1000,], aes(length,counts))+geom_path(color="Red")+xlab("Fragment length (bp)")+ylab("Read counts") + theme_bw() + theme(panel.grid =element_blank()) 
        ggsave(output(.Object)[["lendistrpdfOutput"]])
        
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
        ggsave(output(.Object)[["dnagroovepdfOutput"]])
        g2<-ggplot(rs[rs["check"]==1,]) + 
            geom_vline(xintercept = 186, linetype=2)+ 
            geom_line(aes(x=period,y=strength),color="Red")+ 
            theme_bw() + theme(panel.grid =element_blank()) + 
            annotate("text", x = 186, y = max(rs[rs["check"]==1,2]), 
                     label = "186bp") +xlab("period") + ylab("strength")  
        ggsave(output(.Object)[["histonepdfOutput"]])

        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "FragLenDistr",
    definition = function(.Object,...){
        if(is.null(input(.Object)[["bedInput"]])){
            stop("bedInput is required.")
        }
    }
)



setMethod(
    f = "genReport",
    signature = "FragLenDistr",
    definition = function(.Object, ...){
        readscounts <- read.table(file= output(.Object)[["lendistrtxtOutput"]],header=TRUE)
        report(.Object)$readsCounts <- readscounts
        .Object
    }
)



#' @name FragLenDistr
#' @title Quality control for fragment length distribution
#' @description
#' These functions are used to generate fragment distribution plot.
#' The fourier transform of fragment distribution will be calculated.
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
#' you can use \code{fragLenDistr} instead.
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
#' setTmpDir(td)
#'
#' bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' \dontrun{
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#' fragLenDistr(bedfile)
#' }
#'
#' dir(td)
#' 
#' @importFrom BiocGenerics counts 
#' @importFrom ggplot2 geom_path ggplot geom_vline geom_line theme_bw theme annotate xlab ggsave element_blank


setGeneric("atacFragLenDistr",function(atacProc,reportPrefix=NULL,bedInput=NULL, ...) standardGeneric("atacFragLenDistr"))


#' @rdname FragLenDistr
#' @aliases atacFragLenDistr
#' @export
setMethod(
    f = "atacFragLenDistr",
    signature = "ATACProc",
    definition = function(atacProc,reportPrefix=NULL,bedInput=NULL, ...){
        allpara <- c(list(Class = "FragLenDistr", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)



#' @rdname FragLenDistr
#' @aliases fragLenDistr
#' @export

fragLenDistr<-function(bedInput, reportPrefix=NULL, ...){
    allpara <- c(list(Class = "FragLenDistr", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
