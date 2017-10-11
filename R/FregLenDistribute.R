FregLenDistr <-R6Class(
  classname = "FregLenDistr",
  inherit = ATACProc,
  public = list(
    initialize = function(atacProc,reportPrefix=NULL,bedInput=NULL,editable=FALSE){
      super$initialize("FregLenDistr",editable,list(arg1=atacProc))
      if(private$singleEnd){
          private$writeLog("This process is for pair-end sequencing data.",isWarnning=TRUE)
      }
      if(!is.null(atacProc)){
        private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
        regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
      }else{
          regexProcName<-"(BED|bed|Bed)"
      }
      if(!is.null(bedInput)){
        private$paramlist[["bedInput"]] <- bedInput;
      }
      if(is.null(reportPrefix)){
          if(!is.null(private$paramlist[["bedInput"]])){
              prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
              private$paramlist[["lendistrpdfOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".lendistr.pdf"))
              private$paramlist[["lendistrtxtOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".lendistr.txt"))
              private$paramlist[["dnagroovepdfOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".dnagroove.pdf"))
              private$paramlist[["histonepdfOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".histone.pdf"))
          }
          #private$paramlist[["reportPrefix"]] <- paste(private$paramlist[["bedInput"]],".lenreport",sep="");
      }else{
        private$paramlist[["lendistrpdfOutput"]] <- paste0(reportPrefix,".lendistr.pdf")
        private$paramlist[["lendistrtxtOutput"]] <- paste0(reportPrefix,".lendistr.txt")
        private$paramlist[["dnagroovepdfOutput"]] <- paste0(reportPrefix,".dnagroove.pdf")
        private$paramlist[["histonepdfOutput"]] <- paste0(reportPrefix,".histone.pdf")
      }

        private$paramValidation()
    }

  ),
  private = list(
      processing = function(){
          readslist<-read.table(file = private$paramlist[["bedInput"]],nrows = 1)
          bedcol=length(colnames(readslist))
          if(bedcol>3){
              readslist<-read.table(file = private$paramlist[["bedInput"]],colClasses = c("NULL","integer","integer",rep("NULL",bedcol-3)))
          }else{
              readslist<-read.table(file = private$paramlist[["bedInput"]],colClasses = c("NULL","integer","integer") )
          }



          readlens<-readslist[[2]]-readslist[[1]]

          #read length distribution
#          allreadslen <-as.data.frame(readlens)
#          colnames(allreadslen)<-"length"
#          ggplot(allreadslen)+geom_histogram(bins = max(allreadslen),aes(x=length))
#          ggsave(paste0(private$paramlist[["reportPrefix"]],".lendistr.pdf"))

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

          ggplot(readscounts, aes(length))+geom_ribbon(aes(ymin=0, ymax=counts))
          ggsave(private$paramlist[["lendistrpdfOutput"]])
          write.table(x=readscounts,file = private$paramlist[["lendistrtxtOutput"]],quote = FALSE,row.names = FALSE,sep="\t")

          rs<-Mod(fft(readscounts$counts))/length(readscounts$counts)
          t<-length(readscounts$counts)/(1:(length(readscounts$counts)-1))
          rs<-rs[2:length(rs)]
          # t<-t[1:as.integer(length(t)/2)]
          #rs<-rs[1:as.integer(length(rs)/2)]
          tp<-rep(0,length(rs))
          if(length(t)>15&&sum(t>10&t<11)<=1){
              tp[max(which(t>10))+1]<-1
              tp[min(which(t<11))-1]<-1
          }else{
              tp[t>10&t<11]<-1
          }
          if(length(t)>220&&sum(t>100&t<200)<=1){
              tp[max(which(t>100))+1]<-1
              tp[min(which(t<200))-1]<-1
          }else{
              tp[t>100&t<200]<-2
          }
          if(length(t)>15){
              rs1<-as.data.frame(cbind(t[t<20&t>2],rs[t<20&t>2],tp[t<20&t>2]))
              #rs_1<-as.data.frame(cbind(t[t<20&t>2&tp==0],rs[t<20&t>2&tp==0],tp[t<20&t>2&tp==0]))
              #rs_2<-as.data.frame(cbind(t[t<20&t>2&tp!=0],rs[t<20&t>2&tp!=0],tp[t<20&t>2&tp!=0]))
              #colnames(rs_1)<-c("perior","strength","check")
              #colnames(rs_2)<-c("perior","strength","check")
              colnames(rs1)<-c("perior","strength","check")
              # ggplot(rs1)+geom_line(aes(x=perior,y=strength))+geom_vline(xintercept = 10)+geom_vline(xintercept = 11)
              #ggplot(rs1,aes(x=perior,y=strength))+geom_area(aes(fill="valence",color=check))
              checkdna=1
              ggplot(rs1)+geom_ribbon(data=subset(rs1,perior<=min(rs1$perior[rs1$check==checkdna])),aes(x=perior,ymin=0,ymax=strength),fill="blue")+geom_ribbon(data=subset(rs1,perior>=max(rs1$perior[rs1$check==checkdna])),aes(x=perior,ymin=0,ymax=strength),fill="blue")+geom_ribbon(data=subset(rs1,check==checkdna),aes(x=perior,ymin=0,ymax=strength),fill="red")
              ggsave(private$paramlist[["dnagroovepdfOutput"]])
          }
          if(length(t)>220){
              rs2<-as.data.frame(cbind(t[t<400&t>2],rs[t<400&t>2],tp[t<400&t>2]))
              #rs2<-as.data.frame(cbind(t[t<500&t>2&rs>10&rs<11],rs[t<500&t>2&rs>10&rs<11]))
              colnames(rs2)<-c("perior","strength","check")
              #ggplot(rs2,aes(x=perior,y=strength))+geom_area(aes(fill="valence"))
              #ggplot(rs2)+geom_line(aes(x=perior,y=strength))+geom_vline(xintercept = 150)+geom_vline(xintercept = 200)
              #ggplot(allreadslen)+geom_density(aes(x="strength",fill="clarity"))
              checkhistone=2
              ggplot(rs2)+geom_ribbon(data=subset(rs2,perior<=min(rs2$perior[rs2$check==checkhistone])),aes(x=perior,ymin=0,ymax=strength),fill="blue")+geom_ribbon(data=subset(rs2,perior>=max(rs2$perior[rs2$check==checkhistone])),aes(x=perior,ymin=0,ymax=strength),fill="blue")+geom_ribbon(data=subset(rs2,check==checkhistone),aes(x=perior,ymin=0,ymax=strength),fill="red")
              ggsave(private$paramlist[["histonepdfOutput"]])

          }




      },
    checkRequireParam = function(){
      if(is.null(private$paramlist[["bedInput"]])){
        stop("bedInput is required.")
      }


    },
    checkAllPath = function(){
        private$checkFileExist(private$paramlist[["bedInput"]])
        private$checkFileCreatable(private$paramlist[["lendistrpdfOutput"]])
        private$checkFileCreatable(private$paramlist[["lendistrtxtOutput"]])
        private$checkFileCreatable(private$paramlist[["dnagroovepdfOutput"]])
        private$checkFileCreatable(private$paramlist[["histonepdfOutput"]])

    },
    getReportValImp = function(item){
        readscounts <- read.table(file= private$paramlist[["lendistrtxtOutput"]],header=TRUE)
        if(item == "readsCounts"){
            return(readscounts)
        }
    },
    getReportItemsImp = function(){
        return(c("readsCounts"))
    }
  )


)

#' @name atacFregLenDistr
#' @aliases atacFregLenDistr
#' @aliases fregLenDistr
#' @title Quality control for fregment length distribution
#' @description
#' These functions are used to generate fregment distribution plot.
#' The fourier transform of fregment distribution will be calculated.
#' Strength distribution around period at 10.4bp and 180bp
#' will be shown in another two plots.
#' @param atacProc \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' \code{\link{atacBedUtils}}
#' \code{\link{bedUtils}}
#' @param reportPrefix \code{Character} scalar.
#' The prefix of report files path.
#' @param bedInput \code{Character} scalar.
#' BED file input path.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{fregLenDistr} instead.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
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
#' setConfigure("tmpdir",td)
#'
#' bedbzfile <- system.file(package="ATACpipe", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#' fregLenDistr(bedfile)
#'
#' dir(td)
#'

#' @rdname atacFregLenDistr
#' @export
atacFregLenDistr<-function(atacProc,reportPrefix=NULL,bedInput=NULL){
    atacproc<-FregLenDistr$new(atacProc = atacProc,reportPrefix = reportPrefix,bedInput = bedInput)
    atacproc$process()
    invisible(atacproc)
}
#' @rdname atacFregLenDistr
#' @export
fregLenDistr<-function(bedInput, reportPrefix=NULL){
    atacproc<-FregLenDistr$new(atacProc = NULL,reportPrefix = reportPrefix,bedInput = bedInput)
    atacproc$process()
    invisible(atacproc)
}
