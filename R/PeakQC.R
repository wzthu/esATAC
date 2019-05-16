setClass(Class = "PeakQC",
         contains = "ATACProc",
         slots = list(fixtag = "character"),
         prototype = list(fixtag = "User file")
)



setMethod(
    f = "init",
    signature = "PeakQC",
    definition = function(.Object,prevSteps = list(), ...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        bsgenome <- allparam[["bsgenome"]]
        reportOutput <- allparam[["reportOutput"]]
        qcbedInput <- allparam[["qcbedInput"]]
       
        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                atacProc <- prevSteps[[1]]
                atacProc<-c(unlist(atacProc),list())
                atacProc <- atacProc[[length(atacProc)]]
                input(.Object)[["bedInput"]] <- output(atacProc)[["bedOutput"]]
            }
        }
        print("-------------------")        
        qcbedInput <- qcbedInput[1]
        print("-------------------")
        if(qcbedInput == "DHS"){
            input(.Object)[["qcbedInput"]]<-getRefFiles("DHS");
            .Object@fixtag = "DHS"
        }else if(qcbedInput == "blacklist"){
            input(.Object)[["qcbedInput"]]<-getRefFiles("blacklist");
            .Object@fixtag = "blacklist"
        }else{
            print("-------------------")
            input(.Object)[["qcbedInput"]]<-qcbedInput;
        }
        print("-------------------")
        if(!is.null(bedInput)){
            input(.Object)[["bedInput"]] <- bedInput;
        }
        print("-------------------")
        if(is.null(reportOutput)){
            if(!is.null(input(.Object)[["bedInput"]])){
                print("-------------------")
                output(.Object)[["reportOutput"]] <- getAutoPath(.Object, input(.Object)[["bedInput"]], "BED|Bed|bed","report.txt")
            }
        }else{
            print("-------------------")
            output(.Object)[["reportOutput"]] <- reportOutput;
        }
        print("-------------------")
        param(.Object)[["bsgenome"]] <- bsgenome
        print("-------------------")
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "PeakQC",
    definition = function(.Object,...){
        if(is.null(param(.Object)[["bsgenome"]])){
            genome <- seqinfo(getRefRc("bsgenome"))
        }else{
            genome <- seqinfo(param(.Object)[["bsgenome"]])
        }
print("-------------------")
#        inputbed <- import(con = .Object@paramlist[["bedInput"]], genome = genome,format = "bed")
        inputbed <- import(con = input(.Object)[["bedInput"]], format = "bed")

#        qcbedInput<-import(con = .Object@paramlist[["qcbedInput"]], genome = genome,format = "bed")
        qcbedInput<-import(con = input(.Object)[["qcbedInput"]], format = "bed")


        qcval=list();
        print("-------------------")
        qcval[["totalInput"]]<-length(inputbed)
        print("-------------------")
        qcval[["qcbedInput"]]<-length(subsetByOverlaps(inputbed, qcbedInput,ignore.strand = TRUE))
        print("-------------------")
        qcval[["qcbedRate"]]<-qcval[["qcbedInput"]]/qcval[["totalInput"]]
        print("-------------------")
        print(as.data.frame(qcval))
        print( output(.Object)[["reportOutput"]])
        
        write.table(as.data.frame(qcval),file = output(.Object)[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        print("-------------------")
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "PeakQC",
    definition = function(.Object,...){
        if(is.null(input(.Object)[["bedInput"]])){
            stop("bedInput is required.")
        }
        if(is.null(input(.Object)[["qcbedInput"]])){
            stop("qcbedInput is required.")
        }
    }
)


setMethod(
    f = "genReport",
    signature = "PeakQC",
    definition = function(.Object, ...){
        qcval <- as.list(read.table(file= output(.Object)[["reportOutput"]],header=TRUE))
        cqcval<-as.character(qcval)
        cqcval[3]<-sprintf("%.2f",as.numeric(cqcval[[3]]))
        if(.Object@fixtag=="DHS"){
            report(.Object)$report <- (data.frame(Item=c("Total peaks","DHS regions", "Ratio"),
                                                  Value=cqcval))
        }else if(.Object@fixtag=="blacklist"){
            report(.Object)$report <- (data.frame(Item=c("Total peaks","Blacklist regions", "Ratio"),Value=cqcval))
        }else{
            report(.Object)$report <- (data.frame(Item=names(qcval),Value=as.character(qcval)))
        }
        n <- names(qcval)
        report(.Object) <- c(report(.Object),qcval)
        .Object
    }
)





#' @name PeakQC
#' @title Quality control for peak overlap
#' @description
#' These functions are used to 
#' calculate the overlap ratio in specific quality control rigion.
#' Blacklist and DHS region are provided. 
#' You can also set your own BED file as quality control rigion.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}},
#' \code{\link{atacBedUtils}}.
#' @param reportOutput \code{Character} scalar.
#' The report file path.
#' @param bsgenome \code{BSGenome} object scalar.
#' BSGenome object for specific species.
#' @param qcbedInput \code{Character} scalar.
#' It can be "DHS","blacklist" or
#' Other quality control BED file input path.
#' @param bedInput \code{Character} scalar.
#' BED file input path for quality control.
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{peakQC} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{atacBedUtils}}
#'
#' @examples
#' library(R.utils)
#' library(magrittr)
#' td <- tempdir()
#' setTmpDir(td)
#'
#' bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#' blacklistfile <- system.file(package="esATAC", "extdata", "hg19.blacklist.bed")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' bedUtils(bedInput = bedfile,maxFragLen = 100, chrFilterList = NULL) %>%
#' atacPeakCalling %>% atacPeakQC(qcbedInput = blacklistfile, bsgenome = BSgenome.Hsapiens.UCSC.hg19)
#' dir(td)


setGeneric("atacPeakQC",function(atacProc, bsgenome = NULL,
                                 reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"),
                                 bedInput = NULL, ...) standardGeneric("atacPeakQC"))
#' @rdname PeakQC
#' @aliases atacPeakQC
#' @export
setMethod(
    f = "atacPeakQC",
    signature = "ATACProc",
    definition = function(atacProc, bsgenome = NULL,
                          reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"),
                          bedInput = NULL, ...){
tryCatch(
    {
        allpara <- c(list(Class = "PeakQC", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    },
    error = function(cond){
        if(qcbedInput == "DHS" || qcbedInput == 'blacklist'){
            message('genome is not configured or')
            print(paste(qcbedInput,'is not available for current configured genome'))
            return(NULL)
        }else{
            stop(paste('qcbedInput:', qcbedInput,'does not exist'))
        }
    }
)    
}
)
   

#' @rdname PeakQC
#' @aliases peakQC
#' @export
peakQC<-function(bedInput, bsgenome = NULL, reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"), ...){
tryCatch(
    {
        allpara <- c(list(Class = "PeakQC", prevSteps = list()),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
},
    error = function(cond){
        if(qcbedInput == "DHS" || qcbedInput == 'blacklist'){
            message('genome is not configured or')
            print(paste(qcbedInput,'is not available for current configured genome'))
            return(NULL)
        }else{
            stop(paste('qcbedInput:', qcbedInput,'does not exist'))
        }
    }
)
}
