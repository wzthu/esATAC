setClass(Class = "PeakCallingMACS2",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "PeakCallingMACS2",
    definition = function(.Object,prevSteps = list(), ...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        background <- allparam[["background"]]
        outputPrefix <- allparam[["outputPrefix"]]
        genomeSize <- allparam[["genomeSize"]]
        pvalueThreshold <- allparam[["pvalueThreshold"]]
        extsize <- allparam[["extsize"]]
        shift <- allparam[["shift"]]
       
        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                atacProc <- prevSteps[[1]]
                atacProc<-c(unlist(atacProc),list())
                atacProc <- atacProc[[length(atacProc)]]
                input(.Object)[["bedInput"]] <- output(atacProc)[["bedOutput"]]
            }
        }
        
        if(!is.null(bedInput)){
            input(.Object)[["bedInput"]] <- bedInput
        }
        if(!is.null(background)){
            input(.Object)[["background"]] <- background
        }
        if(!is.null(outputPrefix)){
            output(.Object)[["bedOutput"]] <-getStepWorkDir(.Object,filename = paste0(outputPrefix,"_peaks.narrowPeak.bed"))
            output(.Object)[["narrowPeak"]] <-getStepWorkDir(.Object,filename = paste0(outputPrefix,"_peaks.narrowPeak"))
        }else{
            outputPrefix <- getStepWorkDir(.Object,'MACS')
            output(.Object)[["bedOutput"]] <-paste0(outputPrefix,"_peaks.narrowPeak.bed")
            output(.Object)[["narrowPeak"]] <-paste0(outputPrefix,"_peaks.narrowPeak")
        }
        param(.Object)[['outputPrefix']] <- outputPrefix
        if(is.null(genomeSize)){
            bsgenome<-getRefRc('bsgenome')
            lens<-seqlengths(bsgenome)
            chrs <- names(lens)
            chrs <- chrs[grep('_',chrs,invert = TRUE)]
            chrs <- chrs[chrs!='chrM']
            lens <- lens[chrs]
            genomeSize <- sum(lens)
            param(.Object)[['genomeSizes']] <- lens
        }else{
            param(.Object)[['genomeSizes']] <- NULL
        }
        param(.Object)[['genomeSize']] <- genomeSize
        param(.Object)[['pvalueThreshold']] <- pvalueThreshold
        param(.Object)[['extsize']] <- extsize
        param(.Object)[['shift']] <- shift
        
        
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "PeakCallingMACS2",
    definition = function(.Object,...){
        cmdline <- paste('macs2', 'callpeak', 
                          '-t', input(.Object)[['bedInput']],
                          '-f', 'BED',
                          '-n', param(.Object)[['outputPrefix']],
                          '-g', param(.Object)[['genomeSize']],
                          '-p', param(.Object)[['pvalueThreshold']],
                          '--shift', param(.Object)[['shift']],
                          '--extsize', param(.Object)[['extsize']],
                          '--nomodel -B --SPMR --keep-dup all --call-summits')
        result <- system(command = cmdline, intern = TRUE)
        for(aline in result){
            writeLog(.Object, msg = aline)
        }
        bedfile <- read.table(file = output(.Object)[['narrowPeak']],sep = '\t' ,header = FALSE)
        bedfile <- bedfile[,1:6]
        bedfile$V6 <-'*'
        write.table(bedfile, file = output(.Object)[['bedOutput']],sep = "\t", col.names = FALSE, row.names = FALSE , quote = FALSE)
        if(!is.null(param(.Object)[['genomeSizes']])){
            bedfile <- import.bed(output(.Object)[['bedOutput']])
            bsgenome <- getRefRc('bsgenome')
            sl <- seqlengths(bsgenome)
            chrs <- unique(seqnames(bedfile))
            sl<-sl[chrs]
            for(i in 1:length(sl)){
                a <- names(sl)[i]
                b <- sl[i]
                end(bedfile[(as.character(seqnames(bedfile))==a) & 
                            (end(bedfile)>as.integer(b))])  <- as.integer(b)
            }
            export.bed(bedfile, output(.Object)[['bedOutput']])
        }
       
        .Object
    }
)

setMethod(
    f = "genReport",
    signature = "PeakCallingMACS2",
    definition = function(.Object, ...){
        .Object
    }
)


#' @name PeakCallingMACS2
#' @title Use MACS2 to call peak
#' @description
#' Use  MACS2 installed by to call peak
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}},
#' \code{\link{atacBedUtils}}.
#' @param bedInput \code{Character} scalar.
#' BED file input path.
#' @param background \code{Character} scalar.
#' background directory default: NULL (none)
#' @param outputPrefix \code{Character} scalar.
#' the output bed file path
#' @param threshold \code{Numeric} scalar.
#' threshold (standard deviations) default: NULL (4.0)
#' @param extsize \code{Logical} scalar.
#' verbose output if TRUE.
#' @param shift \code{Character} scalar.
#' wg threshold set default: NULL (calculated)
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{peakCalling} instead.
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
#' library(magrittr)
#' td <- tempdir()
#' setTmpDir(td)
#'
#' bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#'
#' bedUtils(bedInput = bedfile,maxFragLen = 100, chrFilterList = NULL) %>%
#' atacPeakCalling
#'
#' dir(td)


setGeneric("atacPeakCallingMACS2",function(atacProc, bedInput=NULL, background=NULL,outputPrefix=NULL, genomeSize=NULL,
                                      pvalueThreshold=0.01, extsize=150, shift=-round(extsize/2.0), ...) standardGeneric("atacPeakCallingMACS2"))

#' @rdname PeakCallingMACS2
#' @aliases atacPeakCallingMACS2
#' @export
setMethod(
    f = "atacPeakCallingMACS2",
    signature = "ATACProc",
    definition = function(atacProc, bedInput=NULL, background=NULL,outputPrefix=NULL, genomeSize=NULL,
                          pvalueThreshold=0.01, extsize=150, shift=-round(extsize/2.0), ...){
        allpara <- c(list(Class = "PeakCallingMACS2", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname PeakCallingMACS2
#' @aliases peakCallingMACS2
#' @export
peakCallingMACS2 <- function(bedInput, background=NULL,outputPrefix=NULL, genomeSize=NULL,
                             pvalueThreshold=0.01, extsize=150, shift=-round(extsize/2.0), ...){
    allpara <- c(list(Class = "PeakCallingMACS2", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}

#' @rdname PeakCallingMACS2
#' @aliases testPeakCallingMACS2
#' @export
testPeakCallingMACS2 <- function(){
    tryCatch({system(command = 'macs2 --help', intern = T); return(TRUE)}, 
             error = function(e) {message('macs2 is not available for esATAC. Please install MACS2 first');return(FALSE)})
}
