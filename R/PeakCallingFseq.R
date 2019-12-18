setClass(Class = "PeakCallingFseq",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "PeakCallingFseq",
    definition = function(.Object,prevSteps = list(), ...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        background <- allparam[["background"]]
        genomicReadsCount <- allparam[["genomicReadsCount"]]
        fragmentSize <- allparam[["fragmentSize"]]
        featureLength <- allparam[["featureLength"]]
        bedOutput <- allparam[["bedOutput"]]
        ploidyDir <- allparam[["ploidyDir"]]
        wiggleTrackStep <- allparam[["wiggleTrackStep"]]
        threshold <- allparam[["threshold"]]
        verbose <- allparam[["verbose"]]
        wgThresholdSet <- allparam[["wgThresholdSet"]]
        fileformat <- allparam[["fileformat"]]
        
        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                atacProc <- prevSteps[[1]]
                atacProc<-c(unlist(atacProc),list())
                atacProc <- atacProc[[length(atacProc)]]
                input(.Object)[["bedInput"]] <- output(atacProc)[["bedOutput"]]
                param(.Object)[["bedFileList"]] <- basename(output(atacProc)[["bedOutput"]])
                param(.Object)[["inBedDir"]] <- dirname(output(atacProc)[["bedOutput"]])
            }
        }

        if(!is.null(bedInput)){
            input(.Object)[["bedInput"]] <- bedInput;
            input(.Object)[["bedFileList"]] <- basename(bedInput)
            input(.Object)[["inBedDir"]] <- dirname(bedInput)
        }
        if(!is.null(bedOutput)){
            output(.Object)[["bedOutput"]] <-bedOutput
        }else{
            if(!is.null(input(.Object)[["bedInput"]])){
                output(.Object)[["bedOutput"]] <- getAutoPath(.Object, input(.Object)[["bedInput"]], "BED|Bed|bed", "bed")
            }
            #.Object@paramlist[["bedOutput"]] <-paste(.Object@paramlist[["bedInput"]],".peak.bed",sep="");
        }
        param(.Object)[["outTmpDir"]] <- paste0(output(.Object)[["bedOutput"]],".tmp");

        param(.Object)[["background"]] <- background;
        param(.Object)[["genomicReadsCount"]] <- genomicReadsCount;
        param(.Object)[["fragmentSize"]] <- fragmentSize;
        param(.Object)[["featureLength"]] <- featureLength;
        param(.Object)[["fileformat"]] <- fileformat[1];
        param(.Object)[["ploidyDir"]] <- ploidyDir;
        param(.Object)[["wiggleTrackStep"]] <- wiggleTrackStep;
        param(.Object)[["threshold"]] <- threshold;
        param(.Object)[["verbose"]] <- verbose;
        param(.Object)[["wgThresholdSet"]] <- wgThresholdSet;


        .Object
    }
)

setMethod(
    f = "processing",
    signature = "PeakCallingFseq",
    definition = function(.Object,...){
        dir.create(param(.Object)[["outTmpDir"]])
        .fseq_call(bedFileList=param(.Object)[["bedFileList"]],
                   background=param(.Object)[["background"]],
                   genomicReadsCount=param(.Object)[["genomicReadsCount"]],
                   inputDir=param(.Object)[["inBedDir"]],
                   fragmentSize=param(.Object)[["fragmentSize"]],
                   featureLength=param(.Object)[["featureLength"]],
                   outputDir=param(.Object)[["outTmpDir"]],
                   outputFormat=param(.Object)[["fileformat"]],
                   ploidyDir=param(.Object)[["ploidyDir"]],
                   wiggleTrackStep=param(.Object)[["wiggleTrackStep"]],
                   threshold=param(.Object)[["threshold"]],
                   verbose=param(.Object)[["verbose"]],
                   wgThresholdSet=param(.Object)[["wgThresholdSet"]]);

        filename<-list.files(param(.Object)[["outTmpDir"]])
        for(i in 1:length(filename)){
            filename[i]<-strsplit(filename[i],split="\\.")[[1]][1]
        }
        peakfiles <- sort(filename)
        peakfiles <-paste0(peakfiles,".bed")
        peakfiles <- file.path(param(.Object)[["outTmpDir"]],peakfiles)
        file.create(output(.Object)[["bedOutput"]])
        for(i in 1:length(peakfiles)){
            if(file.exists(peakfiles[i])&&file.info(peakfiles[i])$size>0){
                df <- read.table(peakfiles[i], header = FALSE,sep = "\t")
                df <- cbind(df,"*")
                write.table(x=df,file = output(.Object)[["bedOutput"]],append = TRUE,
                            quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
            }
            #file.append(output(.Object)[["bedOutput"]],peakfiles[i])
        }
        #mergeFile(output(.Object)[["bedOutput"]],peakfiles)
#        unlink(.Object@paramlist[["outTmpDir"]],recursive = TRUE,force = TRUE)


        .Object
    }
)

setMethod(
    f = "genReport",
    signature = "PeakCallingFseq",
    definition = function(.Object, ...){
        .Object
    }
)


#' @name PeakCallingFseq
#' @title Use F-seq to call peak
#' @description
#' Use F-seq to call peak
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}},
#' \code{\link{atacBedUtils}}.
#' @param bedInput \code{Character} scalar.
#' BED file input path.
#' @param background \code{Character} scalar.
#' background directory default: NULL (none)
#' @param genomicReadsCount \code{Integer} scalar.
#' genomic count of sequence reads. default: NULL (calculated)
#' @param fragmentSize \code{Integer} scalar.
#' fragment size. set NULL to estimat from data. default:0
#' @param featureLength \code{Character} scalar.
#' feature length default: NULL (600)
#' @param bedOutput \code{Character} scalar.
#' the output bed file path
#' @param ploidyDir \code{Character} scalar.
#' ploidy/input directory. default: NULL
#' @param fileformat \code{Character} scalar.
#' File format of result. default: bed
#' @param wiggleTrackStep \code{Integer} scalar.
#' wiggle track step default: NULL (1)
#' @param threshold \code{Numeric} scalar.
#' threshold (standard deviations) default: NULL (4.0)
#' @param verbose \code{Logical} scalar.
#' verbose output if TRUE.
#' @param wgThresholdSet \code{Character} scalar.
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


setGeneric("atacPeakCalling",function(atacProc,bedInput=NULL,background=NULL,genomicReadsCount=NULL,
                                         fragmentSize=0,featureLength=NULL,bedOutput=NULL,
                                         ploidyDir=NULL,fileformat=c("bed","wig","npf"),
                                         wiggleTrackStep=NULL,threshold=NULL,verbose=TRUE,

                                         wgThresholdSet=NULL, ...) standardGeneric("atacPeakCalling"))

#' @rdname PeakCallingFseq
#' @aliases atacPeakCalling
#' @export
setMethod(
    f = "atacPeakCalling",
    signature = "ATACProc",
    definition = function(atacProc,bedInput=NULL,background=NULL,genomicReadsCount=NULL,
                          fragmentSize=0,featureLength=NULL,bedOutput=NULL,
                          ploidyDir=NULL,fileformat=c("bed","wig","npf"),
                          wiggleTrackStep=NULL,threshold=NULL,verbose=TRUE,
                          wgThresholdSet=NULL, ...){
        allpara <- c(list(Class = "PeakCallingFseq", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname PeakCallingFseq
#' @aliases peakCalling
#' @export
peakCalling <- function(bedInput, background=NULL,genomicReadsCount=NULL,
                            fragmentSize=0,featureLength=NULL,bedOutput=NULL,
                            ploidyDir=NULL,fileformat=c("bed","wig","npf"),
                            wiggleTrackStep=NULL,threshold=NULL,verbose=TRUE,
                            wgThresholdSet=NULL, ...){
    allpara <- c(list(Class = "PeakCallingFseq", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
