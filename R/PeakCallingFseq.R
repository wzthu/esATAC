PeakCallingFseq <-R6Class(
    classname = "PeakCallingFseq",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc,bedInput=NULL,background=NULL,genomicReadsCount=NULL,
                              fragmentSize=0,featureLength=NULL,bedOutput=NULL,
                              fileformat=c("bed","wig","npf"), ploidyDir=NULL,
                              wiggleTrackStep=NULL,threshold=NULL,verbose=TRUE,
                              wgThresholdSet=NULL,editable=FALSE){
            super$initialize("PeakCallingFseq",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                private$paramlist[["bedFileList"]] <- basename(atacProc$getParam("bedOutput"));
                private$paramlist[["inBedDir"]] <- dirname(atacProc$getParam("bedOutput"));
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }

            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
                private$paramlist[["bedFileList"]] <- basename(bedInput)
                private$paramlist[["inBedDir"]] <- dirname(bedInput);
            }
            if(!is.null(bedOutput)){
                private$paramlist[["bedOutput"]] <-bedOutput
            }else{
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["bedOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".bed"))
                }
                #private$paramlist[["bedOutput"]] <-paste(private$paramlist[["bedInput"]],".peak.bed",sep="");
            }
            private$paramlist[["outTmpDir"]] <- paste0(private$paramlist[["bedOutput"]],".tmp");

            private$paramlist[["background"]] <- background;
            private$paramlist[["genomicReadsCount"]] <- genomicReadsCount;
            private$paramlist[["fragmentSize"]] <- fragmentSize;
            private$paramlist[["featureLength"]] <- featureLength;
            private$paramlist[["fileformat"]] <- fileformat[1];
            private$paramlist[["ploidyDir"]] <- ploidyDir;
            private$paramlist[["wiggleTrackStep"]] <- wiggleTrackStep;
            private$paramlist[["threshold"]] <- threshold;
            private$paramlist[["verbose"]] <- verbose;
            private$paramlist[["wgThresholdSet"]] <- wgThresholdSet;



            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){
            dir.create(private$paramlist[["outTmpDir"]])
            .fseq_call(bedFileList=private$paramlist[["bedFileList"]],
                                 background=private$paramlist[["background"]],
                                 genomicReadsCount=private$paramlist[["genomicReadsCount"]],
                                 inputDir=private$paramlist[["inBedDir"]],
                                 fragmentSize=private$paramlist[["fragmentSize"]],
                                 featureLength=private$paramlist[["featureLength"]],
                                 outputDir=private$paramlist[["outTmpDir"]],
                                 outputFormat=private$paramlist[["fileformat"]],
                                 ploidyDir=private$paramlist[["ploidyDir"]],
                                 wiggleTrackStep=private$paramlist[["wiggleTrackStep"]],
                                 threshold=private$paramlist[["threshold"]],
                                 verbose=private$paramlist[["verbose"]],
                                 wgThresholdSet=private$paramlist[["wgThresholdSet"]]);

            filename<-list.files(private$paramlist[["outTmpDir"]])
            for(i in 1:length(filename)){
                filename[i]<-strsplit(filename[i],split="\\.")[[1]][1]
            }
            peakfiles <- sort(filename)
            peakfiles<-paste0(peakfiles,".bed")
            peakfiles <- paste0(private$paramlist[["outTmpDir"]],"/",peakfiles)
            file.create(private$paramlist[["bedOutput"]])
            for(i in 1:length(peakfiles)){
                file.append(private$paramlist[["bedOutput"]],peakfiles[i])
            }
            #mergeFile(private$paramlist[["bedOutput"]],peakfiles)
            unlink(private$paramlist[["outTmpDir"]],recursive = TRUE,force = TRUE)
            private$setFinish()
        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileExist(private$paramlist[["background"]]);
            private$checkPathExist(private$paramlist[["ploidyDir"]]);
            private$checkFileCreatable(private$paramlist[["bedOutput"]]);
        }
    )


)
#' @name atacPeakCalling
#' @aliases atacPeakCalling
#' @aliases peakCalling
#' @title Use F-seq to call peak
#' @description
#' Use F-seq to call peak
#' @param atacProc \code{\link{ATACProc}} object scalar.
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
#' @param wiggleTrackStep \code{Integer} scalar.
#' wiggle track step default: NULL (1)
#' @param threshold \code{Numeric} scalar.
#' threshold (standard deviations) default: NULL (4.0)
#' @param verbose \code{Logical} scalar.
#' verbose output if TRUE.
#' @param wgThresholdSet \code{Character} scalar.
#' wg threshold set default: NULL (calculated)
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{peakCalling} instead.
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
#' library(magrittr)
#' td <- tempdir()
#' setConfigure("tmpdir",td)
#' 
#' bedbzfile <- system.file(package="ATACFlow", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr18.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#' 
#' bedUtils(bedInput = bedfile,maxFregLen = 100, chrFilterList = NULL) %>%
#' atacPeakCalling
#' 
#' dir(td) 

#' @rdname atacPeakCalling
#' @export
atacPeakCalling <- function(atacProc,bedInput=NULL,background=NULL,genomicReadsCount=NULL,
                            fragmentSize=0,featureLength=NULL,bedOutput=NULL,
                             ploidyDir=NULL,#fileformat=c("bed","wig","npf"),
                            wiggleTrackStep=NULL,threshold=NULL,verbose=TRUE,
                            wgThresholdSet=NULL){
    peakcalling <- PeakCallingFseq$new(atacProc,bedInput,background,genomicReadsCount,
                                       fragmentSize,featureLength,bedOutput,fileformat="bed", ploidyDir,
                                       wiggleTrackStep,threshold,verbose,wgThresholdSet)
    peakcalling$process();
    invisible(peakcalling)
}
#' @rdname atacPeakCalling
#' @export
peakCalling <- function(bedInput,background=NULL,genomicReadsCount=NULL,
                            fragmentSize=0,featureLength=NULL,bedOutput=NULL,
                            ploidyDir=NULL,#fileformat=c("bed","wig","npf"),
                            wiggleTrackStep=NULL,threshold=NULL,verbose=TRUE,
                            wgThresholdSet=NULL){
    peakcalling <- PeakCallingFseq$new(atacProc = NULL,bedInput,background,genomicReadsCount,
                                       fragmentSize,featureLength,bedOutput,fileformat="bed", ploidyDir,
                                       wiggleTrackStep,threshold,verbose,wgThresholdSet)
    peakcalling$process();
    invisible(peakcalling)
}
