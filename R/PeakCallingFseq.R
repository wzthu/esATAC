PeakCallingFseq <-R6Class(
    classname = "PeakCallingFseq",
    inherit = BaseProc,
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
                private$paramlist[["inputDir"]] <- dirname(atacProc$getParam("bedOutput"));
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }

            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
                private$paramlist[["bedFileList"]] <- basename(bedInput)
                private$paramlist[["inputDir"]] <- dirname(bedInput);
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
            private$paramlist[["outputDir"]] <- paste0(private$paramlist[["bedOutput"]],".tmp");

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
            dir.create(private$paramlist[["outputDir"]])
            .fseq_call(bedFileList=private$paramlist[["bedFileList"]],
                                 background=private$paramlist[["background"]],
                                 genomicReadsCount=private$paramlist[["genomicReadsCount"]],
                                 inputDir=private$paramlist[["inputDir"]],
                                 fragmentSize=private$paramlist[["fragmentSize"]],
                                 featureLength=private$paramlist[["featureLength"]],
                                 outputDir=private$paramlist[["outputDir"]],
                                 outputFormat=private$paramlist[["fileformat"]],
                                 ploidyDir=private$paramlist[["ploidyDir"]],
                                 wiggleTrackStep=private$paramlist[["wiggleTrackStep"]],
                                 threshold=private$paramlist[["threshold"]],
                                 verbose=private$paramlist[["verbose"]],
                                 wgThresholdSet=private$paramlist[["wgThresholdSet"]]);

            filename<-list.files(private$paramlist[["outputDir"]])
            for(i in 1:length(filename)){
                filename[i]<-strsplit(filename[i],split="\\.")[[1]][1]
            }
            peakfiles <- sort(filename)
            peakfiles<-paste0(peakfiles,".bed")
            peakfiles <- paste0(private$paramlist[["outputDir"]],"/",peakfiles)
            file.create(private$paramlist[["bedOutput"]])
            for(i in 1:length(peakfiles)){
                file.append(private$paramlist[["bedOutput"]],peakfiles[i])
            }
            #mergeFile(private$paramlist[["bedOutput"]],peakfiles)
            unlink(private$paramlist[["outputDir"]],recursive = TRUE,force = TRUE)
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

atacPeakCalling <- function(atacProc,bedInput=NULL,background=NULL,genomicReadsCount=NULL,
                            fragmentSize=0,featureLength=NULL,bedOutput=NULL,
                             ploidyDir=NULL,#fileformat=c("bed","wig","npf"),
                            wiggleTrackStep=NULL,threshold=NULL,verbose=TRUE,
                            wgThresholdSet=NULL){
    peakcalling <- PeakCallingFseq$new(atacProc,bedInput,background,genomicReadsCount,
                                       fragmentSize,featureLength,bedOutput,fileformat="bed", ploidyDir,
                                       wiggleTrackStep,threshold,verbose,wgThresholdSet)
    peakcalling$process();
    return(peakcalling)
}
