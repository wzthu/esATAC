LibComplexQC <-R6Class(
    classname = "LibComplexQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc,reportPrefix=NULL,samInput=NULL,paired = FALSE,subsample=TRUE,subsampleSize=4*10e6,editable=FALSE){
            super$initialize("LibComplexQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["samInput"]] <- atacProc$getParam("samOutput");
            }
            if(!is.null(samInput)){
                private$paramlist[["samInput"]] <- samInput;
            }
            if(is.null(reportPrefix)){
                private$paramlist[["reportPrefix"]] <- paste(private$paramlist[["samInput"]],".libcomplexreport",sep="");
            }else{
                private$paramlist[["reportPrefix"]] <- reportPrefix;
            }

            private$paramlist[["subsample"]] <- subsample;
            private$paramlist[["subsampleSize"]] <- subsampleSize;
            private$paramlist[["paired"]]<-paired

            private$checkFileExist(private$paramlist[["samInput"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
            private$checkRequireParam();
        },
        processing = function(){
            super$processing()
            if(private$paramlist[["paired"]]){
                .sam2bed_merge_call(samfile = private$paramlist[["samInput"]], bedfile = paste0(private$paramlist[["reportPrefix"]],".tmp"),
                                    posOffset = 0, negOffset = 0,sortBed = !private$paramlist[["subsample"]],
                                    uniqueBed = FALSE, filterList = NULL,minFregLen = 0,maxFregLen = 1000000,saveExtLen = FALSE )
            }else{
                .sam2bed_call(samfile = private$paramlist[["samInput"]], bedfile = paste0(private$paramlist[["reportPrefix"]],".tmp"),
                              posOffset = 0, negOffset = 0, sortBed = !private$paramlist[["subsample"]],uniqueBed = FALSE,  filterList = NULL)
            }

            qcval<-.lib_complex_qc_call(bedfile=paste0(private$paramlist[["reportPrefix"]],".tmp"), sortedBed=!private$paramlist[["subsample"]], max_reads=private$paramlist[["subsampleSize"]])

            unlink(paste0(private$paramlist[["reportPrefix"]],".tmp"))
            private$paramlist[["qcval"]]<-qcval
            qcval<-as.matrix(qcval)
            write.table(qcval,file = paste0(private$paramlist[["reportPrefix"]],".txt"),sep="\t",quote = FALSE,col.names = FALSE)
            private$finish <- TRUE
        },
        setResultParam = function(fastqOutput1, fastqOutput2=NULL){
            super$setResultParam();
            private$paramlist[["fastqOutput1"]] <- fastqOutput1
            private$paramlist[["fastqOutput2"]] <- fastqOutput2
        }
    ),
    private = list(
        checkRequireParam = function(){
            if(private$editable){
                return();
            }
            if(is.null(private$paramlist[["samInput"]])){
                stop("samInput is required.")
            }


        }
    )


)
