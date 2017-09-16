LibComplexQC <-R6Class(
    classname = "LibComplexQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc,reportOutput=NULL,samInput=NULL,singleEnd = FALSE,subsample=TRUE,subsampleSize=4e6,editable=FALSE){
            super$initialize("LibComplexQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["samInput"]] <- atacProc$getParam("samOutput");
                regexProcName<-sprintf("(SAM|Sam|sam|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(SAM|Sam|sam)"
                private$singleEnd<-singleEnd
            }
            if(!is.null(samInput)){
                private$paramlist[["samInput"]] <- samInput;
            }
            if(is.null(reportOutput)){
                if(!is.null(private$paramlist[["samInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["samInput"]],regexProcName)
                    private$paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
            }else{
                private$paramlist[["reportOutput"]] <- reportOutput;
            }

            private$paramlist[["subsample"]] <- subsample;
            private$paramlist[["subsampleSize"]] <- subsampleSize;



            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){
            if(!private$singleEnd){
                .sam2bed_merge_call(samfile = private$paramlist[["samInput"]], bedfile = paste0(private$paramlist[["reportOutput"]],".tmp"),
                                    posOffset = 0, negOffset = 0,sortBed = !private$paramlist[["subsample"]],
                                    uniqueBed = FALSE, filterList = NULL,minFregLen = 0,maxFregLen = 1000000,saveExtLen = FALSE ,downSample=private$paramlist[["subsampleSize"]])
                print("test00")
            }else{
                .sam2bed_call(samfile = private$paramlist[["samInput"]], bedfile = paste0(private$paramlist[["reportOutput"]],".tmp"),
                              posOffset = 0, negOffset = 0, sortBed = !private$paramlist[["subsample"]],uniqueBed = FALSE,  filterList = NULL,downSample=private$paramlist[["subsampleSize"]])
                print("test01")
            }
            print("test1")
            qcval<-.lib_complex_qc_call(bedfile=paste0(private$paramlist[["reportOutput"]],".tmp"), sortedBed=!private$paramlist[["subsample"]], max_reads=private$paramlist[["subsampleSize"]])
            print("test2")
            unlink(paste0(private$paramlist[["reportOutput"]],".tmp"))
            print(qcval)
            print(private$paramlist[["reportOutput"]])
            write.table(as.data.frame(qcval),file = private$paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["samInput"]])){
                stop("samInput is required.")
            }


        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["samInput"]]);
            private$checkFileCreatable(private$paramlist[["reportOutput"]]);
        },
        getReportValImp = function(item){
            qcval <- as.list(read.table(file= private$paramlist[["reportOutput"]],header=TRUE))
            if(item == "report"){
                return(data.frame(Item=names(qcval),Value=as.character(qcval)))
            }else{
                return(qcval[[item]])
            }
        },
        getReportItemsImp = function(){
            return(c("report","NRF","PBC1","PBC2","one","two","total","reads"))
        }
    )


)

atacLibComplexQC<-function(atacProc,reportOutput=NULL,samInput=NULL,singleEnd = FALSE,subsample=TRUE,subsampleSize=4*10e6){
    libqc<-LibComplexQC$new(atacProc,reportOutput=reportOutput,samInput=samInput,singleEnd = singleEnd,
                            subsample=subsample,subsampleSize=subsampleSize,editable=FALSE)
    libqc$process()
    invisible(libqc)
}
