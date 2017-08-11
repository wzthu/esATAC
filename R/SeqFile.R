SeqFile <-R6Class(
    classname = "SeqFile",
    inherit = BaseProc,
    public = list(
        initialize = function(fastqInput1, fastqInput2=NULL,editable=FALSE){
            super$initialize("SeqFile",editable,list())
            print("SeqFileInitCall")
            private$paramlist[["fastqInput1"]] <- fastqInput1
            private$paramlist[["fastqInput2"]] <- fastqInput2
            private$checkFileExist(private$paramlist[["fastqInput1"]]);
            private$checkFileExist(private$paramlist[["fastqInput2"]]);
            private$paramlist[["fastqOutput1"]] <- fastqInput1
            private$paramlist[["fastqOutput2"]] <- fastqInput2
            private$checkRequireParam();
            private$finish<-TRUE;
        },
        processing = function(){
            super$processing()
            if(is.null(private$paramlist[["fastqInput2"]])){
                print("Single end data:")
                print(private$paramlist[["fastqInput1"]])
            }else{
                print("Paired end data:")
                print(private$paramlist[["fastqInput1"]])
                print(private$paramlist[["fastqInput2"]])
            }
        },
        setResultParam = function(fastqInput1, fastqInput2=NULL){
            super$setResultParam();
            private$paramlist[["fastqInput1"]] <- fastqInput1
            private$paramlist[["fastqInput2"]] <- fastqInput2
            private$checkFileExist(private$paramlist[["fastqInput1"]])
            private$checkFileExist(private$paramlist[["fastqInput2"]])
            private$paramlist[["fastqOutput1"]] <- fastqInput1
            private$paramlist[["fastqOutput2"]] <- fastqInput2
            }
        ),
    private = list(
        checkRequireParam = function(){
            if(private$editable){
                return();
            }
            if(is.null(private$paramlist[["fastqInput1"]])){
                stop("fastqInput1 is required.")
            }
        }
    )

)
