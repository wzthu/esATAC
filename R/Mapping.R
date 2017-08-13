BowtieMapping <- R6::R6Class(
  classname = "BowtieMapping",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, fileInput = NULL, Reference = NULL, fileOutput = NULL, In_type = NULL, editable = FALSE){
      super$initialize("BowtieMapping", editable, list(arg1 = atacProc))
      print("MappingInitCall")
      # necessary and unchanged parameters, this should be tested
      private$paramlist[["genome_ref"]] <- Reference
      private$paramlist[["seq_type"]] <- In_type

      # from obj or not
      if(!is.null(atacProc)){
        if(In_type == "single"){
          private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1")
        }else{
          private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1")
          private$paramlist[["fastqInput2"]] <- atacProc$getParam("fastqOutput2")
        }
        if(is.null(fileOutput)){
          file_path <- dirname(private$paramlist[["fastqInput1"]])
          private$paramlist[["samOutput"]] <- paste0(file_path, "/mapping_result.sam", collapse = "")
        }else{
          private$paramlist[["samOutput"]] <- fileOutput
        }
        # file check
        private$checkFileExist(private$paramlist[["fastqInput1"]]);
        private$checkPathExist(private$paramlist[["samOutput"]]);
      }else{
        if(In_type == "single"){
          private$paramlist[["fastqInput1"]] <- fileInput[[1]]
        }else{
          private$paramlist[["fastqInput1"]] <- fileInput[[1]]
          private$paramlist[["fastqInput2"]] <- fileInput[[2]]
        }
        if(is.null(fileOutput)){
          file_path <- dirname(private$paramlist[["fastqInput1"]])
          private$paramlist[["samOutput"]] <- paste0(file_path, "/mapping_result.sam", collapse = "")
        }else{
          private$paramlist[["samOutput"]] <- fileOutput
        }
        # file check
        private$checkFileExist(private$paramlist[["fastqInput1"]]);
        private$checkFileExist(private$paramlist[["fastqInput2"]]);
        private$checkPathExist(private$paramlist[["samOutput"]]);
      }
      private$checkRequireParam()
      print("finishMappingInitCall")
    },

    processing = function(){
      if(!super$processing()){
        return()
      }
      if(private$paramlist[["seq_type"]] == "single"){
        seq_file <- as.vector(unlist(private$paramlist[["fastqInput1"]]))
      }else{
        seq_file <- list(private$paramlist[["fastqInput1"]], private$paramlist[["fastqInput2"]])
      }
      Rbowtie::bowtie(sequences = seq_file, index = private$paramlist[["genome_ref"]], S = TRUE, X = 2000, m = 1,threads = getConfigure("threads"),
                      type = private$paramlist[["seq_type"]], outfile = private$paramlist[["samOutput"]], force = TRUE, strict = TRUE)
      private$setFinish()
    },

    setResultParam = function(samOutput){
      super$setResultParam();
      private$paramlist[["samOutput"]] <- samOutput
    }
  ),

  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["genome_ref"]])){
        stop("Genome reference is required!")
      }
      if(is.null(private$paramlist[["seq_type"]])){
        stop("input file type is required!")
      }
    }
  )

)


Bowtie2Mapping <-R6Class(
    classname = "Bowtie2Mapping",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc,samOutput=NULL, bt2Idx=NULL, fastqInput1=NULL, fastqInput2=NULL,paramList="default",editable=FALSE){
            super$initialize("Bowtie2Mapping",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1");
                private$paramlist[["fastqInput2"]] <- atacProc$getParam("fastqOutput2");
                regexProcName<-sprintf("(fastq|fq|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(fastq|fq)"
            }

            if(!is.null(fastqInput1)){
                private$paramlist[["fastqInput1"]] <- fastqInput1;
            }
            if(!is.null(fastqInput2)){
                private$paramlist[["fastqInput2"]] <- fastqInput2;
            }

            if(is.null(private$paramlist[["fastqInput2"]])){
                private$singleEnd<-TRUE
            }else{
                private$singleEnd<-FALSE
            }

            if(is.null(samOutput)){
                if(!is.null(private$paramlist[["fastqInput1"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput1"]],regexProcName)
                    private$paramlist[["samOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".sam"))
                }
            }else{
                private$paramlist[["samOutput"]] <- samOutput;
            }

            if(is.null(bt2Idx)){
                private$paramlist[["bt2Idx"]]=.obtainConfigure("bt2Idx")
            }else{
                private$paramlist[["bt2Idx"]]<-bt2Idx
            }

            if(.obtainConfigure("threads")>1){
                private$paramlist[["paramList"]]<-c("-p",as.character(.obtainConfigure("threads")))
            }
            if(is.null(paramList)){

            }else if(paramList=="default"){
                private$paramlist[["paramList"]]<-c("--no-discordant","--no-unal","--no-mixed","-X","2000",
                                                    private$paramlist[["paramList"]])
            }else{
                private$paramlist[["paramList"]]<-c(paramList,private$paramlist[["paramList"]])
                rejectp<-"-p|--threads|-x|-1|-2|-U|-S"
                private$checkParam(paramlist,rejectp)
            }
            private$paramValidation()
        }


    ),
    private = list(
        processing = function(){
            private$writeLog("start mapping with parameters:")
            private$writeLog(paste0("bowtie2 index:",private$paramlist[["bt2Idx"]]))
            private$writeLog(paste0("samOutput:",private$paramlist[["samOutput"]]))
            private$writeLog(paste0("fastqInput1:",private$paramlist[["fastqInput1"]]))
            private$writeLog(paste0("fastqInput2:",private$paramlist[["fastqInput2"]]))
            private$writeLog(paste0("other parameters:",paste(private$paramlist[["paramList"]],collapse = " ")))
            .bowtie2_call(bowtie2Index=private$paramlist[["bt2Idx"]],samOutput=private$paramlist[["samOutput"]],
                          fastqInput1=private$paramlist[["fastqInput1"]],fastqInput2=private$paramlist[["fastqInput2"]],
                          paramlist=private$paramlist[["paramList"]])



        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["fastqInput1"]])){
                stop("fastqInput1 is required.")
            }
            if(is.null(private$paramlist[["bt2Idx"]])){
                stop("bt2Idx is required")
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["fastqInput1"]]);
            private$checkFileExist(private$paramlist[["fastqInput2"]]);
            private$checkFileCreatable(private$paramlist[["samOutput"]]);
        }
    )


)


atacBowtie2Mapping <- function(atacProc,samOutput=NULL, bt2Idx=NULL,fastqInput1=NULL, fastqInput2=NULL,paramList="default"){
    atacproc<-Bowtie2Mapping$new(atacProc=atacProc,bt2Idx=bt2Idx,samOutput=samOutput, fastqInput1=fastqInput1,
                                   fastqInput2=fastqInput2,paramList=paramList)
    atacproc$process()
    return(atacproc)
}


