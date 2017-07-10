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
      super$processing()
      if(private$paramlist[["seq_type"]] == "single"){
        seq_file <- as.vector(unlist(private$paramlist[["fastqInput1"]]))
      }else{
        seq_file <- list(private$paramlist[["fastqInput1"]], private$paramlist[["fastqInput2"]])
      }
      Rbowtie::bowtie(sequences = seq_file, index = private$paramlist[["genome_ref"]], S = TRUE, X = 2000, m = 1,threads = getConfigure("threads"),
                      type = private$paramlist[["seq_type"]], outfile = private$paramlist[["samOutput"]], force = TRUE, strict = TRUE)
      private$finish <- TRUE
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
        initialize = function(atacProc,bowtie2Index=NULL,samOutput=NULL,threads=NULL,
                              fastqInput1=NULL, fastqInput2=NULL,editable=FALSE){
            super$initialize("Bowtie2Mapping",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1");
                private$paramlist[["fastqInput2"]] <- atacProc$getParam("fastqOutput2");
            }

            if(!is.null(fastqInput1)){
                private$paramlist[["fastqInput1"]] <- fastqInput1;
            }
            if(!is.null(fastqInput2)){
                private$paramlist[["fastqInput2"]] <- fastqInput2;
            }

            if(is.null(samOutput)){
                private$paramlist[["samOutput"]] <- paste(private$paramlist[["fastqInput1"]],".sam",sep="");
            }else{
                private$paramlist[["samOutput"]] <- samOutput;
            }

            if(is.null(bowtie2Index)){
                private$paramlist[["bowtie2Index"]] <- getConfigure("bt2idx");
            }else{
                private$paramlist[["bowtie2Index"]] <- bowtie2Index;
            }

            if(is.null(threads)){
                private$paramlist[["threads"]] <- getConfigure("threads");
            }else{
                private$paramlist[["threads"]] <- as.numeric(threads)
            }

            private$checkFileExist(private$paramlist[["fastqInput1"]]);
            private$checkFileExist(private$paramlist[["fastqInput2"]]);
            private$checkPathExist(private$paramlist[["samOutput"]]);
            private$checkPathExist(private$paramlist[["bowtie2Index"]]);
            private$checkRequireParam();
        },
        processing = function(){
            super$processing()

            if(is.null(private$paramlist[["fastqInput2"]])){
                stop("Single-end-mapping does not support yet!")
            }else{
                .bowtie2_paired_call(bowtie2Index=private$paramlist[["bowtie2Index"]],
                                     samOutput=private$paramlist[["samOutput"]],
                                     fastqInput1=private$paramlist[["fastqInput1"]],
                                     fastqInput2=private$paramlist[["fastqInput2"]],
                                     threads=private$paramlist[["threads"]])
            }

            private$finish <- TRUE
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
            if(is.null(private$paramlist[["fastqInput1"]])){
                stop("fastqInput1 is required.")
            }
            if(is.null(private$paramlist[["bowtie2Index"]])){
                stop("bowtieIdx")
            }
        }
    )


)


