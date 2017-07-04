Mapping <- R6::R6Class(
  classname = "Mapping",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, fileInput = NULL, Reference = NULL, fileOutput = NULL, In_type = NULL, editable = FALSE){
      super$initialize("Mapping", editable, list(arg1 = atacProc))
      print("MappingInitCall")
      # necessary and unchanged parameters, this should be tested
      private$paramlist[["genome_ref"]] <- Reference
      private$paramlist[["seq_type"]] <- In_typ	
      wsq	
      ed
       
       23
        r
        4q23
         rt
         q3
          t
          q3 
          tq
          3 t
          qw
           t
           wdf
           sax
           f
            wq
            r 
            23	
            r
            
            
            
            rq23
            	qer23r3	rwqefr
            	weqa
            	f
            	qwa
            	 
            	 r 
            	 2
            	 	q 
            	 	
            	 	
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

    setResultParam = function(samFilePath){
      super$setResultParam();
      private$paramlist[["samOutput"]] <- samFilePath
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
