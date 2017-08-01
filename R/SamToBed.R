SamToBed <- R6::R6Class(
  classname = "SamToBed",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, merge = TRUE, posOffset = +4, negOffset= -5, chrFilterList= NULL,
                          samInput = NULL, bedOutput = NULL, sortBed = TRUE, uniqueBed = TRUE,
                          minFregLen = 0,maxFregLen = 100, saveExtLen = FALSE, editable=FALSE){
      super$initialize("SamToBed",editable,list(arg1=atacProc))
      print("SamToBedInitCall")
      if(!is.null(atacProc)){
          private$paramlist[["samInput"]] <- atacProc$getParam("samOutput");
      }
      if(!is.null(samInput)){
          private$paramlist[["samInput"]] <- samInput;
      }
      if(is.null(bedOutput)){
          private$paramlist[["bedOutput"]] <- paste0(private$paramlist[["samInput"]],".bed")
      }else{
          private$paramlist[["bedOutput"]] <- bedOutput;
      }

      private$paramlist[["merge"]] <- merge;
      private$paramlist[["posOffset"]] <- posOffset;
      private$paramlist[["negOffset"]] <- negOffset;
      private$paramlist[["filterList"]] <- chrFilterList;
      private$paramlist[["sortBed"]] <- sortBed
      private$paramlist[["uniqueBed"]] <- uniqueBed
      private$paramlist[["minFregLen"]] <- minFregLen
      private$paramlist[["maxFregLen"]] <- maxFregLen
      private$paramlist[["saveExtLen"]] <- saveExtLen

      private$checkFileExist(private$paramlist[["samInput"]]);
      private$checkFileCreatable(private$paramlist[["bedOutput"]]);
      private$checkRequireParam();
      print("finishMappingInitCall")
    },

    processing = function(){
      super$processing()
      if(private$paramlist[["merge"]]){
          .sam2bed_merge_call(samfile = private$paramlist[["samInput"]], bedfile = private$paramlist[["bedOutput"]],
                              posOffset = private$paramlist[["posOffset"]], negOffset = private$paramlist[["negOffset"]],
                              sortBed = private$paramlist[["sortBed"]],uniqueBed = private$paramlist[["uniqueBed"]],
                              filterList = private$paramlist[["filterList"]],minFregLen = private$paramlist[["minFregLen"]],
                              maxFregLen = private$paramlist[["maxFregLen"]],saveExtLen = private$paramlist[["saveExtLen"]] )
      }else{
          .sam2bed_call(samfile = private$paramlist[["SamInput"]], bedfile = private$paramlist[["BedOutput"]])
      }
      private$finish <- TRUE
    },

    setResultParam = function(bedOutput){
      super$setResultParam();
      private$paramlist[["bedOutput"]] <- bedOutput
    }

  ),

  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
        if(is.null(private$paramlist[["samInput"]])){
            stop(paste("samInput is requied"));
        }
    }
  )


)
