BedOverlaps <- R6::R6Class(
  classname = "BedOverlaps",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, BedInput1 = NULL, BedInput2 = NULL, Output = NULL,
                          n.col = NULL, editable=FALSE){
      super$initialize("BedOverlaps",editable,list(arg1=atacProc))
      print("BedOverlapsInitCall")
      if(!is.null(atacProc)){
        # Not using now
      }else{ # atacProc is NULL
        private$paramlist[["BedInput1"]] <- BedInput1
        private$paramlist[["BedInput2"]] <- BedInput2
        if(is.null(Output)){
          tmp_path <- dirname(private$paramlist[["BedInput1"]])
          private$paramlist[["Output"]] <- paste0(tmp_path, "/Output.bede", collapse = "")
        }else{
          private$paramlist[["Output"]] <- Output
        }
        private$paramlist[["col_num"]] <- n.col
      }
      # parameter check
      private$checkRequireParam();
      private$checkFileExist(private$paramlist[["BedInput1"]]);
      private$checkFileExist(private$paramlist[["BedInput2"]]);
      private$checkPathExist(private$paramlist[["Output"]]);
      print("finishBedOverlapsInitCall")
    }, # initialization end

    processing = function(){
      super$processing()
      BedIntersect(Input1 = private$paramlist[["BedInput1"]], Input2 = private$paramlist[["BedInput2"]],
                   bed_output = private$paramlist[["Output"]], n.col = private$paramlist[["col_num"]])
      private$finish <- TRUE
    }, # processing end

    setResultParam = function(Output){
      super$setResultParam();
      private$paramlist[["Output"]] <- Output
    }  # setResultParam end

  ), # public end


  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["BedInput1"]])){
        stop("Parameter BedInput1 is requied!");
      }
      if(is.null(private$paramlist[["BedInput2"]])){
        stop("Parameter BedInput2 is requied!");
      }
      if(is.null(private$paramlist[["col_num"]])){
        stop("Parameter n.col is requied!");
      }
    }

  ) # private end



) #R6 class end
