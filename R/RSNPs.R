RSNPs <- R6::R6Class(
  classname = "RSNPs",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, snp.regions.file = NULL, bio.features.loc = NULL,
                          built.in.biofeatures = NULL, par.threads = NULL,
                          verbose = NULL, method.p = NULL, search.window = NULL,
                          output = NULL, editable = FALSE){
      super$initialize("RSNPs",editable,list(arg1=atacProc))
      print("RSNPsInitCall")

      # necessary and unchanged parameters, this should be tested
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }
      private$paramlist[["snp.regions.file"]] <- snp.regions.file
      private$paramlist[["bio.features.loc"]] <- bio.features.loc
      private$paramlist[["built.in.biofeatures"]] <- built.in.biofeatures
      private$paramlist[["par.threads"]] <- par.threads
      private$paramlist[["verbose"]] <- verbose
      private$paramlist[["method.p"]] <- method.p
      private$paramlist[["search.window"]] <- search.window
      private$paramlist[["output"]] <- output
      # check parameter
      private$checkRequireParam()
      print("finishRSNPsInitCall")
    }, # initialization end

    processing = function(){
      if(!super$processing()){
        return()
      }
      tmp <- FunciSNP::getFSNPs(snp.regions.file = private$paramlist[["snp.regions.file"]],
                                                          bio.features.loc = private$paramlist[["bio.features.loc"]],
                                                          built.in.biofeatures = private$paramlist[["built.in.biofeatures"]],
                                                          par.threads = private$paramlist[["par.threads"]],
                                                          verbose = private$paramlist[["verbose"]],
                                                          method.p = private$paramlist[["method.p"]],
                                                          search.window = private$paramlist[["search.window"]])
      saveRDS(tmp, private$paramlist[["output"]])
      private$setFinish()
    }, # processing end

    setResultParam = function(){
      print("This function is not using!")
    }  # setResultParam end

  ), # piblic end


  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["snp.regions.file"]])){
        stop("Parameter snp.regions.file is required!")
      }
      if(is.null(private$paramlist[["bio.features.loc"]])){
        stop("Parameter bio.features.loc is required!")
      }
      if(is.null(private$paramlist[["output"]])){
        stop("Parameter output is required!")
      }
    }
  ) # private end



) # class end
