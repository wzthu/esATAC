RGo <- R6::R6Class(
  classname = "RGo",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, gene = NULL, OrgDb = NULL, keytype = NULL, ont = NULL,
                          pvalueCutoff = NULL, pAdjustMethod = NULL, universe = NULL, qvalueCutoff = NULL,
                          readable = NULL, pool = NULL, output = NULL, editable = FALSE){

      super$initialize("RGo", editable, list(arg1=atacProc))
      print("RGoInitCall")

      # necessary and unchanged parameters, this should be tested
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }

      private$paramlist[["gene"]] <- gene
      private$paramlist[["OrgDb"]] <- OrgDb
      private$paramlist[["keytype"]] <- keytype
      private$paramlist[["ont"]] <- ont
      private$paramlist[["pvalueCutoff"]] <- pvalueCutoff
      private$paramlist[["pAdjustMethod"]] <- pAdjustMethod
      private$paramlist[["universe"]] <- universe
      private$paramlist[["qvalueCutoff"]] <- qvalueCutoff
      private$paramlist[["readable"]] <- readable
      private$paramlist[["pool"]] <- pool
      private$paramlist[["output"]] <- output
      private$checkRequireParam()

      print("finishRGoInitCall")
    }, # initialization end

    processing = function(){
      if(!super$processing()){
        return()
      }
      tmp <- clusterProfiler::enrichGO(gene = private$paramlist[["gene"]],
                                       OrgDb = private$paramlist[["OrgDb"]],
                                       keytype = private$paramlist[["keytype"]],
                                       ont = private$paramlist[["ont"]],
                                       pvalueCutoff = private$paramlist[["pvalueCutoff"]],
                                       pAdjustMethod = private$paramlist[["pAdjustMethod"]],
                                       universe = private$paramlist[["universe"]],
                                       qvalueCutoff = private$paramlist[["qvalueCutoff"]],
                                       readable = private$paramlist[["readable"]],
                                       pool = private$paramlist[["pool"]])
      write.table(x = tmp, file = private$paramlist[["output"]], append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
      private$setFinish()
    }, # processing end


    setResultParam = function(output_path){
      private$paramlist[["output"]] <- output_path
    }  # setResultParam end

  ), # public end


  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["gene"]])){
        stop("Parameter gene is required!")
      }
      if(is.null(private$paramlist[["OrgDb"]])){
        stop("Parameter OrgDb is required!")
      }
      if(is.null(private$paramlist[["output"]])){
        stop("Parameter output is required!")
      }
    }
  ) # private end

) # R6 class end
