RGoDavid <- R6::R6Class(
  classname = "RGoDavid",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, gene = NULL, idType = NULL, listType = NULL,
                          annotation = NULL, pvalueCutoff = NULL, pAdjustMethod = NULL,
                          qvalueCutoff = NULL, species = NULL, david.user = NULL, output = NULL, editable = FALSE){

      super$initialize("RGoDavid",editable,list(arg1=atacProc))
      print("RGoDavidInitCall")

      # necessary and unchanged parameters, this should be tested
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }

      private$paramlist[["gene"]] <- gene
      private$paramlist[["idType"]] <- idType
      private$paramlist[["listType"]] <- listType
      private$paramlist[["annotation"]] <- annotation
      private$paramlist[["pvalueCutoff"]] <- pvalueCutoff
      private$paramlist[["pAdjustMethod"]] <- pAdjustMethod
      private$paramlist[["qvalueCutoff"]] <- qvalueCutoff
      private$paramlist[["species"]] <- species
      private$paramlist[["david.user"]] <- david.user
      private$paramlist[["output"]] <- output
      private$checkRequireParam()

      print("finishRGoDavidInitCall")
    }, # initialization end

    processing = function(){
      super$processing()
      tmp <- clusterProfiler::enrichDAVID(gene = private$paramlist[["gene"]],
                                                                        idType = private$paramlist[["idType"]],
                                                                        listType = private$paramlist[["listType"]],
                                                                        annotation = private$paramlist[["annotation"]],
                                                                        pvalueCutoff = private$paramlist[["pvalueCutoff"]],
                                                                        pAdjustMethod = private$paramlist[["pAdjustMethod"]],
                                                                        qvalueCutoff = private$paramlist[["qvalueCutoff"]],
                                                                        species = private$paramlist[["species"]],
                                                                        david.user = private$paramlist[["david.user"]])
      write.table(x = tmp, file = private$paramlist[["output"]], append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
      private$finish <- TRUE
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
      if(is.null(private$paramlist[["david.user"]])){
        stop("Parameter david.user is required!")
      }
      if(is.null(private$paramlist[["output"]])){
        stop("Parameter output is required!")
      }
    }
  ) # private end

) # R6 class end
