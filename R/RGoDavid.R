RGoDavid <- R6::R6Class(
  classname = "RGoDavid",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, gene = NULL, idType = NULL,annotation = NULL,
                          pvalueCutoff = NULL, pAdjustMethod = NULL,qvalueCutoff = NULL,
                          species = NULL, david.user = NULL, output = NULL, editable = FALSE){
      super$initialize("RGoDavid",editable,list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }

      private$paramlist[["gene"]] <- gene
      private$paramlist[["idType"]] <- idType
      private$paramlist[["annotation"]] <- annotation
      private$paramlist[["pvalueCutoff"]] <- pvalueCutoff
      private$paramlist[["pAdjustMethod"]] <- pAdjustMethod
      private$paramlist[["qvalueCutoff"]] <- qvalueCutoff
      private$paramlist[["species"]] <- species
      private$paramlist[["david.user"]] <- david.user
      private$paramlist[["output"]] <- output
      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("destination:%s", private$paramlist[["output"]]))
      tmp <- clusterProfiler::enrichDAVID(gene = private$paramlist[["gene"]],
                                          idType = private$paramlist[["idType"]],
                                          annotation = private$paramlist[["annotation"]],
                                          pvalueCutoff = private$paramlist[["pvalueCutoff"]],
                                          pAdjustMethod = private$paramlist[["pAdjustMethod"]],
                                          qvalueCutoff = private$paramlist[["qvalueCutoff"]],
                                          species = private$paramlist[["species"]],
                                          david.user = private$paramlist[["david.user"]])
      write.table(x = tmp, file = private$paramlist[["output"]], append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["gene"]])){
        stop("Parameter gene is required!")
      }
      if(is.null(private$paramlist[["david.user"]])){
        stop("Parameter david.user is required!")
      }
      if(is.null(private$paramlist[["output"]])){
        stop("Parameter output is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkPathExist(private$paramlist[["output"]])
    } # checkAllPath end

  ) # private end

) # R6 class end

#' Using clusterProfiler to do GO analysis(David method).
#'
#'
#'
GODavid <- function(atacProc = NULL, gene, idType = "ENTREZ_GENE_ID",
                    annotation = "GOTERM_BP_FAT", pvalueCutoff = 0.05,
                    pAdjustMethod = "BH", qvalueCutoff = 0.2, species = NA, david.user = NULL, output = NULL){
  tmp <- RGoDavid$new(atacProc, gene, idType, annotation, pvalueCutoff,
                      pAdjustMethod, qvalueCutoff, species, david.user, output)
  tmp$process()
  return(tmp)
}
