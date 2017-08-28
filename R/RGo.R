RGo <- R6::R6Class(
  classname = "RGo",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, gene = NULL, OrgDb = NULL, keytype = NULL, ont = NULL,
                          pvalueCutoff = NULL, pAdjustMethod = NULL, universe = NULL, qvalueCutoff = NULL,
                          readable = NULL, pool = NULL, output = NULL, editable = FALSE){
      super$initialize("RGo", editable, list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){
        tmp <- read.table(file = atacProc$getParam("annoOutput"), header = TRUE,
                          sep = "\t", quote = "");
        private$paramlist[["gene"]] <- as.character(base::unique(tmp$geneId))
      }else{
        private$paramlist[["gene"]] <- gene
      }
      private$paramlist[["OrgDb"]] <- OrgDb
      private$paramlist[["output"]] <- output

      # unnecessary parameters
      private$paramlist[["keytype"]] <- keytype
      private$paramlist[["ont"]] <- ont
      private$paramlist[["pvalueCutoff"]] <- pvalueCutoff
      private$paramlist[["pAdjustMethod"]] <- pAdjustMethod
      private$paramlist[["universe"]] <- universe
      private$paramlist[["qvalueCutoff"]] <- qvalueCutoff
      private$paramlist[["readable"]] <- readable
      private$paramlist[["pool"]] <- pool

      # parameter check
      private$paramValidation()
    } # initialization end

  ), # public end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("Database:%s", private$paramlist[["OrgDb"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["output"]]))
      private$writeLog(sprintf("keytype of input gene:%s", private$paramlist[["keytype"]]))
      tmp <- clusterProfiler::enrichGO(gene = private$paramlist[["gene"]],
                                       OrgDb = private$paramlist[["OrgDb"]],
                                       keyType = private$paramlist[["keytype"]],
                                       ont = private$paramlist[["ont"]],
                                       pvalueCutoff = private$paramlist[["pvalueCutoff"]],
                                       pAdjustMethod = private$paramlist[["pAdjustMethod"]],
                                       universe = private$paramlist[["universe"]],
                                       qvalueCutoff = private$paramlist[["qvalueCutoff"]],
                                       readable = private$paramlist[["readable"]],
                                       pool = private$paramlist[["pool"]])
      write.table(x = tmp, file = private$paramlist[["output"]],
                  append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
    }, # processing end


    checkRequireParam = function(){
      if(is.null(private$paramlist[["gene"]])){
        stop("Parameter atacProc or gene is required!")
      }
      if(is.null(private$paramlist[["OrgDb"]])){
        stop("Parameter OrgDb is required!")
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

#' Using clusterProfiler to do GO analysis.
#' @param atacProc Result from function "PeakAnno", extract all the gene ID.
#'
#'
GOAnalysis <- function(atacProc = NULL, gene = NULL, OrgDb = NULL, keytype = "ENTREZID", ont = "MF",
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, qvalueCutoff = 0.2,
                       readable = FALSE, pool = FALSE, output = NULL){
  tmp <- RGo$new(atacProc, gene, OrgDb, keytype, ont, pvalueCutoff,
                 pAdjustMethod , universe, qvalueCutoff, readable, pool, output)
  tmp$process()
  return(tmp)
}
