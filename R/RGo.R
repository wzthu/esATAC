RGo <- R6::R6Class(
  classname = "RGo",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, gene = NULL, OrgDb = NULL, keytype = NULL, ont = NULL,
                          pvalueCutoff = NULL, pAdjustMethod = NULL, universe = NULL, qvalueCutoff = NULL,
                          readable = NULL, pool = NULL, goOutput = NULL, editable = FALSE){
      super$initialize("RGo", editable, list(arg1 = atacProc))

      # necessary parameters
      if(!is.null(atacProc)){
        tmp <- read.table(file = atacProc$getParam("annoOutput.df"), header = TRUE,
                          sep = "\t", quote = "")
        private$paramlist[["gene"]] <- as.character(base::unique(tmp$geneId))
        regexProcName <- sprintf("(df|%s)", atacProc$getProcName())
      }else{
        private$paramlist[["gene"]] <- gene
      }
      private$paramlist[["OrgDb"]] <- OrgDb

      # unnecessary parameters
      if(is.null(goOutput)){
        if(!is.null(atacProc)){
          prefix <- private$getBasenamePrefix(atacProc$getParam("annoOutput.df"), regexProcName)
          private$paramlist[["goOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                        paste(prefix, ".", self$getProcName(), ".df", sep = ""))
        }else{
          private$paramlist[["goOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                        paste("GOanalysis.", self$getProcName(), ".df", sep = ""))
        }
      }else{
        name_split <- unlist(base::strsplit(x = goOutput, split = ".", fixed = TRUE))
        suffix <- tail(name_split, 1)
        if(suffix == "df"){
          private$paramlist[["goOutput"]] <- goOutput
        }else{
          private$paramlist[["goOutput"]] <- paste(goOutput, ".df", sep = "")
        }
      }

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
      private$writeLog(sprintf("destination:%s", private$paramlist[["goOutput"]]))
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
      write.table(x = tmp, file = private$paramlist[["goOutput"]],
                  append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
    }, # processing end


    checkRequireParam = function(){
      if(is.null(private$paramlist[["gene"]])){
        stop("Parameter atacProc or gene is required!")
      }
      if(is.null(private$paramlist[["OrgDb"]])){
        stop("Parameter OrgDb is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkPathExist(private$paramlist[["goOutput"]])
    } # checkAllPath end

  ) # private end

) # R6 class end

#' Using clusterProfiler to do GO analysis.
#' @param atacProc Result from function "PeakAnno", extract all the gene ID.
#' @param ont One of "MF", "BP", and "CC" subontologies. "MF" for molecular function,
#' "BP" for biological process, "CC" for cellular component.
#'
GOAnalysis <- function(atacProc = NULL, gene = NULL, OrgDb = NULL, keytype = "ENTREZID", ont = "MF",
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, qvalueCutoff = 0.2,
                       readable = FALSE, pool = FALSE, goOutput = NULL){
  tmp <- RGo$new(atacProc, gene, OrgDb, keytype, ont, pvalueCutoff,
                 pAdjustMethod , universe, qvalueCutoff, readable, pool, goOutput)
  tmp$process()
  invisible(tmp)
}
