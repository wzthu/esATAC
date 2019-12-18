setClass(Class = "RGo",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "RGo",
    definition = function(.Object, prevSteps = list(), ...){
        allparam <- list(...)
        gene <- allparam[["gene"]]
        OrgDb <- allparam[["OrgDb"]]
        keytype <- allparam[["keytype"]]
        ont <- allparam[["ont"]]
        pvalueCutoff <- allparam[["pvalueCutoff"]]
        pAdjustMethod <- allparam[["pAdjustMethod"]]
        universe <- allparam[["universe"]]
        qvalueCutoff <- allparam[["qvalueCutoff"]]
        readable <- allparam[["readable"]]
        pool <- allparam[["pool"]]
        goOutput <- allparam[["goOutput"]]
       
        atacProc <- NULL
        if(length(prevSteps) > 0){
            atacProc <- prevSteps[[1]]
        }

        if(!is.null(atacProc)){
            tmp <- read.table(file = getParam(atacProc, "annoOutput.txt"), header = TRUE,
                              sep = "\t", quote = "")
            # .Object@paramlist[["gene"]] <- as.character(base::unique(tmp$geneId))
            param(.Object)[["gene"]] <- as.character(tmp$geneId)
        }else{
            param(.Object)[["gene"]] <- gene
        }
        if(is.null(OrgDb)){
            param(.Object)[["OrgDb"]] <- getRefRc("annoDb")
        }else{
            param(.Object)[["OrgDb"]] <- OrgDb
        }

        if(is.null(goOutput)){
            if(!is.null(atacProc)){
                output(.Object)[["goOutput"]] <- getAutoPath(.Object, getParam(atacProc, "annoOutput.txt"),
                                                             "txt", "txt")
            }else{
                output(.Object)[["goOutput"]] <- file.path(getStepWorkDir(.Object),
                                                             "GOanalysis.txt")
            }
        }else{
            output(.Object)[["goOutput"]] <- addFileSuffix(goOutput, ".txt")
        }

        param(.Object)[["keytype"]] <- keytype
        param(.Object)[["ont"]] <- ont
        param(.Object)[["pvalueCutoff"]] <- pvalueCutoff
        param(.Object)[["pAdjustMethod"]] <- pAdjustMethod
        param(.Object)[["universe"]] <- universe
        param(.Object)[["qvalueCutoff"]] <- qvalueCutoff
        param(.Object)[["readable"]] <- readable
        param(.Object)[["pool"]] <- pool

        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RGo",
    definition = function(.Object,...){
        tmp <- clusterProfiler::enrichGO(gene = param(.Object)[["gene"]],
                                         OrgDb = param(.Object)[["OrgDb"]],
                                         keyType = param(.Object)[["keytype"]],
                                         ont = param(.Object)[["ont"]],
                                         pvalueCutoff = param(.Object)[["pvalueCutoff"]],
                                         pAdjustMethod = param(.Object)[["pAdjustMethod"]],
                                         universe = param(.Object)[["universe"]],
                                         qvalueCutoff = param(.Object)[["qvalueCutoff"]],
                                         readable = param(.Object)[["readable"]],
                                         pool = param(.Object)[["pool"]])
        write.table(x = tmp, file = output(.Object)[["goOutput"]],
                    append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "RGo",
    definition = function(.Object, ...){
        report(.Object)$goOutput <- output(.Object)[["goOutput"]]
        .Object
    }
)


#' @name RGo
#' @title Gene Ontology Analysis
#' @description
#' Ranking functional groups based on a set of genes. For more information,
#' please see \link[clusterProfiler]{enrichGO}.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakAnno}}.
#' @param gene A vector of entrez gene id.
#' @param OrgDb Genome wide annotation databese.
#' @param keytype Keytype of input gene.
#' @param ont One of "MF", "BP", and "CC" subontologies.
#' "MF" for molecular function,
#' "BP" for biological process, "CC" for cellular component.
#' @param pvalueCutoff pvalueCutoff.
#' @param pAdjustMethod One of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' @param universe Background genes.
#' @param qvalueCutoff qvalue cutoff.
#' @param readable whether mapping gene ID to gene Name.
#' @param pool If ont=’ALL’, whether pool 3 GO sub-ontologies.
#' @param goOutput \code{Character} scalar.
#' Output file path. Defult:in the same folder as your input file with the
#' suffix "df".
#' @param ... Additional arguments, currently unused.
#' @details This function using \link[clusterProfiler]{enrichGO} to do GO
#' analysis but fixed some parameters. If atacProc is not NULL, it will read
#' the gene ID from the output of \code{\link{atacPeakAnno}}.
#' @return An invisible \code{\link{ATACProc-class}} object scalar.
#' @author Wei Zhang
#' @examples
#'
#' \dontrun{
#' library(org.Hs.eg.db)
#' # generate simulated geneID
#' geneId <- as.character(sample(seq(10000), 100))
#' goanalysis(gene = geneId, OrgDb = 'org.Hs.eg.db')
#' }
#'
#' @references Guangchuang Yu., Li-Gen Wang, Yanyan Han, Qing-Yu He.
#' clusterProfiler: an R package for comparing biological themes among gene
#' clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287
#' @seealso
#' \code{\link{atacPeakAnno}}
#' \link[clusterProfiler]{enrichGO} function enrichGO in package
#' "clusterProfiler"
#' @importFrom clusterProfiler enrichGO


setGeneric("atacGOAnalysis",function(atacProc, gene = NULL, OrgDb = NULL, keytype = "ENTREZID", ont = "MF",
                                     pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, qvalueCutoff = 0.2,
                                     readable = FALSE, pool = FALSE, goOutput = NULL, ...) standardGeneric("atacGOAnalysis"))

#' @rdname RGo
#' @aliases atacGOAnalysis
#' @export
setMethod(
    f = "atacGOAnalysis",
    signature = "ATACProc",
    definition = function(atacProc, gene = NULL, OrgDb = NULL, keytype = "ENTREZID", ont = "MF",
                          pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, qvalueCutoff = 0.2,
                          readable = FALSE, pool = FALSE, goOutput = NULL, ...){
        allpara <- c(list(Class = "RGo", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname RGo
#' @aliases goanalysis
#' @export
goanalysis <- function(gene, OrgDb = NULL, keytype = "ENTREZID", ont = "MF",
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, qvalueCutoff = 0.2,
                       readable = FALSE, pool = FALSE, goOutput = NULL, ...){
    allpara <- c(list(Class = "RGo", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
