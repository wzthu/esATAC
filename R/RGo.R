setClass(Class = "RGo",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "RGo",
    definition = function(.Object, atacProc, ..., gene = NULL, OrgDb = NULL, keytype = NULL, ont = NULL,
                          pvalueCutoff = NULL, pAdjustMethod = NULL, universe = NULL, qvalueCutoff = NULL,
                          readable = NULL, pool = NULL, goOutput = NULL, editable = FALSE){
        .Object <- init(.Object, "RGo", editable, list(arg1 = atacProc))

        if(!is.null(atacProc)){
            tmp <- read.table(file = getParam(atacProc, "annoOutput.df"), header = TRUE,
                              sep = "\t", quote = "")
            # .Object@paramlist[["gene"]] <- as.character(base::unique(tmp$geneId))
            .Object@paramlist[["gene"]] <- as.character(tmp$geneId)
            regexProcName <- sprintf("(df|%s)", getProcName(atacProc))
        }else{
            .Object@paramlist[["gene"]] <- gene
        }
        if(is.null(OrgDb)){
            .Object@paramlist[["OrgDb"]] <- .obtainConfigure("annoDb")
        }else{
            .Object@paramlist[["OrgDb"]] <- OrgDb
        }

        if(is.null(goOutput)){
            if(!is.null(atacProc)){
                prefix <- getBasenamePrefix(.Object, getParam(atacProc, "annoOutput.df"), regexProcName)
                .Object@paramlist[["goOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                             paste(prefix, ".", getProcName(.Object), ".df", sep = ""))
            }else{
                .Object@paramlist[["goOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                             paste("GOanalysis.", getProcName(.Object), ".df", sep = ""))
            }
        }else{
            name_split <- unlist(base::strsplit(x = goOutput, split = ".", fixed = TRUE))
            suffix <- tail(name_split, 1)
            if(suffix == "df"){
                .Object@paramlist[["goOutput"]] <- goOutput
            }else{
                .Object@paramlist[["goOutput"]] <- paste(goOutput, ".df", sep = "")
            }
        }

        .Object@paramlist[["keytype"]] <- keytype
        .Object@paramlist[["ont"]] <- ont
        .Object@paramlist[["pvalueCutoff"]] <- pvalueCutoff
        .Object@paramlist[["pAdjustMethod"]] <- pAdjustMethod
        .Object@paramlist[["universe"]] <- universe
        .Object@paramlist[["qvalueCutoff"]] <- qvalueCutoff
        .Object@paramlist[["readable"]] <- readable
        .Object@paramlist[["pool"]] <- pool

        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RGo",
    definition = function(.Object,...){
        .Object <- writeLog(.Object,paste0("processing file:"))
        .Object <- writeLog(.Object,sprintf("Database:%s",.Object@paramlist[["OrgDb"]]))
        .Object <- writeLog(.Object,sprintf("destination:%s",.Object@paramlist[["goOutput"]]))
        .Object <- writeLog(.Object,sprintf("dkeytype of input gene:%s",.Object@paramlist[["keytype"]]))
        tmp <- clusterProfiler::enrichGO(gene = .Object@paramlist[["gene"]],
                                         OrgDb = .Object@paramlist[["OrgDb"]],
                                         keyType = .Object@paramlist[["keytype"]],
                                         ont = .Object@paramlist[["ont"]],
                                         pvalueCutoff = .Object@paramlist[["pvalueCutoff"]],
                                         pAdjustMethod = .Object@paramlist[["pAdjustMethod"]],
                                         universe = .Object@paramlist[["universe"]],
                                         qvalueCutoff = .Object@paramlist[["qvalueCutoff"]],
                                         readable = .Object@paramlist[["readable"]],
                                         pool = .Object@paramlist[["pool"]])
        write.table(x = tmp, file = .Object@paramlist[["goOutput"]],
                    append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "RGo",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["gene"]])){
            stop("Parameter atacProc or gene is required!")
        }
        if(is.null(.Object@paramlist[["OrgDb"]])){
            stop("Parameter OrgDb is required!")
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "RGo",
    definition = function(.Object,...){
        checkPathExist(.Object, .Object@paramlist[["goOutput"]])
    }
)


setMethod(
    f = "getReportValImp",
    signature = "RGo",
    definition = function(.Object, item){
        if(item == "goOutput"){
            return(.Object@paramlist[["goOutput"]])
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "RGo",
    definition = function(.Object){
        return(c("goOutput"))
    }
)


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
#' library(clusterProfiler)
#' data(geneList)
#' geneId <- names(geneList)[1:100]
#' #goanalysis(gene = geneId, OrgDb = 'org.Hs.eg.db')
#'
#'
#' @references Guangchuang Yu., Li-Gen Wang, Yanyan Han, Qing-Yu He.
#' clusterProfiler: an R package for comparing biological themes among gene
#' clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287
#' @seealso
#' \code{\link{atacPeakAnno}}
#' \link[clusterProfiler]{enrichGO} function enrichGO in package
#' "clusterProfiler"
#' @importFrom clusterProfiler enrichGO
#' @name atacGOAnalysis
#' @export
#' @docType methods
#' @rdname atacGOAnalysis-methods
setGeneric("atacGOAnalysis",function(atacProc = NULL, gene = NULL, OrgDb = NULL, keytype = "ENTREZID", ont = "MF",
                                     pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, qvalueCutoff = 0.2,
                                     readable = FALSE, pool = FALSE, goOutput = NULL, ...) standardGeneric("atacGOAnalysis"))
#' @rdname atacGOAnalysis-methods
#' @aliases atacGOAnalysis
setMethod(
    f = "atacGOAnalysis",
    signature = "ATACProc",
    definition = function(atacProc = NULL, gene = NULL, OrgDb = NULL, keytype = "ENTREZID", ont = "MF",
                          pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, qvalueCutoff = 0.2,
                          readable = FALSE, pool = FALSE, goOutput = NULL, ...){
        atacproc <- new(
            "RGo",
            atacProc = atacProc,
            gene = gene,
            OrgDb = OrgDb,
            keytype = keytype,
            ont = ont,
            pvalueCutoff = pvalueCutoff,
            pAdjustMethod = pAdjustMethod,
            universe = universe,
            qvalueCutoff = qvalueCutoff,
            readable = readable,
            pool = pool,
            goOutput = goOutput)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname atacGOAnalysis-methods
#' @export
goanalysis <- function(gene, OrgDb = NULL, keytype = "ENTREZID", ont = "MF",
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, qvalueCutoff = 0.2,
                       readable = FALSE, pool = FALSE, goOutput = NULL, ...){
    atacproc <- new(
        "RGo",
        atacProc = NULL,
        gene = gene,
        OrgDb = OrgDb,
        keytype = keytype,
        ont = ont,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        universe = universe,
        qvalueCutoff = qvalueCutoff,
        readable = readable,
        pool = pool,
        goOutput = goOutput)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
