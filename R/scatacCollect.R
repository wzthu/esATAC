setClass(Class = "scatacCollect",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "scatacCollect",
    definition = function(.Object, prevSteps = list(), ...){
        allparam <- list(...)
        fragmentFile <- allparam[["fragment"]]
        fragmentIdxFile <- allparam[["fragmentIdx"]]
        bamFile <- allparam[["bam"]]
        bamIdxFile <- allparam[["bamIdx"]]
        csvFile <- allparam[["csv"]]

        atacProc <- NULL
        if(length(prevSteps) > 0){
            atacProc <- prevSteps[[1]]
        }

        # necessary parameters
        if(!is.null(atacProc)){
            param(.Object)[["fragmentFile"]] <- getParam(atacProc, "fragOutput")
            param(.Object)[["fragmentIdxFile"]] <- getParam(atacProc, "fragIdxOutput")
            param(.Object)[["bamFile"]] <- getParam(atacProc, "bamOutput")
            param(.Object)[["bamIdxFile"]] <- getParam(atacProc, "bamIdxOutput")
            param(.Object)[["csvFile"]] <- getParam(atacProc, "csvOutput")
        }else{
            param(.Object)[["fragmentFile"]] <- fragmentFile
            param(.Object)[["fragmentIdxFile"]] <- fragmentIdxFile
            param(.Object)[["bamFile"]] <- bamFile
            param(.Object)[["bamIdxFile"]] <- bamIdxFile
            param(.Object)[["csvFile"]] <- csvFile
        }

        if (is.null(param(.Object)[["fragmentFile"]])) {
            stop("Fragment file is missing! Please check parameter 'atacProc' or 'fragment'!")
        }

        if (is.null(param(.Object)[["fragmentIdxFile"]])) {
            fragmentIdxFile_path <- paste0(param(.Object)[["fragmentFile"]], ".tbi")
            mess <- paste0("Fragment file index is not specified, using The following default path:",
                           fragmentIdxFile_path)
            warning(mess)
            param(.Object)[["fragmentIdxFile"]] <- fragmentIdxFile_path
        }

        if (!file.exists(param(.Object)[["csvFile"]])) {
            stop("Single Cell CSV file index is missing!")
        }

        # init output
        output(.Object)[["fragOutput"]] <- param(.Object)[["fragmentFile"]]
        output(.Object)[["fragIdxOutput"]] <- param(.Object)[["fragmentIdxFile"]]
        output(.Object)[["bamOutput"]] <- param(.Object)[["bamFile"]]
        output(.Object)[["bamIdxOutput"]] <- param(.Object)[["bamIdxFile"]]
        output(.Object)[["csvOutput"]] <- param(.Object)[["csvFile"]]
        output(.Object)[["fragOB.rds"]] <- getStepWorkDir(.Object, filename = "fragmentObjects.rds")
        
        property(.Object)[["fragOutput"]] <- output(.Object)[["fragOutput"]]
        property(.Object)[["fragIdxOutput"]] <- output(.Object)[["fragIdxOutput"]]
        property(.Object)[["bamOutput"]] <- output(.Object)[["bamOutput"]]
        property(.Object)[["bamIdxOutput"]] <- output(.Object)[["bamIdxOutput"]]
        property(.Object)[["csvOutput"]] <- output(.Object)[["csvOutput"]]
        property(.Object)[["fragOB.rds"]] <- output(.Object)[["fragOB.rds"]]

        .Object
    }
)


setMethod(
    f = "processing",
    signature = "scatacCollect",
    definition = function(.Object,...){

        mess <- paste0("Now, reading ",
                       param(.Object)[["csvFile"]],
                       "......")

        metadata <- read.csv(file = param(.Object)[["csvFile"]],
                             header = TRUE,
                             row.names = 1)

        print("Processing cell barcode.......")
        cells <- rownames(metadata)

        mess <- paste0("Creating Fragment Object from ",
                       param(.Object)[["fragmentFile"]],
                       "......")
        print(mess)
        fragments <- CreateFragmentObject(param(.Object)[["fragmentFile"]],
                                          cells = cells)

        saveRDS(object = fragments,
                file = output(.Object)[["fragOB.rds"]])

        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "scatacCollect",
    definition = function(.Object, ...){
        report(.Object)[["fragOB.rds"]] <- output(.Object)[["fragOB.rds"]]
        .Object
    }
)


#' @name scatacCollect
#'
#' @title Get Single Cell Pre-processing Information
#'
#' @description
#' Get scATAC-seq pre-processing results and create Fragment Object.
#'
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' ###########################################################.
#' @param fragment A fragment file path.
#' @param fragmentIdx Genome wide annotation databese.
#' @param bam Keytype of input gene.
#' @param bamIdx One of "MF", "BP", and "CC" subontologies.
#' "MF" for molecular function,
#' "BP" for biological process, "CC" for cellular component.
#' @param csv pvalueCutoff.
#' @param ... Additional arguments, currently unused.
#'
#' @return An invisible \code{\link{ATACProc-class}} object scalar.
#'
#' @author Wei Zhang
#'
#' @examples
#' print(123)


setGeneric("scCollect", function(atacProc, fragment = NULL,
                                         fragmentIdx = NULL,
                                         bam = NULL, bamIdx = NULL,
                                         csv = NULL, ...) standardGeneric("scCollect"))

#' @rdname scatacCollect
#' @aliases scCollect
#' @export
setMethod(
    f = "scCollect",
    signature = "ATACProc",
    definition = function(atacProc, fragment = NULL, fragmentIdx = NULL,
                          bam = NULL, bamIdx = NULL,
                          csv = NULL, ...){
        allpara <- c(list(Class = "scatacCollect", prevSteps = list(atacProc)),
                     as.list(environment()), list(...))
        step <- do.call(new, allpara)
        invisible(step)
    }
)

#' @rdname scatacCollect
#' @aliases sccollect
#' @export
sccollect <- function(fragment = NULL, fragmentIdx = NULL,
                       bam = NULL, bamIdx = NULL,
                       csv = NULL, ...){

    allpara <- c(list(Class = "scatacCollect", prevSteps = list()),
                 as.list(environment()),
                 list(...))

    step <- do.call(new, allpara)

    invisible(step)
}
