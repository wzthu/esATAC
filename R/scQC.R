setClass(Class = "SCQC",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "SCQC",
    definition = function(.Object, prevSteps = list(), ...){
        # necessary parameters
        allparam <- list(...)
        print(1)
        param(.Object)[["n"]] <- allparam[["n"]]
        param(.Object)[["gene.annotation"]] <- allparam[["gene.annotation"]]
        param(.Object)[["peak"]] <- allparam[["peak"]]
        param(.Object)[["blacklist"]] <- allparam[["blacklist"]]
        print(2)
        if (is.null(param(.Object)[["gene.annotation"]])) {
            stop("Please input gene.annotation!")
        }
        
        if (is.null(param(.Object)[["peak"]])) {
            stop("Please input peak!")
        }
        
        if (is.null(param(.Object)[["blacklist"]])) {
            stop("Please input blacklist!")
        }
        
        print(3)
        atacProc <- NULL
        if(length(prevSteps) > 0){
            atacProc <- prevSteps[[1]]
        }
        
        print(4)
        print(prevSteps)
        print(atacProc)
        if(!is.null(atacProc)){
            stop("The upstream of this step must from atacSCCollect!")
        }else{
            input(.Object)[["fragOB.rds"]] <- output(atacProc)[["fragOB.rds"]]
        }
        
        # init output
        print(5)
        
        output(.Object)[["fragInPeaks.rds"]] <- getStepWorkDir(.Object, filename = "fragInPeaks.rds")
        output(.Object)[["fragInBlacklist.rds"]] <- getStepWorkDir(.Object, filename = "fragInBlacklist.rds")
        output(.Object)[["tssQC.rds"]] <- getStepWorkDir(.Object, filename = "tssQC.rds")
        output(.Object)[["nucleosomeQC.rds"]] <- getStepWorkDir(.Object, filename = "nucleosomeQC.rds")
        
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "SCQC",
    definition = function(.Object,...){
        # reading fragment file
        print("Reading fragment object......")
        fragment <- readRDS(input(.Object)[["fragOB.rds"]])
        
        
        # nucleosome QC
        print("Now, processing nucleosome QC......")
        
        out_nsQC <- scNucleosomeQC(frags = fragment, n = param(.Object)[["n"]])
        
        print("Saving results......")
        saveRDS(object = out_nsQC,
                file = output(.Object)[["nucleosomeQC.rds"]])
        
        
        # TSS QC
        
        
        
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "SCQC",
    definition = function(.Object, ...){
        report(.Object)[["nucleosomeQC.rds"]] <- output(.Object)[["nucleosomeQC.rds"]]
        .Object
    }
)




#' @name SCQC
#'
#' @title Get Single Cell Pre-processing Information
#'
#' @description
#' Get scATAC-seq pre-processing results and create Fragment Object.
#'
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSCCollect}}.
#' @param n Number of lines to read from the fragment file.
#' Default: 2000.
#' @param gene.annotation gene annotation from \link{GetGRangesFromEnsDb},
#' must in UCSC style.
#' @param peak peak regions in GRanges.
#' @param blacklist blacklist regions in GRanges.
#' @param ... Additional arguments, currently unused.
#'
#' @return An invisible \code{\link{ATACProc-class}} object scalar.
#'
#' @author Wei Zhang
#'
#' @examples
#' print(123)
setGeneric("atacSCQC", function(atacProc, ...) standardGeneric("atacSCQC"))

123131
#' @rdname SCQC
#' @aliases atacSCQC
#' @export
setMethod(
    f = "atacSCQC",
    signature = "ATACProc",
    definition = function(atacProc, n = 2000, gene.annotation = NULL, 
                          peak = NULL, blacklist = NULL, ...){
        allpara <- c(list(Class = "SCQC", prevSteps = list(atacProc)),
                     as.list(environment()), list(...))
        print(allpara)
        
        step <- do.call(new, allpara)
        invisible(step)
    }
)
