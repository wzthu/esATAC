setClass(Class = "SCCollect",
         contains = "ATACProc"
)


#' @importFrom tools file_path_as_absolute
setMethod(
    f = "init",
    signature = "SCCollect",
    definition = function(.Object, prevSteps = list(), ...){
        allparam <- list(...)
        fragmentFile <- allparam[["fragment"]]
        csvFile <- allparam[["csv"]]
        
        atacProc <- NULL
        if(length(prevSteps) > 0){
            atacProc <- prevSteps[[1]]
        }
        
        # necessary parameters
        if(!is.null(atacProc)){
            input(.Object)[["fragmentFile"]] <- getParam(atacProc, "fragOutput")
            input(.Object)[["csvFile"]] <- getParam(atacProc, "csvOutput")
        }else{
            input(.Object)[["fragmentFile"]] <- file_path_as_absolute(fragmentFile)
            input(.Object)[["csvFile"]] <- file_path_as_absolute(csvFile)
        }
        
        if (is.null(input(.Object)[["fragmentFile"]])) {
            stop("Fragment file is missing! Please check parameter 'atacProc' or 'fragment'!")
        }
        
        if (!file.exists(input(.Object)[["csvFile"]])) {
            stop("Single Cell CSV file index is missing!")
        }
        
        # init output
        output(.Object)[["fragOB.rds"]] <- getStepWorkDir(.Object, filename = "fragmentObjects.rds")
        
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "SCCollect",
    definition = function(.Object,...){
        
        mess <- paste0("Now, reading ",
                       input(.Object)[["csvFile"]],
                       "......")
        print(mess)
        
        metadata <- read.csv(file = input(.Object)[["csvFile"]],
                             header = TRUE,
                             row.names = 1)
        
        print("Processing cell barcode.......")
        cells <- rownames(metadata)
        
        mess <- paste0("Creating Fragment Object from ",
                       input(.Object)[["fragmentFile"]],
                       "......")
        print(mess)
        fragments <- CreateFragmentObject(input(.Object)[["fragmentFile"]],
                                          cells = cells)
        
        saveRDS(object = fragments,
                file = output(.Object)[["fragOB.rds"]])
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "SCCollect",
    definition = function(.Object, ...){
        report(.Object)[["fragOB.rds"]] <- output(.Object)[["fragOB.rds"]]
        .Object
    }
)


#' @name SCCollect
#'
#' @title Get Single Cell Pre-processing Information
#'
#' @description
#' Get scATAC-seq pre-processing results and create Fragment Object.
#'
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSCBam2Frag}}.
#' @param fragment A fragment file path.
#' @param csv CSV file contains the cell information for each barcode.
#' @param ... Additional arguments, currently unused.
#'
#' @return An invisible \code{\link{ATACProc-class}} object scalar.
#'
#' @author Wei Zhang
#'
#' @examples
#' print(123)
setGeneric("atacSCCollect", function(atacProc, fragment = NULL,
                                     csv = NULL, ...) standardGeneric("atacSCCollect"))

#' @rdname SCCollect
#' @aliases atacSCCollect
#' @export
setMethod(
    f = "atacSCCollect",
    signature = "ATACProc",
    definition = function(atacProc, fragment = NULL,
                          csv = NULL, ...){
        allpara <- c(list(Class = "SCCollect", prevSteps = list(atacProc)),
                     as.list(environment()), list(...))
        step <- do.call(new, allpara)
        invisible(step)
    }
)


#' @rdname SCCollect
#' @aliases scCollect
#' @export
scCollect <- function(fragment = NULL, csv = NULL, ...){
    
    allpara <- c(list(Class = "SCCollect", prevSteps = list()),
                 as.list(environment()),
                 list(...))
    
    step <- do.call(new, allpara)
    
    invisible(step)
}
