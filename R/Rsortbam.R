Rsortbam <- R6::R6Class(
    classname = "Rsortbam",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, bamInput = NULL, bamOutput = NULL, editable=FALSE){
            super$initialize("Rsortbam", editable, list(arg1 = atacProc))

            # necessary parameters
            if(!is.null(atacProc)){ # atacproc from SamToBam
                private$paramlist[["bamInput"]] <- atacProc$getParam("bamOutput")
                regexProcName <- sprintf("(bam|%s)", atacProc$getProcName())
            }else if(is.null(atacProc)){ # input
                private$paramlist[["bamInput"]] <- bamInput
                regexProcName <- "(bam)"
            }
            # unnecessary parameters
            if(is.null(bamOutput)){
                prefix <- private$getBasenamePrefix(private$paramlist[["bamInput"]], regexProcName)
                private$paramlist[["bamOutput_tmp"]] <- file.path(.obtainConfigure("tmpdir"),
                                                                  paste0(prefix, ".", self$getProcName()))
                private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["bamOutput_tmp"]], ".bam", collapse = "")
            }else{
                name_split <- unlist(base::strsplit(x = bamOutput, split = ".", fixed = TRUE))
                suffix <- tail(name_split, 1)
                if(suffix == "bam"){
                    private$paramlist[["bamOutput_tmp"]] <- paste0(name_split[-length(name_split)], collapse = ".")
                    private$paramlist[["bamOutput"]] <- bamOutput
                }else{
                    private$paramlist[["bamOutput_tmp"]] <- bamOutput
                    private$paramlist[["bamOutput"]] <- paste0(private$paramlist[["bamOutput_tmp"]], ".bam", collapse = "")
                }
            }

            # parameter check
            private$paramValidation()
        } # initialization end

    ), # public end

    private = list(
        processing = function(){
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("source:%s", private$paramlist[["bamInput"]]))
            private$writeLog(sprintf("destination:%s", private$paramlist[["bamOutput"]]))
            Rsamtools::sortBam(file = private$paramlist[["bamInput"]], destination = private$paramlist[["bamOutput_tmp"]])
            Rsamtools::indexBam(files = private$paramlist[["bamOutput"]])
        }, # processing end

        checkRequireParam = function(){
            if(is.null(private$paramlist[["bamInput"]])){
                stop("Parameter bamInput or atacProc is required!")
            }
        }, # checkRequireParam end

        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bamInput"]])
            private$checkPathExist(private$paramlist[["bamOutput"]])
        } # checkAllPath end

    ) # private end

) # class end

#' @name atacBamSort
#' @aliases atacBamSort
#' @aliases bamsort
#' @title Sort bam file and rebuild bai index.
#' @param atacProc \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSam2Bam}}.
#' @param bamInput \code{Character} scalar.
#' Input bam file path.
#' @param bamOutput \code{Character} scalar.
#' Output bam file path.
#' @return An invisible \code{\link{ATACProc}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(Rsamtools)
#' ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' bamsort(bamInput = ex1_file)
#'
#' @seealso
#' \link[Rsamtools]{sortBam}
#' \link[Rsamtools]{indexBam}
#' \code{\link{atacSam2Bam}}
#' \code{\link{atacBam2Bed}}

#' @rdname atacBamSort
#' @export
atacBamSort <- function(atacProc = NULL, bamInput = NULL, bamOutput = NULL){
    tmp <- Rsortbam$new(atacProc, bamInput, bamOutput)
    tmp$process()
    invisible(tmp)
}

#' @rdname atacBamSort
#' @export
bamsort <- function(bamInput = NULL, bamOutput = NULL){
    tmp <- Rsortbam$new(atacProc = NULL, bamInput, bamOutput)
    tmp$process()
    invisible(tmp)
}
