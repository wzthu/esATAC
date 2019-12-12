
setClass(Class = "SingleRepReport",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "SingleRepReport",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        htmlOutput <- allparam[["htmlOutput"]]
        
        
        if(is.null(htmlOutput)){
            output(.Object)$htmlOutput <- getStepWorkDir(.Object, filename = "Report.html")
        }else{
            output(.Object)$htmlOutput <- getStepWorkDir(.Object, filename = htmlOutput)
        }
        
        .Object
    }
)

#' @importFrom rmarkdown render
#' @importFrom pipeFrame stepType

setMethod(
    f = "processing",
    signature = "SingleSampleReport",
    definition = function(.Object, ...){
        htmlOutput <- output(.Object)[["htmlOutput"]]
        prevSteps <- list(...)[["prevSteps"]]
        prevStepsType <- lapply(prevSteps, function(step){
            return(stepType(step))
        })
        names(prevSteps) <- unlist(prevStepsType)
        
        save(prevSteps, file = getStepWorkDir(.Object = .Object, filename = "PrevSteps.Rdata"))
        
        reportmkd <- getStepWorkDir(.Object = .Object, filename = "Report.Rmd")
        
        reportmkd1 <- getStepWorkDir(.Object = .Object, filename = "Report.code.Rmd")
        
        file.copy(from = system.file(package = "esATAC", "extdata","Report.Rmd"),
                  to = reportmkd,overwrite = TRUE)
        
  #      file.copy(from = system.file(package = "enrichTF", "extdata","Report.code.Rmd"),
  #                to = reportmkd1,overwrite = TRUE)
        
        render(reportmkd)
        
  #      render(reportmkd1, quiet = TRUE)
        
        
        .Object
    }
)


#' @name SingleSampleReport
#' @importFrom rtracklayer import
#' @importFrom rtracklayer import.bed
#' @title Final report for single group of regions
#' @description
#' When user call all steps in the pipeline, the final report can be generated.
#' @param prevStep \code{\link{Step-class}} object scalar.
#' Any steps object in this package is acceptable when the pipeline is ready.
#' @param htmlOutput \code{Character} scalar.
#' HTML report file directory
#' Default: NULL ("Report.html")
#' @param ... Additional arguments, currently unused.
#' @details
#' The report is HTML format. All link in HTML file is the relative directory
#' in report step folder and other step folder
#' If user want to move HTML file and keep all link access available,
#' they should move the whole pipeline folder at the same time.
#' @return An invisible \code{\link{ATACProc-class}}
#' object (\code{\link{Step-class}} based) scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacPipe}}



setGeneric("atacSingleRepReport",
           function(prevStep, htmlOutput = NULL,...)
               standardGeneric("atacSingleRepReport"))



#' @rdname SingleRepReport
#' @aliases enrichSingleSampleReport
#' @export
setMethod(
    f = "atacSingleRepReport",
    signature = "Step",
    definition = function(prevStep, htmlOutput = NULL, ...){
        allpara <- c(list(Class = "SingleRepReport",
                          prevSteps = list(prevStep), isReportStep = TRUE),
                     as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
