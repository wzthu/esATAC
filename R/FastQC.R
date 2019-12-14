setClass(Class = "FastQC",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "FastQC",
    definition = function(.Object,prevSteps = list(), ...){
        atacProc <- NULL
        if(length(prevSteps)>0){
            atacProc <- prevSteps[[1]]
        }
        allparam <- list(...)
        input_file <- allparam[["input_file"]]
        output_file <- allparam[["output_file"]]
        
        if((!is.null(atacProc)) && (class(atacProc)[1] == "UnzipAndMerge")){ # atacproc from UnzipAndMerge
            if(is.null(getParam(atacProc,"fastqOutput2"))){ # single end
                input(.Object)[["Input"]] <- c(as.vector(unlist(getParam(atacProc, "fastqOutput1"))))
            }else{ # paired end
                input(.Object)[["Input"]] <- c(as.vector(unlist(getParam(atacProc, "fastqOutput1"))),
                                                  as.vector(unlist(getParam(atacProc, "fastqOutput2"))))
            }
        }else if((!is.null(atacProc)) && (class(atacProc)[1] == "Renamer")){ # atacproc from renamer
            if(is.null(getParam(atacProc,"fastqOutput2"))){ # single end
                input(.Object)[["Input"]] <- c(as.vector(unlist(getParam(atacProc, "fastqOutput1"))))
            }else{ # paired end
                input(.Object)[["Input"]] <- c(as.vector(unlist(getParam(atacProc, "fastqOutput1"))),
                                                  as.vector(unlist(getParam(atacProc, "fastqOutput2"))))
            }
        }else if(((!is.null(atacProc)) && (class(atacProc)[1] != "UnzipAndMerge")) ||
                 ((!is.null(atacProc)) && (class(atacProc)[1] != "Renamer"))){
            stop("Input class must be got from 'UnzipAndMerge' or 'Renamer'!")
        }else{
            input(.Object)[["Input"]] <- input_file
        }

        if(is.null(output_file)){
            output(.Object)[["Output"]] <- getAutoPath(.Object,input(.Object)[["Input"]][1],"fastq|fq","FastQC.pdf")
        }else{
            output(.Object)[["Output"]] <- addFileSuffix(output_file,".pdf")
        }

        .Object
    }
)


setMethod(
    f = "processing",
    signature = "FastQC",
    definition = function(.Object,...){
        qQCReport(input = input(.Object)[["Input"]], pdfFilename = output(.Object)[["Output"]])

        .Object
    }
)

setMethod(
    f = "genReport",
    signature = "FastQC",
    definition = function(.Object, ...){
        report(.Object)$pdf <- output(.Object)[["Output"]]
        .Object
    }
)



#' @name FastQC
#' @title Quality control for ATAC-seq data.
#' @description
#' Generate quality control plots from fastq of ATAC-seq data.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacUnzipAndMerge}},
#' \code{\link{atacRenamer}}
#' @param input_file \code{Character} scalar.
#' Input file path. One or more(\code{vector}) fastq file path.
#' @param output_file \code{Character} scalar.
#' output file path. Defult:"input_file_QC.pdf" in the same
#' folder as your input file.
#' @param ... Additional arguments, currently unused.
#' @details Every highthroughput sequencing need quality control analysis, this
#' function provide QC for ATAC-seq, such as GC content.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream
#' analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' fra_path <- system.file("extdata", "chr20_1.2.fq.bz2", package="esATAC")
#' fq1 <- as.vector(bunzip2(filename = fra_path,
#' destname = file.path(getwd(), "chr20_1.fq"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' fra_path <- system.file("extdata", "chr20_2.2.fq.bz2", package="esATAC")
#' fq2 <- as.vector(bunzip2(filename = fra_path,
#' destname = file.path(getwd(), "chr20_2.fq"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' \dontrun{
#' qcreport(input_file = c(fq1, fq2))
#' }
#'
#'
#' @seealso
#' \code{\link{atacUnzipAndMerge}},
#' \code{\link{atacRenamer}}


setGeneric("atacQCReport",function(atacProc,
                                   input_file = NULL,
                                   output_file = NULL, ...) standardGeneric("atacQCReport"))


#' @rdname FastQC
#' @aliases atacQCReport
#' @export
setMethod(
    f = "atacQCReport",
    signature = "ATACProc",
    definition = function(atacProc,
                          input_file = NULL,
                          output_file = NULL, ...){
        allpara <- c(list(Class = "FastQC", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)


#' @rdname FastQC
#' @aliases qcreport
#' @export
qcreport <- function(input_file, output_file = NULL, ...){
    allpara <- c(list(Class = "FastQC", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
