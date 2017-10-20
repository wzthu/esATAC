setClass(Class = "FastQC",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "FastQC",
    definition = function(.Object, atacProc, ..., input_file = NULL, output_file = NULL, editable = FALSE){
        .Object <- init(.Object,"FastQC",editable,list(arg1=atacProc))

        if((!is.null(atacProc)) && (class(atacProc)[1] == "UnzipAndMerge")){ # atacproc from UnzipAndMerge
            if(is.null(atacProc$getParam("fastqOutput2"))){ # single end
                .Object@paramlist[["Input"]] <- c(as.vector(unlist(getParam(atacProc, "fastqOutput1"))))
            }else{ # paired end
                .Object@paramlist[["Input"]] <- c(as.vector(unlist(getParam(atacProc, "fastqOutput1"))),
                                                  as.vector(unlist(getParam(atacProc, "fastqOutput2"))))
            }
        }else if((!is.null(atacProc)) && (class(atacProc)[1] == "Renamer")){ # atacproc from renamer
            if(is.null(atacProc$getParam("fastqOutput2"))){ # single end
                .Object@paramlist[["Input"]] <- c(as.vector(unlist(getParam(atacProc, "fastqOutput1"))))
            }else{ # paired end
                .Object@paramlist[["Input"]] <- c(as.vector(unlist(getParam(atacProc, "fastqOutput1"))),
                                                  as.vector(unlist(getParam(atacProc, "fastqOutput2"))))
            }
        }else if((!is.null(atacProc)) && (class(atacProc)[1] == "SamToBam")){
            .Object@paramlist[["Input"]] <- c(as.vector(unlist(getParam(atacProc, "bamOutput"))))
        }else if(((!is.null(atacProc)) && (class(atacProc)[1] != "UnzipAndMerge")) ||
                 ((!is.null(atacProc)) && (class(atacProc)[1] != "Renamer")) ||
                 ((!is.null(atacProc)) && (class(atacProc)[1] != "SamToBam"))){
            stop("Input class must be got from 'UnzipAndMerge' or 'Renamer' or 'SamToBam'!")
        }else{
            .Object@paramlist[["Input"]] <- input_file
        }

        if(is.null(output_file)){
            output_name <- paste(basename(tools::file_path_sans_ext(.Object@paramlist[["Input"]][1])),
                                 "_FastQC.pdf", sep = "")
            .Object@paramlist[["Output"]] <- file.path(.obtainConfigure("tmpdir"), output_name)
        }else{
            .Object@paramlist[["Output"]] <- output_file
        }

        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "FastQC",
    definition = function(.Object,...){
        .Object <- writeLog(.Object,paste0("processing file:"))
        .Object <- writeLog(.Object,sprintf("source:%s", .Object@paramlist[["Input"]]))
        .Object <- writeLog(.Object,sprintf("destination:%s", .Object@paramlist[["Output"]]))
        QuasR::qQCReport(input = .Object@paramlist[["Input"]], pdfFilename = .Object@paramlist[["Output"]])
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "FastQC",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["Input"]])){
            stop("Parameter input_file is required!")
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "FastQC",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["Input"]])
        checkPathExist(.Object,.Object@paramlist[["Output"]])
    }
)


setMethod(
    f = "getReportValImp",
    signature = "FastQC",
    definition = function(.Object, item){
        if(item == "pdf"){
            return(.Object@paramlist[["Output"]])
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "FastQC",
    definition = function(.Object){
        return(c("pdf"))
    }
)



FastQC <- R6::R6Class(
    classname = "FastQC",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, input_file = NULL, output_file = NULL, editable = FALSE){
            super$initialize("FastQC", editable, list(arg1 = atacProc))

            # necessary and unchanged parameters, this should be tested
            # in this class, there is no necessary and unchanged parameters
            # input processing
            if((!is.null(atacProc)) && (class(atacProc)[1] == "UnzipAndMerge")){ # atacproc from UnzipAndMerge
                if(is.null(atacProc$getParam("fastqOutput2"))){ # single end
                    private$paramlist[["Input"]] <- c(as.vector(unlist(atacProc$getParam("fastqOutput1"))))
                }else{ # paired end
                    private$paramlist[["Input"]] <- c(as.vector(unlist(atacProc$getParam("fastqOutput1"))),
                                                      as.vector(unlist(atacProc$getParam("fastqOutput2"))))
                }
            }else if((!is.null(atacProc)) && (class(atacProc)[1] == "Renamer")){ # atacproc from renamer
                if(is.null(atacProc$getParam("fastqOutput2"))){ # single end
                    private$paramlist[["Input"]] <- c(as.vector(unlist(atacProc$getParam("fastqOutput1"))))
                }else{ # paired end
                    private$paramlist[["Input"]] <- c(as.vector(unlist(atacProc$getParam("fastqOutput1"))),
                                                      as.vector(unlist(atacProc$getParam("fastqOutput2"))))
                }
            }else if((!is.null(atacProc)) && (class(atacProc)[1] == "SamToBam")){
                private$paramlist[["Input"]] <- c(as.vector(unlist(atacProc$getParam("bamOutput"))))
            }else if(((!is.null(atacProc)) && (class(atacProc)[1] != "UnzipAndMerge")) ||
                     ((!is.null(atacProc)) && (class(atacProc)[1] != "Renamer")) ||
                     ((!is.null(atacProc)) && (class(atacProc)[1] != "SamToBam"))){
                stop("Input class must be got from 'UnzipAndMerge' or 'Renamer' or 'SamToBam'!")
            }else{
                private$paramlist[["Input"]] <- input_file
            }
            #output processing
            if(is.null(output_file)){
                output_name <- paste(basename(tools::file_path_sans_ext(private$paramlist[["Input"]][1])),
                                     "_FastQC.pdf", sep = "")
                private$paramlist[["Output"]] <- file.path(.obtainConfigure("tmpdir"), output_name)
            }else{
                private$paramlist[["Output"]] <- output_file
            }

            # parameter check
            private$paramValidation()
        } # initialization end

    ),

    private = list(
        processing = function(){
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("source:%s", private$paramlist[["Input"]]))
            private$writeLog(sprintf("destination:%s", private$paramlist[["Output"]]))
            print(private$paramlist[["Input"]])
            QuasR::qQCReport(input = private$paramlist[["Input"]], pdfFilename = private$paramlist[["Output"]])
        }, # processing end

        checkRequireParam = function(){
            if(is.null(private$paramlist[["Input"]])){
                stop("Parameter input_file is required!")
            }
        }, # checkRequireParam end

        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["Input"]])
            private$checkPathExist(private$paramlist[["Output"]])
        }, # checkAllPath end

        getReportValImp = function(item){
            if(item == "pdf"){
                return(private$paramlist[["Output"]])
            }
        }, # getReportValImp end

        getReportItemsImp = function(){
            return(c("pdf"))
        } # getReportItemsImp end
    ) # private end

) # class end



#' @name atacQCReport
#' @aliases atacQCReport
#' @aliases qcreport
#' @title Quality control for ATAC-seq data.
#' @description
#' Generate quality control plots from fastq/fasta/bam of ATAC-seq data.
#' @param atacProc \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacUnzipAndMerge}},
#' \code{\link{atacRenamer}},
#' \code{\link{atacSam2Bam}}.
#' @param input_file \code{Character} scalar.
#' Input file path. One or two(\code{vector}) fastq/a or a bam file path.
#' @param output_file \code{Character} scalar.
#' output file path. Defult:"input_file_FastQC.pdf" in the same
#' folder as your input file.
#' @details Every highthroughput sequencing need quality control analysis, this
#' function provide QC for ATAC-seq, such as GC content.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream
#' analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' fra_path <- system.file("extdata", "chr20_1.fq.bz2", package="ATACpipe")
#' fq1 <- as.vector(bunzip2(filename = fra_path,
#' destname = file.path(getwd(), "chr20_1.fq"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' fra_path <- system.file("extdata", "chr20_2.fq.bz2", package="ATACpipe")
#' fq2 <- as.vector(bunzip2(filename = fra_path,
#' destname = file.path(getwd(), "chr20_2.fq"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' qcreport(input_file = c(fq1, fq2))
#'
#'
#' @seealso
#' \code{\link{atacUnzipAndMerge}},
#' \code{\link{atacRenamer}},
#' \code{\link{atacSam2Bam}}.
#' @importFrom QuasR qQCReport

#' @rdname atacQCReport
#' @exportMethod atacQCReport
setGeneric("atacQCReport",function(atacProc = NULL,
                                   input_file = NULL,
                                   output_file = NULL) standardGeneric("atacQCReport"))
setMethod(
    f = "atacQCReport",
    signature = "ATACProc",
    definition = function(atacProc = NULL,
                          input_file = NULL,
                          output_file = NULL){
        atacproc <- new(
            "FastQC",
            atacProc = atacProc,
            input_file = input_file,
            output_file = output_file)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)


#' @rdname atacQCReport
#' @export
qcreport <- function(input_file, output_file = NULL){
    atacproc <- new(
        "FastQC",
        atacProc = NULL,
        input_file = input_file,
        output_file = output_file)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
