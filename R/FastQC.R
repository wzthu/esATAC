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
#' @seealso
#' \code{\link{atacUnzipAndMerge}},
#' \code{\link{atacRenamer}},
#' \code{\link{atacSam2Bam}}.
#' @importFrom QuasR qQCReport 

#' @rdname atacQCReport
#' @export
atacQCReport <- function(atacProc = NULL, input_file = NULL, output_file = NULL){
    tmp <- FastQC$new(atacProc, input_file, output_file)
    tmp$process()
    invisible(tmp)
}

#' @rdname atacQCReport
#' @export
qcreport <- function(input_file, output_file = NULL){
    tmp <- FastQC$new(atacProc = NULL, input_file, output_file)
    tmp$process()
    invisible(tmp)
}
