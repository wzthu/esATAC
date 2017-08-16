FastQC <- R6::R6Class(
  classname = "FastQC",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, input_file = NULL, output_file = NULL, editable = FALSE){
      super$initialize("FastQC", editable, list(arg1 = atacProc))

      # necessary and unchanged parameters, this should be tested
      # in this class, there is no necessary and unchanged parameters
      # input processing
      if((!is.null(atacProc)) && (class(atacProc)[1] == "UnzipAndMerge")){ # atacproc from UnzipAndMerge
        print("UnzipAndMerge")
        if(is.null(atacProc$getParam("fastqOutput2"))){ # single end
          private$paramlist[["Input1"]] <- as.vector(unlist(atacProc$getParam("fastqOutput1")))
          private$paramlist[["Input"]] <- c(private$paramlist[["Input1"]])
        }else{ # paired end
          private$paramlist[["Input1"]] <- as.vector(unlist(atacProc$getParam("fastqOutput1")))
          private$paramlist[["Input2"]] <- as.vector(unlist(atacProc$getParam("fastqOutput2")))
          private$paramlist[["Input"]] <- c(private$paramlist[["Input1"]], private$paramlist[["Input2"]])
        }
      }else if((!is.null(atacProc)) && (class(atacProc)[1] == "Renamer")){ # atacproc from renamer
        print("Renamer")
        if(is.null(atacProc$getParam("fastqOutput2"))){ # single end
          private$paramlist[["Input1"]] <- as.vector(unlist(atacProc$getParam("fastqOutput1")))
          private$paramlist[["Input"]] <- c(private$paramlist[["Input1"]])
        }else{ # paired end
          private$paramlist[["Input1"]] <- as.vector(unlist(atacProc$getParam("fastqOutput1")))
          private$paramlist[["Input2"]] <- as.vector(unlist(atacProc$getParam("fastqOutput2")))
          private$paramlist[["Input"]] <- c(private$paramlist[["Input1"]], private$paramlist[["Input2"]])
        }
      }else if((!is.null(atacProc)) && (class(atacProc)[1] == "SamToBam")){
        print("SamToBam")
        private$paramlist[["Input1"]] <- as.vector(unlist(atacProc$getParam("bamOutput")))
        private$paramlist[["Input"]] <- c(private$paramlist[["Input1"]])
      }else if(((!is.null(atacProc)) && (class(atacProc)[1] != "UnzipAndMerge")) ||
               ((!is.null(atacProc)) && (class(atacProc)[1] != "Renamer")) ||
               ((!is.null(atacProc)) && (class(atacProc)[1] != "SamToBam"))){
        stop("Input class must be got from 'UnzipAndMerge' or 'Renamer' or 'SamToBam'!")
      }else if(is.null(atacProc) && (length(input_file) == 1)){ # input one file
        private$paramlist[["Input1"]] <- as.vector(unlist(input_file))
        private$paramlist[["Input"]] <- c(private$paramlist[["Input1"]])
      }else if(is.null(atacProc) && (length(input_file) == 2)){ # input two file
        private$paramlist[["Input1"]] <- as.vector(unlist(input_file))[1]
        private$paramlist[["Input2"]] <- as.vector(unlist(input_file))[2]
        private$paramlist[["Input"]] <- c(private$paramlist[["Input1"]], private$paramlist[["Input2"]])
      }
      #output processing
      if(is.null(output_file)){
        file_path <- dirname(private$paramlist[["Input1"]])[1]
        private$paramlist[["Output"]] <- paste0(file_path, "/FastQC.pdf", collapse = "")
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
      if(is.null(private$paramlist[["Input1"]])){
        stop("Parameter input_file is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["Input1"]])
      private$checkPathExist(private$paramlist[["Output"]])
    } # checkAllPath end

  ) # private end

) # class end



#' Quality control using Quasr::qQCreport
#' @param input_file a vector.
#' @export
atacQCReport <- function(atacProc = NULL, input_file = NULL, output_file = NULL){
  tmp <- FastQC$new(atacProc, input_file, output_file)
  tmp$process()
  return(tmp)
}
