setClass(Class = "SamToBed",
         contains = "ATACProc"
)

setMethod(
    f = "initialize",
    signature = "SamToBed",
    definition = function(.Object,atacProc, ..., merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= NULL,
                          reportOutput =NULL,samInput = NULL, bedOutput = NULL, sortBed = TRUE, uniqueBed = TRUE,
                          minFregLen = 0,maxFregLen = 100, saveExtLen = FALSE, editable=FALSE){
        .Object <- init(.Object,"SamToBed",editable,list(arg1=atacProc))
        merge=merge[1]
        if(!is.null(atacProc)){
            .Object@paramlist[["samInput"]] <- getParam(atacProc, "samOutput");
            regexProcName<-sprintf("(SAM|Sam|sam|%s)",getProcName(atacProc))
            if(merge=="auto"){
                if(.Object@singleEnd){
                    merge=FALSE
                }else{
                    merge=TRUE
                }

            }else if(merge=="yes"){
                merge=TRUE
            }else if(merge=="no"){
                merge=FALSE
            }else{
                stop(paste0("Invalid value of merge: ",merge))
            }

        }else{
            regexProcName<-"(SAM|Sam|sam)"
            if(!editable){
                if(merge=="auto"){
                    if(is.null(samInput)||!file.exists(samInput)){
                        stop(paste0("samInput file does not exist! ",samInput))
                    }
                    asamfile <- file(samInput, "r")
                    lines<-readLines(asamfile,n=1000)
                    close(asamfile)
                    lines<-lines[!grepl("^@",lines)]
                    merge=FALSE
                    for(i in 2:length(lines)){
                        code1=strsplit(lines[i-1],"\t")[[1]][1]
                        code2=strsplit(lines[i],"\t")[[1]][1]
                        if(code1==code2){
                            merge=TRUE
                        }
                    }
                }else if(merge=="yes"){
                    merge=TRUE
                }else if(merge=="no"){
                    merge=FALSE
                }else{
                    stop(paste0("Invalid value of merge: ",merge))
                }
            }else{
                merge=TRUE
            }


        }

        if(!is.null(samInput)){
            .Object@paramlist[["samInput"]] <- samInput;
        }
        if(is.null(bedOutput)){
            if(!is.null(.Object@paramlist[["samInput"]])){
                prefix <- getBasenamePrefix(.Object, .Object@paramlist[["samInput"]],regexProcName)
                .Object@paramlist[["bedOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".bed"))
            }
        }else{
            .Object@paramlist[["bedOutput"]] <- bedOutput;
        }
        if(is.null(reportOutput)){
            if(!is.null(.Object@paramlist[["samInput"]])){
                prefix<-getBasenamePrefix(.Object, .Object@paramlist[["samInput"]],regexProcName)
                .Object@paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".report"))
            }
        }else{
            .Object@paramlist[["reportOutput"]] <- reportOutput;
        }


        .Object@paramlist[["merge"]] <- merge;
        .Object@paramlist[["posOffset"]] <- posOffset;
        .Object@paramlist[["negOffset"]] <- negOffset;
        .Object@paramlist[["filterList"]] <- chrFilterList;
        .Object@paramlist[["sortBed"]] <- sortBed
        .Object@paramlist[["uniqueBed"]] <- uniqueBed
        .Object@paramlist[["minFregLen"]] <- minFregLen
        .Object@paramlist[["maxFregLen"]] <- maxFregLen
        .Object@paramlist[["saveExtLen"]] <- saveExtLen

        paramValidation(.Object)
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "SamToBed",
    definition = function(.Object,...){
        if(.Object@paramlist[["merge"]]){
            qcval<-.sam2bed_merge_call(samfile = .Object@paramlist[["samInput"]], bedfile = .Object@paramlist[["bedOutput"]],
                                       posOffset = .Object@paramlist[["posOffset"]], negOffset = .Object@paramlist[["negOffset"]],
                                       sortBed = .Object@paramlist[["sortBed"]],uniqueBed = .Object@paramlist[["uniqueBed"]],
                                       filterList = .Object@paramlist[["filterList"]],minFregLen = .Object@paramlist[["minFregLen"]],
                                       maxFregLen = .Object@paramlist[["maxFregLen"]],saveExtLen = .Object@paramlist[["saveExtLen"]] )
        }else{
            qcval<-.sam2bed_call(samfile = .Object@paramlist[["samInput"]], bedfile = .Object@paramlist[["bedOutput"]],
                                 posOffset = .Object@paramlist[["posOffset"]], negOffset = .Object@paramlist[["negOffset"]],
                                 sortBed = .Object@paramlist[["sortBed"]],uniqueBed = .Object@paramlist[["uniqueBed"]],
                                 filterList = .Object@paramlist[["filterList"]])
        }

        write.table(as.data.frame(qcval),file = .Object@paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "SamToBed",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["samInput"]])){
            stop(paste("samInput is requied"));
        }
    }
)



setMethod(
    f = "checkAllPath",
    signature = "SamToBed",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["samInput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["bedOutput"]]);
    }
)


setMethod(
    f = "getReportValImp",
    signature = "SamToBed",
    definition = function(.Object, item){
        qcval <- as.list(read.table(file= .Object@paramlist[["reportOutput"]],header=TRUE))
        if(item == "report"){
            data.frame(
                Item = c("Total mapped reads",
                         sprintf("Chromasome %s filted reads",paste(.Object@paramlist[["filterList"]],collapse = "/")),
                         "Filted multimap reads",
                         "Removed fragment size out of range",
                         "Removed duplicate reads"
                ),
                Retain = c(qcval[["total"]],
                           as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]])),
                           as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]) -as.integer(qcval[["multimap"]])),
                           as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]) -as.integer(qcval[["multimap"]] - as.integer(qcval[["extlen"]]))),
                           qcval[["save"]]

                ),
                Filted = c("/",
                           qcval[["filted"]],
                           qcval[["multimap"]],
                           qcval[["unique"]],
                           qcval[["extlen"]]
                )

            )
            return(data.frame(Item=names(qcval),Value=as.character(qcval)))
        }else if(item == "non-mitochondrial")(
            return(as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]])))
        )else if(item == "non-mitochondrial-multimap"){
            return(as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]) -as.integer(qcval[["multimap"]])))
        }else{
            return(qcval[[item]])
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "SamToBed",
    definition = function(.Object){
        return(c("report","total","save","filted","extlen","unique","multimap","non-mitochondrial","non-mitochondrial-multimap"))
    }
)




SamToBed <- R6::R6Class(
  classname = "SamToBed",
  inherit = ATACProc,
  public = list(
    initialize = function(atacProc, merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= NULL,
                          reportOutput =NULL,samInput = NULL, bedOutput = NULL, sortBed = TRUE, uniqueBed = TRUE,
                          minFregLen = 0,maxFregLen = 100, saveExtLen = FALSE, editable=FALSE){
      super$initialize("SamToBed",editable,list(arg1=atacProc))
      merge=merge[1]
      if(!is.null(atacProc)){
          private$paramlist[["samInput"]] <- atacProc$getParam("samOutput");
          regexProcName<-sprintf("(SAM|Sam|sam|%s)",atacProc$getProcName())
          if(merge=="auto"){
              if(private$singleEnd){
                  merge=FALSE
              }else{
                  merge=TRUE
              }

          }else if(merge=="yes"){
              merge=TRUE
          }else if(merge=="no"){
              merge=FALSE
          }else{
              stop(paste0("Invalid value of merge: ",merge))
          }

      }else{
          regexProcName<-"(SAM|Sam|sam)"
          if(!editable){
              if(merge=="auto"){
                  if(is.null(samInput)||!file.exists(samInput)){
                      stop(paste0("samInput file does not exist! ",samInput))
                  }
                  asamfile <- file(samInput, "r")
                  lines<-readLines(asamfile,n=1000)
                  close(asamfile)
                  lines<-lines[!grepl("^@",lines)]
                  merge=FALSE
                  for(i in 2:length(lines)){
                      code1=strsplit(lines[i-1],"\t")[[1]][1]
                      code2=strsplit(lines[i],"\t")[[1]][1]
                      if(code1==code2){
                          merge=TRUE
                      }
                  }
              }else if(merge=="yes"){
                    merge=TRUE
              }else if(merge=="no"){
                    merge=FALSE
              }else{
                  stop(paste0("Invalid value of merge: ",merge))
              }
          }else{
              merge=TRUE
          }


      }

      if(!is.null(samInput)){
          private$paramlist[["samInput"]] <- samInput;
      }
      if(is.null(bedOutput)){
          if(!is.null(private$paramlist[["samInput"]])){
              prefix<-private$getBasenamePrefix(private$paramlist[["samInput"]],regexProcName)
              private$paramlist[["bedOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".bed"))
          }
      }else{
          private$paramlist[["bedOutput"]] <- bedOutput;
      }
      if(is.null(reportOutput)){
          if(!is.null(private$paramlist[["samInput"]])){
              prefix<-private$getBasenamePrefix(private$paramlist[["samInput"]],regexProcName)
              private$paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
          }
      }else{
          private$paramlist[["reportOutput"]] <- reportOutput;
      }


      private$paramlist[["merge"]] <- merge;
      private$paramlist[["posOffset"]] <- posOffset;
      private$paramlist[["negOffset"]] <- negOffset;
      private$paramlist[["filterList"]] <- chrFilterList;
      private$paramlist[["sortBed"]] <- sortBed
      private$paramlist[["uniqueBed"]] <- uniqueBed
      private$paramlist[["minFregLen"]] <- minFregLen
      private$paramlist[["maxFregLen"]] <- maxFregLen
      private$paramlist[["saveExtLen"]] <- saveExtLen


      private$paramValidation()
    }
  ),

    private = list(
        processing = function(){
            if(private$paramlist[["merge"]]){
                qcval<-.sam2bed_merge_call(samfile = private$paramlist[["samInput"]], bedfile = private$paramlist[["bedOutput"]],
                                                                       posOffset = private$paramlist[["posOffset"]], negOffset = private$paramlist[["negOffset"]],
                                                                       sortBed = private$paramlist[["sortBed"]],uniqueBed = private$paramlist[["uniqueBed"]],
                                                                       filterList = private$paramlist[["filterList"]],minFregLen = private$paramlist[["minFregLen"]],
                                                                       maxFregLen = private$paramlist[["maxFregLen"]],saveExtLen = private$paramlist[["saveExtLen"]] )
            }else{
                qcval<-.sam2bed_call(samfile = private$paramlist[["samInput"]], bedfile = private$paramlist[["bedOutput"]],
                                                                 posOffset = private$paramlist[["posOffset"]], negOffset = private$paramlist[["negOffset"]],
                                                                 sortBed = private$paramlist[["sortBed"]],uniqueBed = private$paramlist[["uniqueBed"]],
                                                                 filterList = private$paramlist[["filterList"]])
            }

            write.table(as.data.frame(qcval),file = private$paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["samInput"]])){
                stop(paste("samInput is requied"));
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["samInput"]]);
            private$checkFileCreatable(private$paramlist[["bedOutput"]]);
        },
        getReportValImp = function(item){
            qcval <- as.list(read.table(file= private$paramlist[["reportOutput"]],header=TRUE))
            if(item == "report"){
                data.frame(
                    Item = c("Total mapped reads",
                        sprintf("Chromasome %s filted reads",paste(private$paramlist[["filterList"]],collapse = "/")),
                        "Filted multimap reads",
                        "Removed fragment size out of range",
                        "Removed duplicate reads"
                    ),
                    Retain = c(qcval[["total"]],
                              as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]])),
                              as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]) -as.integer(qcval[["multimap"]])),
                              as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]) -as.integer(qcval[["multimap"]] - as.integer(qcval[["extlen"]]))),
                              qcval[["save"]]

                    ),
                    Filted = c("/",
                               qcval[["filted"]],
                               qcval[["multimap"]],
                               qcval[["unique"]],
                               qcval[["extlen"]]
                    )

                )
                return(data.frame(Item=names(qcval),Value=as.character(qcval)))
            }else if(item == "non-mitochondrial")(
                return(as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]])))
            )else if(item == "non-mitochondrial-multimap"){
                return(as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]) -as.integer(qcval[["multimap"]])))
            }else{
                return(qcval[[item]])
            }
        },
        getReportItemsImp = function(){
            return(c("report","total","save","filted","extlen","unique","multimap","non-mitochondrial","non-mitochondrial-multimap"))
        }
  )


)





#' @name atacSamToBed
#' @aliases atacSamToBed
#' @aliases samToBed
#' @title Convert SAM file to BED file
#' @description
#' This function is used to convert SAM file to BED file and
#' merge interleave paired end reads,
#' shift reads,
#' filter reads according to chromosome,
#' filter reads according to fregment size,
#' sort,
#' remove duplicates reads before generating BED file.
#' @param atacProc \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacBowtie2Mapping}}
#' \code{\link{bowtie2Mapping}}
#' @param samInput \code{Character} scalar.
#' SAM file input path.
#' @param bedOutput \code{Character} scalar.
#' Bed file output path.
#' @param reportOutput \code{Character} scalar
#' report file path
#' @param merge \code{Logical} scalar
#' Merge paired end reads.
#' @param posOffset \code{Integer} scalar
#' The offset that positive strand reads will shift.
#' @param negOffset \code{Integer} scalar
#' The offset that negative strand reads will shift.
#' @param chrFilterList \code{Character} vector
#' The chromatin(or regex of chromatin) will be discard
#' @param sortBed \code{Logical} scalar
#' Sort bed file in the order of chromatin, start, end
#' @param uniqueBed \code{Logical} scalar
#' Remove duplicates reads in bed if TRUE. default: FALSE
#' @param minFregLen \code{Integer} scalar
#' The minimum fregment size will be retained.
#' @param maxFregLen \code{Integer} scalar
#' The maximum fregment size will be retained.
#' @param saveExtLen \code{Logical} scaler
#' Save the fregment that are not in the range of minFregLen and MaxFregLen
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{samToBed} instead.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacBowtie2Mapping}}
#' \code{\link{bowtie2Mapping}}
#' \code{\link{atacFregLenDistr}}
#' \code{\link{atacExtractCutSite}}
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacBedUtils}}
#' \code{\link{atacTSSQC}}
#' \code{\link{atacBedToBigWig}}
#'
#' @examples
#' library(R.utils)
#' library(magrittr)
#' td <- tempdir()
#' setConfigure("tmpdir",td)
#'
#' sambzfile <- system.file(package="ATACpipe", "extdata", "Example.sam.bz2")
#' samfile <- file.path(td,"Example.sam")
#' bunzip2(sambzfile,destname=samfile,overwrite=TRUE,remove=FALSE)
#' samToBed(samInput = samfile)
#'
#' @rdname atacSamToBed
#' @exportMethod atacSamToBed
setGeneric("atacSamToBed",function(atacProc, reportOutput =NULL,merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= "chrM",#chrUn.*|chrM|.*random.*
                                  samInput = NULL, bedOutput = NULL, sortBed = TRUE, minFregLen = 0,maxFregLen = 100,
                                  saveExtLen = FALSE,uniqueBed = TRUE) standardGeneric("atacSamToBed"))
setMethod(
    f = "atacSamToBed",
    signature = "ATACProc",
    definition = function(atacProc, reportOutput =NULL,merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= "chrM",#chrUn.*|chrM|.*random.*
                          samInput = NULL, bedOutput = NULL, sortBed = TRUE, minFregLen = 0,maxFregLen = 100,
                          saveExtLen = FALSE,uniqueBed = TRUE){
        atacproc <- new(
            "SamToBed",
            atacProc = atacProc,
            reportOutput = reportOutput,
            fastqOutput1 = fastqOutput1,
            merge = merge,
            posOffset = posOffset,
            negOffset = negOffset,
            chrFilterList = chrFilterList,
            samInput = samInput,
            bedOutput = bedOutput,
            sortBed = sortBed,
            minFregLen = minFregLen,
            maxFregLen = maxFregLen,
            saveExtLen = saveExtLen,
            uniqueBed = uniqueBed)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)
#' @rdname atacSamToBed
#' @export
samToBed <- function(samInput, reportOutput =NULL,merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= "chrM",#chrUn.*|chrM|.*random.*
                         bedOutput = NULL, sortBed = TRUE, minFregLen = 0,maxFregLen = 100,
                        saveExtLen = FALSE,uniqueBed = TRUE){
    atacproc <- new(
        "SamToBed",
        atacProc = NULL,
        reportOutput = reportOutput,
        fastqOutput1 = fastqOutput1,
        merge = merge,
        posOffset = posOffset,
        negOffset = negOffset,
        chrFilterList = chrFilterList,
        samInput = samInput,
        bedOutput = bedOutput,
        sortBed = sortBed,
        minFregLen = minFregLen,
        maxFregLen = maxFregLen,
        saveExtLen = saveExtLen,
        uniqueBed = uniqueBed)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
