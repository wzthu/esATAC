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
#' \code{\link{atacSamToBed}},
#' \code{\link{atacBamToBed}}.
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
#' The chromatin(or regex of chromatin) will be retain/discard
#' if \code{select} is TRUE/FALSE
#' @param select \code{Logical} scalar
#' The chromatin in \code{chrFilterList} will be retain if TRUE. default: FALSE
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
#' @examples
#'
#' sam_bz <- system.file("extdata", "Example.sam.bz2", package="ATACFlow")
#' sam_path <- as.vector(bunzip2(filename = sam_bz,
#' destname = file.path(getwd(), "Example.sam"),
#' ext="bz2", FUN=bzfile, remove = FALSE))
#' samToBed(samInput = sam_path)
#'
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{atacBamToBed}}
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacBamToBed}}
#' \code{\link{atacBamToBed}}
#' \code{\link{atacBamToBed}}
#' \code{\link{atacBamToBed}}
#' \code{\link{atacBamToBed}}

#' @rdname atacSamToBed
#' @export
atacSamToBed <- function(atacProc, reportOutput =NULL,merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= "chrM",#chrUn.*|chrM|.*random.*
                        samInput = NULL, bedOutput = NULL, sortBed = TRUE, minFregLen = 0,maxFregLen = 100,
                        saveExtLen = FALSE,uniqueBed = TRUE){
    atacproc <- SamToBed$new(atacProc=atacProc, merge=merge, posOffset=posOffset, negOffset=negOffset, chrFilterList=chrFilterList,
                             samInput=samInput, bedOutput=bedOutput, sortBed=sortBed, uniqueBed=uniqueBed, minFregLen=minFregLen,
                             maxFregLen=maxFregLen, saveExtLen=saveExtLen,reportOutput=reportOutput)
    atacproc$process()
    invisible(atacproc)
}
#' @rdname atacSamToBed
#' @export
samToBed <- function(samInput, reportOutput =NULL,merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= "chrM",#chrUn.*|chrM|.*random.*
                         bedOutput = NULL, sortBed = TRUE, minFregLen = 0,maxFregLen = 100,
                        saveExtLen = FALSE,uniqueBed = TRUE){
    atacproc <- SamToBed$new(atacProc=NULL, merge=merge, posOffset=posOffset, negOffset=negOffset, chrFilterList=chrFilterList,
                             samInput=samInput, bedOutput=bedOutput, sortBed=sortBed, uniqueBed=uniqueBed, minFregLen=minFregLen,
                             maxFregLen=maxFregLen, saveExtLen=saveExtLen,reportOutput=reportOutput)
    atacproc$process()
    invisible(atacproc)
}
