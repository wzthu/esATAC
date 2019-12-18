setClass(Class = "SamToBed",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "SamToBed",
    definition = function(.Object, prevSteps = list(), ...){
        allparam <- list(...)
        samInput <- allparam[["samInput"]]
        reportOutput <- allparam[["reportOutput"]]
        merge <- allparam[["merge"]]
        posOffset <- allparam[["posOffset"]]
        negOffset <- allparam[["negOffset"]]
        chrFilterList <- allparam[["chrFilterList"]]
        bedOutput <- allparam[["bedOutput"]]
        sortBed <- allparam[["sortBed"]]
        minFragLen <- allparam[["minFragLen"]]
        maxFragLen <- allparam[["maxFragLen"]]
        saveExtLen <- allparam[["saveExtLen"]]
        uniqueBed <- allparam[["uniqueBed"]]
        
        merge <- match.arg(merge, c("auto","yes","no"))
        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                atacProc <- prevSteps[[1]]
                atacProc <- c(unlist(atacProc),list())
                atacProc <- atacProc[[length(atacProc)]]
                input(.Object)[["samInput"]] <- output(atacProc)[["samOutput"]]
                if(merge=="auto"){
                    if(property(.Object)$singleEnd){
                        merge=FALSE
                    }else{
                        merge=TRUE
                    }
                }else if(merge=="yes"){
                    merge=TRUE
                }else if(merge=="no"){
                    merge=FALSE
                }
            }

        }else{
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


        }

        if(!is.null(samInput)){
            input(.Object)[["samInput"]] <- samInput;
        }
        if(is.null(bedOutput)){
            if(!is.null(input(.Object)[["samInput"]])){
                output(.Object)[["bedOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]], "SAM|sam|Sam","bed")
            }
        }else{
            output(.Object)[["bedOutput"]] <- bedOutput;
        }
        if(is.null(reportOutput)){
            if(!is.null(input(.Object)[["samInput"]])){
                output(.Object)[["reportOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]], "SAM|sam|Sam",".report.txt")
            }
        }else{
            output(.Object)[["reportOutput"]] <- reportOutput;
        }


        param(.Object)[["merge"]] <- merge;
        param(.Object)[["posOffset"]] <- posOffset;
        param(.Object)[["negOffset"]] <- negOffset;
        param(.Object)[["filterList"]] <- chrFilterList;
        param(.Object)[["sortBed"]] <- sortBed
        param(.Object)[["uniqueBed"]] <- uniqueBed
        param(.Object)[["minFragLen"]] <- minFragLen
        param(.Object)[["maxFragLen"]] <- maxFragLen
        param(.Object)[["saveExtLen"]] <- saveExtLen

        .Object
    }
)

setMethod(
    f = "processing",
    signature = "SamToBed",
    definition = function(.Object,...){
        if(param(.Object)[["merge"]]){
            qcval<-.sam2bed_merge_call(samfile = input(.Object)[["samInput"]], bedfile = output(.Object)[["bedOutput"]],
                                       posOffset = param(.Object)[["posOffset"]], negOffset = param(.Object)[["negOffset"]],
                                       sortBed = param(.Object)[["sortBed"]],uniqueBed = param(.Object)[["uniqueBed"]],
                                       filterList = param(.Object)[["filterList"]],minFragLen = param(.Object)[["minFragLen"]],
                                       maxFragLen = param(.Object)[["maxFragLen"]],saveExtLen = param(.Object)[["saveExtLen"]] )
        }else{
            qcval<-.sam2bed_call(samfile = input(.Object)[["samInput"]], bedfile = output(.Object)[["bedOutput"]],
                                 posOffset = param(.Object)[["posOffset"]], negOffset = param(.Object)[["negOffset"]],
                                 sortBed = param(.Object)[["sortBed"]],uniqueBed = param(.Object)[["uniqueBed"]],
                                 filterList = param(.Object)[["filterList"]])
        }

        write.table(as.data.frame(qcval),file = output(.Object)[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "SamToBed",
    definition = function(.Object, ...){
        qcval <- as.list(read.table(file= output(.Object)[["reportOutput"]],header=TRUE))
        report(.Object)$table <- data.frame(
            Item = c("Total mapped reads",
                     sprintf("Chromasome %s filted reads",paste(param(.Object)[["filterList"]],collapse = "/")),
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
        report(.Object)$report <- data.frame(Item=names(qcval),Value=as.character(qcval))
        report(.Object)[["non-mitochondrial"]] <- as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]))
        report(.Object)[["non-mitochondrial-multimap"]] <- as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]) -as.integer(qcval[["multimap"]]))
        for(n in names(qcval)){
            report(.Object)[[n]] <- qcval[[n]]
        }  
        .Object
    }
)


#' @name SamToBed
#' @title Convert SAM file to BED file
#' @description
#' This function is used to convert SAM file to BED file and
#' merge interleave paired end reads,
#' shift reads,
#' filter reads according to chromosome,
#' filter reads according to fragment size,
#' sort,
#' remove duplicates reads before generating BED file.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
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
#' @param minFragLen \code{Integer} scalar
#' The minimum fragment size will be retained.
#' @param maxFragLen \code{Integer} scalar
#' The maximum fragment size will be retained.
#' @param saveExtLen \code{Logical} scaler
#' Save the fragment that are not in the range of minFragLen and maxFragLen
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{samToBed} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacBowtie2Mapping}}
#' \code{\link{bowtie2Mapping}}
#' \code{\link{atacFragLenDistr}}
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
#' setTmpDir(td)
#'
#' sambzfile <- system.file(package="esATAC", "extdata", "Example.sam.bz2")
#' samfile <- file.path(td,"Example.sam")
#' bunzip2(sambzfile,destname=samfile,overwrite=TRUE,remove=FALSE)
#' samToBed(samInput = samfile)

setGeneric("atacSamToBed",function(atacProc, reportOutput =NULL,merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= "chrM",#chrUn.*|chrM|.*random.*
                                  samInput = NULL, bedOutput = NULL, sortBed = TRUE, minFragLen = 0,maxFragLen = 100,
                                  saveExtLen = FALSE,uniqueBed = TRUE, ...) standardGeneric("atacSamToBed"))

#' @rdname SamToBed
#' @aliases atacSamToBed
#' @export
setMethod(
    f = "atacSamToBed",
    signature = "ATACProc",
    definition = function(atacProc, reportOutput =NULL,merge = c("auto","yes","no"), posOffset = +4, negOffset= -5, chrFilterList= "chrM",#chrUn.*|chrM|.*random.*
                          samInput = NULL, bedOutput = NULL, sortBed = TRUE, minFragLen = 0,maxFragLen = 100,
                          saveExtLen = FALSE,uniqueBed = TRUE, ...){
        allpara <- c(list(Class = "SamToBed", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname SamToBed
#' @aliases samToBed
#' @export
samToBed <- function(samInput, reportOutput =NULL,
                     merge = c("auto","yes","no"), 
                     posOffset = +4, negOffset= -5, 
                     chrFilterList= "chrM",#chrUn.*|chrM|.*random.*
                     bedOutput = NULL, sortBed = TRUE, 
                     minFragLen = 0,maxFragLen = 100,
                     saveExtLen = FALSE,uniqueBed = TRUE, ...){
    allpara <- c(list(Class = "SamToBed", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
