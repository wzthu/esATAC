setClass(Class = "SCSamToBam",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "SCSamToBam",
    definition = function(.Object,prevSteps = list(),...){
        atacProc <- NULL
        if(length(prevSteps)>0){
            atacProc <- prevSteps[[1]]
        }
        allparam <- list(...)
        samInput <- allparam[["samInput"]]
        bamOutput <- allparam[["bamOutput"]]
        uniqueBamOutput <- allparam[["uniqueBamOutput"]]
        statCsvOutput <- allparam[["statCsvOutput"]]
        tsvOutput <-allparam[["tsvOutput"]]
        isSort <- allparam[["isSort"]]
        
        # necessary parameters
        if((!is.null(atacProc)) ){
            input(.Object)[["samInput"]] <- getParam(atacProc, "samOutput")
        }else if(is.null(atacProc)){ # input
            input(.Object)[["samInput"]] <- samInput
        }
        # unnecessary parameters
        if(is.null(bamOutput)){
            param(.Object)[["name_tmp"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","" )
        
            output(.Object)[["bamOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","bam")
            if(isSort){
                output(.Object)[["baiOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","bam.bai")
            }
        }else{
            bamOutput <- addFileSuffix(bamOutput,".bam")
            param(.Object)[["name_tmp"]] <- substring(bamOutput,first = 1,last = nchar(bamOutput) - 4)
            output(.Object)[["bamOutput"]] <- bamOutput
            if(isSort){
                output(.Object)[["baiOutput"]] <- addFileSuffix(bamOutput,".bai")
            }
            
        }
        if(is.null(uniqueBamOutput)){
            param(.Object)[["unique_name_tmp"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","unique" )

            output(.Object)[["uniqueBamOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","unique.bam")
            if(isSort){
                output(.Object)[["uniqueBaiOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","unique.bam.bai")
            }
        }else{
            uniqueBamOutput <- addFileSuffix(uniqueBamOutput,".bam")
            param(.Object)[["unique_name_tmp"]] <- substring(uniqueBamOutput,first = 1,last = nchar(uniqueBamOutput) - 4)
            output(.Object)[["uniqueBamOutput"]] <- uniqueBamOutput
            if(isSort){
                output(.Object)[["uniqueBaiOutput"]] <- addFileSuffix(bamOutput,".bai")
            }
        }
        if(is.null(statCsvOutput)){
            output(.Object)[["statCsvOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","csv")
        }
        if(is.null(tsvOutput)){
            output(.Object)[["tsvOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","stat.tsv")
        }
        param(.Object)[['isSort']] <- isSort
        .Object

    }
)




setMethod(
    f = "processing",
    signature = "SCSamToBam",
    definition = function(.Object,...){
        samInput <- input(.Object)[["samInput"]]
        bamOutput <- output(.Object)[["bamOutput"]]
        desOutput <- param(.Object)[["name_tmp"]]
        desOutputByBarcode <- paste0(desOutput,'.byBarcode')
        uniqueBamOutput <- output(.Object)[["uniqueBamOutput"]]
        uniqueDesOutput <- param(.Object)[["unique_name_tmp"]]
        statCsvOutput <- output(.Object)[["statCsvOutput"]]
        tsvOutput <- output(.Object)[["tsvOutput"]]
        isSort <- param(.Object)[["isSort"]]        
        if(file.exists(bamOutput)){
            file.remove(bamOutput)
        }
        if(file.exists(uniqueBamOutput)){
            file.remove(uniqueBamOutput)
        }
        if(file.exists(tsvOutput)){
            file.remove(tsvOutput)
        }
        if(file.exists(statCsvOutput)){
            file.remove(statCsvOutput)
        }
        Rsamtools::asBam(file = samInput, 
                         destination = desOutput,
                         overwrite = TRUE, indexDestination = isSort)
        Rsamtools::sortBam(file = bamOutput, destination = desOutputByBarcode, byQname=TRUE, maxMemory=15 * 1024)
        Rsamtools::asSam(file=paste0(desOutputByBarcode,'.bam'),overwrite=TRUE)
        
        samfile <- file(paste0(desOutputByBarcode,'.sam'),'r')
        headLine <- c()
        while(TRUE){
            lines <- readLines(samfile,n=1)
            if(length(lines)==0||!startsWith(lines,'@')){
                break
            }else{
                headLine <- c(headLine,lines)
            }
        }
        lines <- c(lines, readLines(samfile,n=10000000+1))
        barcode <- NULL
        vec_left_list <- list()
        left_lines <- c()
        endflag <- FALSE
        stat_list <- list()
        count <- 1
        nb_barcode <- 0
        while(TRUE){
            vec_list <- c(vec_left_list,strsplit(lines, '\t'))
            lines <- c(left_lines, lines)
            barcodes <- lapply(vec_list,function(x) x[1])
            barcodes <- unlist(lapply(strsplit(unlist(barcodes),':'), function(x) x[1]))
            bs <- unique(barcodes)
            ed <- min(length(bs) - 1,10000)
            if(endflag){
                ed <- length(bs)
            }
            nb_barcode <- nb_barcode + ed
            print(nb_barcode)
            vec_left_list <- vec_list[!(barcodes %in% bs[1:ed])]
            left_lines <- lines[!(barcodes %in% bs[1:ed])]
            if(ed>=1){
                stat_rs <- lapply(1:ed, function(i){
                    b <- bs[i]
                    if(!file.exists(paste0(desOutput,'.',b,'.sam'))){
                        write(headLine,file = paste0(desOutput,'.',b,'.sam'), append = TRUE, sep = "\n")
                    }
                    write(lines[barcodes==b],file = paste0(desOutput,'.',b,'.sam'), append = TRUE, sep = "\n")
                    asBam(paste0(desOutput,'.',b,'.sam'),overwrite=TRUE)
                    file.remove(paste0(desOutput,'.',b,'.sam'))
                    rs <- readGAlignmentPairs(paste0(desOutput,'.',b,'.bam'),param=ScanBamParam(what=scanBamWhat()),use.names = T)
                    file.remove(paste0(desOutput,'.',b,'.bam'))
                    file.remove(paste0(desOutput,'.',b,'.bam.bai'))
                    st1 <- start(first(rs))
                    ed1 <- end(first(rs))
                    st2 <- start(second(rs))
                    ed2 <- end(second(rs))
                    st <- st1
                    ed <- ed2
                    sel1 <- as.logical(strand(first(rs)) == '-')
                    st[sel1] <- st2[sel1]
                    ed[sel1] <- st1[sel1]
                    bed <- paste(rname(first(rs)),st,ed,sep='\t')
                    sel <- duplicated(bed)
                    rtracklayer::export(rs[!sel],BamFile(paste0(desOutput,'.',b,'.unique.sam')))
                    bed1 <- as.data.frame(table(bed))
                    return(list(tsv=paste(bed1$bed,b, bed1$Freq,sep='\t'),
           #                 bed= rs[!sel],
           #                 bf = BamFile(uniqueBams[i]),
                            stat=data.frame(barcode=b,
                                  total = length(rs),
                                  duplicate = sum(bed1$Freq>1),
                                  chimeric = sum(mcols(first(rs))$flag==2048 | mcols(second(rs))$flag==2048),
                                  unmapped = sum(mcols(first(rs))$flag==4 | mcols(second(rs))$flag==4),
                                  lowmapq = sum(mcols(first(rs))$mapq < 30 | mcols(second(rs))$mapq < 30),
                                  mitochondrial = sum(rname(first(rs))=='chrM' | rname(second(rs))=='chrM'),
                                  nonprimary = sum(mcols(first(rs))$flag==256 | mcols(second(rs))$flag==256),
                                  passed_filters = sum(!sel))))
                })
                if(ed==1){
                    file.rename(from=paste0(desOutput,'.',bs[1:ed],'.unique.bam'),to=paste0(desOutput,'.',count,'.unique.bam'))
                }else{
                    mergeBam(files=paste0(desOutput,'.',bs[1:ed],'.unique.bam'),destination=paste0(desOutput,'.',count,'.unique.bam'))
                }
                file.remove(paste0(desOutput,'.',bs[1:ed],'.unique.bam'))
                file.remove(paste0(desOutput,'.',bs[1:ed],'.unique.bam.bai'))
                count <- count + 1
                tsv <- write(unlist(lapply(stat_rs, function(v){v$tsv})),file=tsvOutput, append = TRUE, sep = "\n")
                stat_list <- c(stat_list, lapply(stat_rs, function(v){v$stat}))
            }
            if(endflag){
                break
            }
            if(length(bs) < 2 * 10000){
                lines <- readLines(samfile,n=10000000)
            }else{
                lines <- readLines(samfile,n=100)
            }
                 if(length(lines)==0 && length(bs) == 0){
                    break
                }
                if(length(lines)==0 && length(bs) < 2*10000){
                    endflag <- TRUE
                }
        }
        Rsamtools::mergeBam(files=paste0(desOutput,'.',1:(count-1),'.unique.bam'), destination=uniqueBamOutput,overwrite = TRUE)
        Rsamtools::indexBam(file=uniqueBamOutput)
        write.csv(do.call(rbind,stat_list),file=statCsvOutput,quote = FALSE, row.names=FALSE)
        file.remove(paste0(desOutput,'.',1:(count-1),'.unique.bam'))
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "SCSamToBam",
    definition = function(.Object, ...){
        .Object
    }
)


#' @name SCSamToBam
#' @title Convert sam format to bam format.
#' @description
#' This function convert a sam file into a bam file.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacBowtie2Mapping}}.
#' @param samInput \code{Character} scalar.
#' Sam file input path.
#' @param bamOutput \code{Character} scalar.
#' Bam file output path. If ignored, bed file will be put in the same path as
#' the sam file.
#' @param isSort \code{Logical} scalar.
#' Sort bam.
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{bamToBed} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' sam_bz <- system.file("extdata", "Example.sam.bz2", package="esATAC")
#' sam_path <- as.vector(bunzip2(filename = sam_bz,
#' destname = file.path(getwd(), "Example.sam"),
#' ext="bz2", FUN=bzfile, remove = FALSE))
#' sam2bam(samInput = sam_path)
#'
#' @seealso
#' \code{\link{atacBowtie2Mapping}}
#' \code{\link{atacBam2Bed}}
#' \code{\link{atacBamSort}}


setGeneric("atacSCSam2Bam",function(atacProc,
                                  samInput = NULL, bamOutput = NULL,uniqueBamOutput = NULL, tsvOutput = NULL, statCsvOutput = NULL, isSort=TRUE, ...) standardGeneric("atacSCSam2Bam"))
#' @rdname SCSamToBam
#' @aliases atacSCSam2Bam
#' @export
setMethod(
    f = "atacSCSam2Bam",
    signature = "ATACProc",
    definition = function(atacProc,
                          samInput = NULL, bamOutput = NULL, uniqueBamOutput = NULL, tsvOutput = NULL, statCsvOutput = NULL,isSort=TRUE, ...){
        allpara <- c(list(Class = "SCSamToBam", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname SCSamToBam
#' @aliases scSam2bam
#' @export
scSam2bam <- function(samInput, bamOutput = NULL, uniqueBamOutput = NULL, tsvOutput = NULL, statCsvOutput = NULL, isSort=TRUE, ...){
    allpara <- c(list(Class = "SCSamToBam", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
