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
        samfile <- file(samInput,'r')
        headLine <- c()
        while(TRUE){
            lines <- readLines(samfile,n=1)
            if(length(lines)==0||!startsWith(lines,'@')){
                break
            }else{
                headLine <- c(headLine,lines)
            }
        }
        lines <- c(lines, readLines(samfile,n=10000000))
        barcode <- NULL
        while(length(lines)>0){
            lb <- unlist(lapply(strsplit(lines,':'), function(x) x[1]))
            barcode<-c(barcode,unique(lb))
            lapply(barcode, function(b){
                if(!file.exists(paste0(desOutput,'.',b,'.sam'))){
                    write(headLine,file = paste0(desOutput,'.',b,'.sam'), append = TRUE, sep = "\n")
                }
                write(lines[lb==b],file = paste0(desOutput,'.',b,'.sam'), append = TRUE, sep = "\n")
            })
            lines <- readLines(samfile,n=10000000)
        }
        lapply(barcode, function(b){
            asBam(paste0(desOutput,'.',b,'.sam'),overwrite=TRUE)
        })
        bams <- paste0(desOutput,'.',barcode,'.bam')
        uniqueBams <- paste0(desOutput,'.',barcode,'.uniqued.bam')
        batch_size <- 10000
        cl <- makeCluster(28)
        statrs <- lapply(seq(1,length(barcode), batch_size),function(batch){
            tsv_csv <- parLapply(cl = cl, X=batch:min(batch+batch_size-1,length(barcode)), 
                fun=function(i,bams,readGAlignmentPairs ,ScanBamParam, scanBamWhat,start,end ,first,second,strand,mcols,rname,BamFile,uniqueBams,barcode){
                Sys.sleep(runif(10))
                rs <- readGAlignmentPairs(bams[i],param=ScanBamParam(what=scanBamWhat()),use.names = T)
            #    rtracklayer::export(rs,BamFile(bams[i]))
            #    asSam(file=bams[i])
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
          #      rtracklayer::export(rs[!sel],BamFile(uniqueBams[i]))
                bed1 <- as.data.frame(table(bed))
#                write(paste(bed1$bed,barcode[i], bed1$Freq,sep='\t'),file=tsvOutput, append = TRUE, sep = "\n")
                return(list(tsv=paste(bed1$bed,barcode[i], bed1$Freq,sep='\t'),
                            bed= rs[!sel],
                            bf = BamFile(uniqueBams[i]),
                            stat=data.frame(barcode=barcode[i],
                                  total = length(rs),
                                  duplicate = sum(bed1$Freq>1),
                                  chimeric = sum(mcols(first(rs))$flag==2048 | mcols(second(rs))$flag==2048),
                                  unmapped = sum(mcols(first(rs))$flag==4 | mcols(second(rs))$flag==4),
                                  lowmapq = sum(mcols(first(rs))$mapq < 30 | mcols(second(rs))$mapq < 30),
                                  mitochondrial = sum(rname(first(rs))=='chrM' | rname(second(rs))=='chrM'),
                                  nonprimary = sum(mcols(first(rs))$flag==256 | mcols(second(rs))$flag==256),
                                  passed_filters = sum(!sel))))
            },bams = bams,
                readGAlignmentPairs = readGAlignmentPairs,
                ScanBamParam = ScanBamParam,
                scanBamWhat = scanBamWhat,
                start = start,
                end = end,
                first = first,
                second = second,
                strand = strand,
                mcols = mcols,
                rname = rname,
                BamFile = BamFile,
                uniqueBams = uniqueBams,
                barcode = barcode)
            tsv <- write(unlist(lapply(tsv_csv, function(v){v$tsv})),file=tsvOutput, append = TRUE, sep = "\n")
            csv <- do.call(rbind,lapply(tsv_csv, function(v){v$stat}))
            lapply(tsv_csv, function(v){ rtracklayer::export(v$bed,v$bf)})
            return(csv)
        })
        write.csv(do.call(rbind,statrs),file=statCsvOutput,quote = FALSE, row.names=FALSE)
        mergeBam(files=uniqueBams, destination=uniqueBamOutput,overwrite = TRUE)
        indexBam(file=uniqueBamOutput)
        file.remove(uniqueBams)
        file.remove(paste0(uniqueBams,'.bai'))
        file.remove(bams)
        file.remove(paste0(bams,'.bai'))
        file.remove(paste0(desOutput,'.',barcode,'.sam'))
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
