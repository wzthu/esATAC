setClass(Class = "BamToBed",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "BamToBed",
    definition = function(.Object, prevSteps = list(),...){
        allparam <- list(...)
        bamInput <- allparam[["bamInput"]]
        bedOutput <- allparam[["bedOutput"]]
        reportOutput <- allparam[["reportOutput"]]
        mergePairIntoFrag <- allparam[["mergePairIntoFrag"]]
        posOffset <- allparam[["posOffset"]]
        negOffset <- allparam[["negOffset"]]
        chrFilterList <- allparam[["chrFilterList"]]
        sortBed <- allparam[["sortBed"]]
        rmMultiMap <- allparam[["rmMultiMap"]]
        minFragLen <- allparam[["minFragLen"]]
        maxFragLen <- allparam[["maxFragLen"]]
        saveExtLen <- allparam[["saveExtLen"]]
        uniqueBed <- allparam[["uniqueBed"]]
        bsgenome <- allparam[["bsgenome"]]
       
        if(length(prevSteps) > 0){
            prevSteps <- prevSteps[[1]]
            input(.Object)[["bamInput"]] <- output(prevSteps)[["bamOutput"]]
        }else{
            input(.Object)[["bamInput"]] <- bamInput
        }

        if(is.null(bedOutput)){
            output(.Object)[["bedOutput"]] <- getAutoPath(.Object, input(.Object)[["bamInput"]],"bam","bed")
        }else{
            output(.Object)[["bedOutput"]] <- addFileSuffix(bedOutput, ".bed")
        }
        
        if(is.null(reportOutput)){
            if(!is.null(input(.Object)[["bamInput"]])){
                output(.Object)[["reportOutput"]] <- getAutoPath(.Object, input(.Object)[["bamInput"]], "BAM|bam|Bam",".report.txt")
            }
        }else{
            output(.Object)[["reportOutput"]] <- reportOutput;
        }
        
        param(.Object)[["mergePairIntoFrag"]] <-  match.arg(mergePairIntoFrag, c("auto","yes","no"))
        param(.Object)[["posOffset"]] <- posOffset
        param(.Object)[["negOffset"]] <- negOffset
        param(.Object)[["filterList"]] <- chrFilterList
        param(.Object)[["sortBed"]] <- sortBed
        param(.Object)[['rmMultiMap']] <- rmMultiMap
        param(.Object)[["uniqueBed"]] <- match.arg(uniqueBed, c("auto","yes","no"))
        param(.Object)[["minFragLen"]] <- minFragLen
        param(.Object)[["maxFragLen"]] <- maxFragLen
        param(.Object)[["saveExtLen"]] <- saveExtLen
        if(is.null(bsgenome)){
            param(.Object)[['bsgenome']] <- getRefRc('bsgenome')
        }else{
            param(.Object)[['bsgenome']] <- bsgenome
        }
        
        .Object
    } # definition end
) # setMethod initialize end


#' @importFrom Rsamtools testPairedEndBam scanBam ScanBamParam
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom IRanges IRanges
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignmentPairs
#' @importFrom rtracklayer import export import.bed export.bed

setMethod(
    f = "processing",
    signature = "BamToBed",
    definition = function(.Object, ...){
        bamInput <- input(.Object)[["bamInput"]]
        mergePairIntoFrag <- param(.Object)[['mergePairIntoFrag']]
        uniqueBed <- param(.Object)[['uniqueBed']]
        isPaired <- Rsamtools::testPairedEndBam(bamInput)
        if(isPaired){
            if(mergePairIntoFrag == 'auto' || mergePairIntoFrag == 'yes'){
                mergePairIntoFrag <- TRUE
            }else{
                mergePairIntoFrag <- FALSE
            }
            if(uniqueBed == 'auto' || uniqueBed == 'yes'){
                uniqueBed <- TRUE
            }else{
                uniqueBed <- FALSE
            }
        }else{
            if(mergePairIntoFrag == 'auto' || mergePairIntoFrag == 'no'){
                mergePairIntoFrag <- FALSE
            }else{
                stop('mergePairIntoFrag can not be yes: single end reads can not be merged!')
            }
            if(uniqueBed == 'auto' || uniqueBed == 'no'){
                uniqueBed <- FALSE
            }else{
                uniqueBed <- TRUE
            }
        }
        
        allchr <- seqnames(param(.Object)[['bsgenome']])
        chrlen <- seqlengths(param(.Object)[['bsgenome']])
        
        allchr <- allchr[grep(param(.Object)[["filterList"]],allchr,invert = TRUE)]
        chrlen <- chrlen[allchr]
        
        totalReads <- length(scanBam(file=bamInput, param =ScanBamParam(what=c('strand')))[[1]]$strand)
        mutiMapReads <- 0
        extLenReads <- 0
        pcrReads <- 0
        chrMReads <- 0
        filteredReads <- 0
        saveReads <- 0
        
        readsAfterFiltered <- 0
        
        allchrWithM <- c(allchr,'chrM')
        
        originBedDir <- getStepWorkDir(.Object,'origin')
        cleanBedDir <- getStepWorkDir(.Object,'clean')
        
        dir.create(originBedDir)
        dir.create(cleanBedDir)
        
        if(mergePairIntoFrag){
            XS <- lapply(allchrWithM, function(x){
                s<-list()
                s[[x]] <- IRanges(start = 1,end = 536870912)
                if(param(.Object)[['rmMultiMap']]){
                    p <- ScanBamParam(which=do.call(IRangesList,s), tag = c('XS'))
                    obj<-readGAlignmentPairs(file = bamInput,param = p)
                    rtracklayer::export.bed(obj,
                                            con = file.path(originBedDir,paste0(x,'.bed')))
                    f <- mcols(first(obj))$XS
                    s <- mcols(second(obj))$XS
                    return(is.na(f) | is.na(s))
                }else{
                    p <- ScanBamParam(which=do.call(IRangesList,s))
                    rtracklayer::export.bed(readGAlignmentPairs(file = bamInput,param = p),
                                            con = file.path(originBedDir,paste0(x,'.bed')))
                }
                
            })
            
            chrMReads <- length(import.bed(con = file.path(originBedDir,"chrM.bed")))
            
            names(XS) <- allchrWithM
            for(x in allchr){
                gr <-import.bed(con = file.path(originBedDir,paste0(x,'.bed')))
                
                readsAfterFiltered <- readsAfterFiltered + length(gr)
                
                if(param(.Object)[['rmMultiMap']]){
                    gr <- gr[XS[[x]]]
                    mutiMapReads <- mutiMapReads + sum(!XS[[x]])
                }
                start(gr[strand(gr)=='+']) <- start(gr[strand(gr)=='+']) + param(.Object)[["posOffset"]]
                end(gr[strand(gr)=='+']) <- end(gr[strand(gr)=='+']) + param(.Object)[["posOffset"]]
                start(gr[strand(gr)=='-']) <- start(gr[strand(gr)=='-']) + param(.Object)[["negOffset"]]
                end(gr[strand(gr)=='-']) <- end(gr[strand(gr)=='-']) + param(.Object)[["negOffset"]]
                
                extLenReads <- extLenReads + sum(width(gr)<param(.Object)[["minFragLen"]] | 
                                                 width(gr)>param(.Object)[["maxFragLen"]])
                
                gr <- gr[width(gr)>=param(.Object)[["minFragLen"]] & 
                             width(gr)<=param(.Object)[["maxFragLen"]] ]
                
                
                
                if(uniqueBed){
                    beforeUnique <- length(gr)
                    gr <- unique(gr)
                    pcrReads <- pcrReads + beforeUnique - length(gr)
                }
                
                
                saveReads <-saveReads + length(gr)
                export.bed(gr, con = file.path(cleanBedDir,paste0(x,'.bed')))
            }
            
            mutiMapReads <- mutiMapReads * 2
            extLenReads <- extLenReads * 2
            pcrReads <- pcrReads * 2
            chrMReads <- chrMReads * 2
            filteredReads <- totalReads - readsAfterFiltered * 2
            saveReads <- saveReads * 2
        }else{
            XS <- lapply(allchrWithM, function(x){
                s<-list()
                s[[x]] <- IRanges(start = 1,end = 536870912)
                if(param(.Object)[['rmMultiMap']]){
                    p <- ScanBamParam(which=do.call(IRangesList,s), tag = c('XS'))
                    obj<-readGAlignments(file = bamInput,param = p)
                    rtracklayer::export.bed(obj,
                                            con = file.path(originBedDir,paste0(x,'.bed')))
                    f <- mcols(first(obj))$XS
                    s <- mcols(second(obj))$XS
                    return(is.na(f) | is.na(s))
                }else{
                    p <- ScanBamParam(which=do.call(IRangesList,s))
                    rtracklayer::export.bed(readGAlignments(file = bamInput,param = p),
                                            con = file.path(originBedDir,paste0(x,'.bed')))
                }
                
            })
            
            chrMReads <- length(import.bed(con = file.path(originBedDir,"chrM.bed")))
            
            names(XS) <- allchrWithM
            for(x in allchr){
                gr <-import.bed(con =  file.path(originBedDir,paste0(x,'.bed')))
                
                readsAfterFiltered <- readsAfterFiltered + length(gr)
                
                if(param(.Object)[['rmMultiMap']]){
                    gr <- gr[XS[[x]]]
                    mutiMapReads <- mutiMapReads + sum(!XS[[x]])
                }
                start(gr[strand(gr)=='+']) <- start(gr[strand(gr)=='+']) + param(.Object)[["posOffset"]]
                end(gr[strand(gr)=='-']) <- end(gr[strand(gr)=='-']) + param(.Object)[["negOffset"]]
                
                extLenReads <- extLenReads + sum(width(gr)<param(.Object)[["minFragLen"]] | 
                                                     width(gr)>param(.Object)[["maxFragLen"]])
                
                gr <- gr[width(gr)>=param(.Object)[["minFragLen"]] & 
                             width(gr)<=param(.Object)[["maxFragLen"]] ]
                
                
                
                if(uniqueBed){
                    beforeUnique <- length(gr)
                    gr <- unique(gr)
                    pcrReads <- pcrReads + beforeUnique - length(gr)
                }
                
                
                saveReads <-saveReads + length(gr)
                export.bed(gr, con =  file.path(cleanBedDir,paste0(x,'.bed')))
            }
            filteredReads <- totalReads - readsAfterFiltered 
        }
        if(file.exists(output(.Object)[['bedOutput']])){
            file.remove(output(.Object)[['bedOutput']])
        }
        if(param(.Object)[["sortBed"]]){
            tmpbed <- paste0(output(.Object)[['bedOutput']],'.tmp.bed')
            if(file.exists(tmpbed)){
                file.remove(tmpbed)
            }
            lapply(allchr, function(x){
                file.append(tmpbed, file.path(cleanBedDir,paste0(x,'.bed')))
            })
            export.bed(sort(import.bed(tmpbed)),output(.Object)[['bedOutput']])
            file.remove(tmpbed)
        }else{
            lapply(allchr, function(x){
                file.append(output(.Object)[['bedOutput']], file.path(cleanBedDir,paste0(x,'.bed')))
            })
        }
        
        
        
        write.table(data.frame(total=totalReads,
                               save=saveReads,
                               filted=filteredReads,
                               extlen=extLenReads,
                               unique = pcrReads,
                               multimap=mutiMapReads,
                               chrM = chrMReads),
                    file = output(.Object)[["reportOutput"]],
                    quote=FALSE,sep="\t",row.names=FALSE)
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "BamToBed",
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
        report(.Object)[["non-mitochondrial"]] <- as.character(as.integer(qcval[["total"]])-as.integer(qcval[["chrM"]]))
        report(.Object)[["non-filtered-chr"]] <- as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]))
        report(.Object)[["non-mitochondrial-multimap"]] <- as.character(as.integer(qcval[["total"]])-as.integer(qcval[["filted"]]) -as.integer(qcval[["multimap"]]))
        for(n in names(qcval)){
            report(.Object)[[n]] <- qcval[[n]]
        }  
        .Object
    }
)




#' @name BamToBed
#' @title Convert bam format to bed format.
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
#' \code{\link{atacBamSort}},
#' \code{\link{atacSam2Bam}}.
#' @param bamInput \code{Character} scalar.
#' Bam file input path.
#' @param bedOutput \code{Character} scalar.
#' Bed file output path. If ignored, bed file will be put in the same path as
#' the bam file.
#' @param reportOutput \code{Character} scalar.
#' Report file path.
#' @param bsgenome \code{BSgenome} object.
#' This object from bioconductor
#' @param mergePairIntoFrag \code{Logical} scalar
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
#' @param saveExtLen \code{Logical} scaler.
#' Save the fragment that are not in the range of minFragLen and maxFragLen
#' @param rmMultiMap \code{Logical} scalar.
#' Remove multi-map reads.
#' @param ... Additional arguments, currently unused.
#' @details The bam file wiil be automatically obtained from
#' object(\code{atacProc}) or input by hand. Output can be ignored.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Zheng Wei, Wei Zhang
#' @examples
#'
#' library(Rsamtools)
#' ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' bam2bed(bamInput = ex1_file)
#'
#' @seealso
#' \code{\link{atacBamSort}}
#' \code{\link{atacSam2Bam}}
#' @importFrom rtracklayer export



setGeneric("atacBam2Bed", function(atacProc, bamInput = NULL, bedOutput = NULL, 
                                   reportOutput =NULL, bsgenome = NULL,
                                   mergePairIntoFrag = c("auto","yes","no"), 
                                   posOffset = +4, negOffset= -5, 
                                   chrFilterList= "chrM|_",
                                   sortBed = TRUE, rmMultiMap=TRUE,
                                   minFragLen = 0,maxFragLen = 2000,
                                   saveExtLen = FALSE,uniqueBed = c("auto","yes","no"), ...) standardGeneric("atacBam2Bed"))

#' @rdname BamToBed
#' @aliases atacBam2Bed
#' @export
setMethod(
    f = "atacBam2Bed",
    signature = "ATACProc",

    definition = function(atacProc, bamInput = NULL, bedOutput = NULL, 
                          reportOutput =NULL, bsgenome = NULL,
                          mergePairIntoFrag = c("auto","yes","no"), 
                          posOffset = +4, negOffset= -5, 
                          chrFilterList= "chrM|_",#chrUn.*|chrM|.*random.*
                          sortBed = TRUE, rmMultiMap=TRUE,
                          minFragLen = 0,maxFragLen = 2000,
                          saveExtLen = FALSE,uniqueBed = c("auto","yes","no"),...){
        allpara <- c(list(Class = "BamToBed", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname BamToBed
#' @aliases bam2bed
#' @export

bam2bed <- function(bamInput, bedOutput = NULL, 
                    reportOutput =NULL, bsgenome = NULL,
                    mergePairIntoFrag = c("auto","yes","no"), 
                    posOffset = +4, negOffset= -5, 
                    chrFilterList= "chrM|_",#chrUn.*|chrM|.*random.*
                    sortBed = TRUE, rmMultiMap=TRUE,
                    minFragLen = 0,maxFragLen = 2000,
                    saveExtLen = FALSE,uniqueBed = c("auto","yes","no"), ...){
    allpara <- c(list(Class = "BamToBed", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
