setClass(Class = "SCCellFilter",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "SCCellFilter",
    definition = function(.Object, prevSteps = list(), ...){
        allparam <- list(...)
        csvInput <- allparam[["csvInput"]]
        
        threshold_fragQC1 <- allparam[["threshold_fragQC1"]]
        threshold_fragQC2 <- allparam[["threshold_fragQC2"]]
        threshold_nsQC <- allparam[["threshold_nsQC"]]
        threshold_tssQC <- allparam[["threshold_tssQC"]]
        threshold_peakQC1 <- allparam[["threshold_peakQC1"]]
        threshold_peakQC2 <- allparam[["threshold_peakQC2"]]
        threshold_peakQC3 <- allparam[["threshold_peakQC3"]]
        threshold_blacklistQC <- allparam[["threshold_blacklistQC"]]
        
        atacProc <- NULL
        if(length(prevSteps) > 0){
            atacProc <- prevSteps[[1]]
        }
        
        # necessary parameters
        if(is.null(atacProc)){
            input(.Object)[["csvInput"]] <- csvInput
        }else{
            input(.Object)[["csvInput"]] <- getParam(atacProc, "csvOutput")
            property(.Object)[["peak_bc_matrix.h5"]] <- getParam(atacProc, "peak_bc_matrix.h5")
        }
        
        if (is.null(threshold_fragQC1)) {
            param(.Object)[["threshold_fragQC1"]] <- 1000
        }else{
            param(.Object)[["threshold_fragQC1"]] <- threshold_fragQC1
        }
        
        if (is.null(threshold_fragQC2)) {
            param(.Object)[["threshold_fragQC2"]] <- 45000
        }else{
            param(.Object)[["threshold_fragQC2"]] <- threshold_fragQC2
        }
        
        if (is.null(threshold_nsQC)) {
            param(.Object)[["threshold_nsQC"]] <- 4
        }else{
            param(.Object)[["threshold_nsQC"]] <- threshold_nsQC
        }
        
        if (is.null(threshold_tssQC)) {
            param(.Object)[["threshold_tssQC"]] <- 2
        }else{
            param(.Object)[["threshold_tssQC"]] <- threshold_tssQC
        }
        
        if (is.null(threshold_peakQC1)) {
            param(.Object)[["threshold_peakQC1"]] <- 3000
        }else{
            param(.Object)[["threshold_peakQC1"]] <- threshold_peakQC1
        }
        
        if (is.null(threshold_peakQC2)) {
            param(.Object)[["threshold_peakQC2"]] <- 20000
        }else{
            param(.Object)[["threshold_peakQC2"]] <- threshold_peakQC2
        }
        
        if (is.null(threshold_peakQC3)) {
            param(.Object)[["threshold_peakQC3"]] <- 0.15
        }else{
            param(.Object)[["threshold_peakQC3"]] <- threshold_peakQC3
        }
        
        if (is.null(threshold_blacklistQC)) {
            param(.Object)[["threshold_blacklistQC"]] <- 0.05
        }else{
            param(.Object)[["threshold_blacklistQC"]] <- threshold_blacklistQC
        }
        
        # output
        output(.Object)[["csvOutput"]] <- getStepWorkDir(.Object, filename = "singlecell.csv")
        
        .Object
    }
)



#' @importFrom data.table fread data.table fwrite
setMethod(
    f = "processing",
    signature = "SCCellFilter",
    definition = function(.Object, ...){
        
        ## reading fragment file and create fragment object
        print("Reading single cell information......")
        metadata <- fread(file = input(.Object)[["csvInput"]])
        
        ## filtering total counts
        valid_frag <- as.numeric(metadata$passed_filters)
        tmp_meta <- metadata[which((valid_frag > param(.Object)[["threshold_fragQC1"]]) & (valid_frag < param(.Object)[["threshold_fragQC2"]])), ]
        print(paste0(nrow(tmp_meta), " cells passed fragment QC."))
        
        ## filtering nsQC
        nsQC <- as.numeric(tmp_meta$nucleosome_signal)
        tmp_meta <- tmp_meta[which(nsQC < param(.Object)[["threshold_nsQC"]]), ]
        print(paste0(nrow(tmp_meta), " cells passed nucleosome QC."))
        
        ## filtering tssQC
        tssQC <- as.numeric(tmp_meta$TSS.enrichment)
        tmp_meta <- tmp_meta[which(tssQC > param(.Object)[["threshold_tssQC"]]), ]
        print(paste0(nrow(tmp_meta), " cells passed TSS QC."))
        
        ## filtering peakQC
        peakCount <- as.numeric(tmp_meta$frag_in_peaks)
        tmp_meta <- tmp_meta[which((peakCount > param(.Object)[["threshold_peakQC1"]]) & (peakCount < param(.Object)[["threshold_peakQC2"]])), ]
        peakCount <- as.numeric(tmp_meta$frag_in_peaks)
        validCount <- as.numeric(tmp_meta$passed_filters)
        ratioPeak <- peakCount / validCount
        tmp_meta <- tmp_meta[which(ratioPeak > param(.Object)[["threshold_peakQC3"]]), ]
        print(paste0(nrow(tmp_meta), " cells passed Peak QC."))
        
        
        ## filtering blacklistQC
        blackCount <- as.numeric(tmp_meta$frag_in_blacklist)
        peakCount <- as.numeric(tmp_meta$frag_in_peaks)
        ratioBlack <- blackCount / peakCount
        tmp_meta <- tmp_meta[which(ratioBlack < param(.Object)[["threshold_blacklistQC"]]), ]
        print(paste0(nrow(tmp_meta), " cells passed blacklist QC."))
        
        ## update csv
        fwrite(x = tmp_meta, file = output(.Object)[["csvOutput"]], na = "Na")
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "SCCellFilter",
    definition = function(.Object, ...){
        report(.Object)[["csvOutput"]] <- output(.Object)[["csvOutput"]]
        .Object
    }
)

#' @name SCCellFilter
#'
#' @title Filtering valid single cells
#'
#' @description
#' Filtering valid single cells for downstream aalysis.
#'
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{SCQC}}
#' 
#' @param csvInput scATAC-seq csv record file.
#' 
#' @param threshold_fragQC1 Lower bound for fragment count per cell.
#'
#' @param threshold_fragQC2 Upper bound for fragment count per cell.
#' 
#' @param threshold_nsQC Nucleosome QC threshold (less than).
#' 
#' @param threshold_tssQC TSS QC threshold (great than).
#' 
#' @param threshold_peakQC1 Lower bound for fragment count in peaks.
#' 
#' @param threshold_peakQC2 Upper bound for fragment count in peaks.
#' 
#' @param threshold_peakQC3 peak counts ratio (great than).
#' 
#' @param threshold_blacklistQC blacklist counts ratio (less than).
#' 
#' @param ... Additional arguments, currently unused.
#'
#' @return An invisible \code{\link{ATACProc-class}} object scalar.
#'
#' @author Wei Zhang
#'
#' @examples
#' print(123)
setGeneric("atacSCCellFilter", function(atacProc, 
                                        csvInput = NULL,
                                        threshold_fragQC1 = 1000,
                                        threshold_fragQC2 = 45000,
                                        threshold_nsQC = 4,
                                        threshold_tssQC = 2,
                                        threshold_peakQC1 = 3000, 
                                        threshold_peakQC2 = 20000,
                                        threshold_peakQC3 = 0.15,
                                        threshold_blacklistQC = 0.05, ...) standardGeneric("atacSCCellFilter"))

#' @rdname SCCellFilter
#' @aliases atacSCCellFilter
#' @export
setMethod(
    f = "atacSCCellFilter",
    signature = "ATACProc",
    definition = function(atacProc, 
                          csvInput = NULL,
                          threshold_fragQC1 = 1000,
                          threshold_fragQC2 = 45000,
                          threshold_nsQC = 4,
                          threshold_tssQC = 2,
                          threshold_peakQC1 = 3000, 
                          threshold_peakQC2 = 20000,
                          threshold_peakQC3 = 0.15,
                          threshold_blacklistQC = 0.05, ...){
        allpara <- c(list(Class = "SCCellFilter", prevSteps = list(atacProc)),
                     as.list(environment()), list(...))
        step <- do.call(new, allpara)
        invisible(step)
    }
)


#' @rdname SCCellFilter
#' @aliases atacSCCellFilter
#' @export
atacscCellFilter <- function(csvInput = NULL,
                             threshold_fragQC1 = 1000,
                             threshold_fragQC2 = 45000,
                             threshold_nsQC = 4,
                             threshold_tssQC = 2,
                             threshold_peakQC1 = 3000, 
                             threshold_peakQC2 = 20000,
                             threshold_peakQC3 = 0.15,
                             threshold_blacklistQC = 0.05, ...){
    
    allpara <- c(list(Class = "SCCellFilter", prevSteps = list()), as.list(environment()), list(...))
    step <- do.call(new,allpara)
    invisible(step)
}




























