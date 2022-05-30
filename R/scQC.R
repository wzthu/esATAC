setClass(Class = "SCQC",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "SCQC",
    definition = function(.Object, prevSteps = list(), ...){
        allparam <- list(...)
        fragInput <- allparam[["fragInput"]]
        csvInput <- allparam[["csvInput"]]
        
        peak <- allparam[["peak"]]
        gene.annotation <- allparam[["gene.annotation"]]
        blacklist <- allparam[["blacklist"]]
        n <- allparam[["n"]]
        
        
        
        param(.Object)[["n"]] <- allparam[["n"]]
        param(.Object)[["gene.annotation"]] <- allparam[["gene.annotation"]]
        param(.Object)[["peak"]] <- allparam[["peak"]]
        param(.Object)[["blacklist"]] <- allparam[["blacklist"]]
        
        
        atacProcFrag <- NULL
        if(length(prevSteps) > 0){
            atacProcFrag <- prevSteps[[1]]
        }
        
        atacProcPeak <- NULL
        if(length(prevSteps) > 0){
            atacProcPeak <- prevSteps[[2]]
        }
        
        
        # necessary parameters
        if(is.null(atacProcFrag)){
            input(.Object)[["fragInput"]] <- fragInput
            input(.Object)[["csvInput"]] <- csvInput
        }else{
            input(.Object)[["fragInput"]] <- getParam(atacProcFrag, "fragOutput")
            input(.Object)[["csvInput"]] <- getParam(atacProcFrag, "csvOutput")
        }
        
        if(is.null(atacProcPeak)){
            param(.Object)[["peak"]] <- peak
        }else{
            param(.Object)[["peak"]] <- getParam(atacProcPeak, "peakOUtput")
        }
        
        if (is.null(gene.annotation)) {
            param(.Object)[["gene.annotation"]] <- getRefFiles("EnsDb.Hsapiens")
        }else{
            param(.Object)[["gene.annotation"]] <- gene.annotation
        }
        
        if (is.null(blacklist)) {
            param(.Object)[["blacklist"]] <- getRefFiles("blacklist")
        }else{
            param(.Object)[["blacklist"]] <- blacklist
        }
        
        param(.Object)[["n"]] <- n
        
        # output
        output(.Object)[["nucleosomeQC.rds"]] <- getStepWorkDir(.Object, filename = "nucleosomeQC.rds")
        output(.Object)[["nucleosomeQC.pdf"]] <- getStepWorkDir(.Object, filename = "nucleosomeQC.pdf")
        
        output(.Object)[["tssQC.rds"]] <- getStepWorkDir(.Object, filename = "tssQC.rds")
        output(.Object)[["tssQC.pdf"]] <- getStepWorkDir(.Object, filename = "tssQC.pdf")
        
        output(.Object)[["peak_bc_matrix.h5"]] <- getStepWorkDir(.Object, filename = "peak_bc_matrix.h5")
        output(.Object)[["peak_bc_matrix.rds"]] <- getStepWorkDir(.Object, filename = "peak_bc_matrix.rds")
        output(.Object)[["fragInPeaks.rds"]] <- getStepWorkDir(.Object, filename = "fragInPeaks.rds")
        
        output(.Object)[["fragInBlacklist.rds"]] <- getStepWorkDir(.Object, filename = "fragInBlacklist.rds")
        
        output(.Object)[["Violin.pdf"]] <- getStepWorkDir(.Object, filename = "Violin.pdf")
        
        .Object
    }
)

#' @importFrom rtracklayer import
#' @importFrom ggplot2 ggsave
#' @importFrom Matrix colSums
setMethod(
    f = "processing",
    signature = "SCQC",
    definition = function(.Object,...){
        ## reading fragment file and create fragment object
        print("Reading scATAC-seq fragment file......")
        fragment <- fragCreate(fragment = input(.Object)[["fragInput"]], 
                               csv = input(.Object)[["csvInput"]])
        
        cells <- as.character(fragment@cells)
        
        # nucleosome QC
        print("Now, processing scATAC-seq nucleosome QC......")
        nucleosomeQC <- scNucleosomeQC(frags = fragment, n = NULL)
        
        print("Saving scATAC-seq nucleosome QC results......")
        saveRDS(object = nucleosomeQC,
                file = output(.Object)[["nucleosomeQC.rds"]])
        
        p_NucleosomeQC <- scPlotNucleosomeQC(fragment = fragment,
                                             nsQC = nucleosomeQC)
        
        ggsave(output(.Object)[["nucleosomeQC.pdf"]], plot = p_NucleosomeQC)
        
        ## TSS QC
        print("Now, processing scATAC-seq TSS QC......")
        if (class(param(.Object)[["gene.annotation"]]) == "EnsDb") {
            print("Gene annotation detected in 'EnsDb' class, converting to GRanges......")
            annotations <- GetGRangesFromEnsDb(ensdb = param(.Object)[["gene.annotation"]])
            seqlevelsStyle(annotations) <- 'UCSC'
        }else if (class(param(.Object)[["gene.annotation"]]) == "character") {
            print("Gene annotation detected in file, reading and converting to GRanges......")
            annotations <- import(param(.Object)[["gene.annotation"]])
        }else{
            stop("Parameter 'gene.annotation' error!")
        }
        
        tssQC <- scTssQC(object = fragment, gene.annotation = annotations)
        
        print("Saving scATAC-seq TSS QC results......")
        saveRDS(object = tssQC,
                file = output(.Object)[["tssQC.rds"]])
        
        p_TssQC <- scPlotTssQC(tssInfo = tssQC)
        
        ggsave(output(.Object)[["tssQC.pdf"]], plot = p_TssQC)
        
        ## fragInPeaks QC
        print("Now, processing scATAC-seq FRiP QC......")
        print("Fetching single cell data in peaks......")
        
        if (param(.Object)[["peak"]] == "GRanges") {
            print("Gene annotation detected in 'GRanges' class......")
            peaks_GR <- param(.Object)[["peak"]]
        }else if (class(param(.Object)[["peak"]]) == "character") {
            print("Peaks detected in file, reading and converting to GRanges......")
            peaks_GR <- import(param(.Object)[["peak"]])
        }else{
            stop("Parameter 'peak' error!")
        }
        
        peak_BC_matrix <- FeatureMatrix(fragments = fragment, 
                                        features = peaks_GR, 
                                        process_n = param(.Object)[["n"]])
        
        print("Saving Peak x Barcode Matrix............")
        .write10xCounts(path = output(.Object)[["peak_bc_matrix.h5"]], x = peak_BC_matrix)
        saveRDS(object = peak_BC_matrix,
                file = output(.Object)[["peak_bc_matrix.rds"]])
        
        frag_in_peaks <- colSums(peak_BC_matrix)
        
        saveRDS(object = frag_in_peaks, 
                file = output(.Object)[["fragInPeaks.rds"]])
        
        
        ## fragInBlacklist QC
        print("Now, processing scATAC-seq fragInBlacklist QC......")
        print("Fetching single cell data in blacklist regions......")
        if (param(.Object)[["blacklist"]] == "GRanges") {
            print("Gene annotation detected in 'GRanges' class......")
            blacklist_GR <- param(.Object)[["blacklist"]]
        }else if (class(param(.Object)[["blacklist"]]) == "character") {
            print("Peaks detected in file, reading and converting to GRanges......")
            blacklist_GR <- import(param(.Object)[["blacklist"]])
        }else{
            stop("Parameter 'blacklist' error!")
        }
        
        blacklist_GR <- import(param(.Object)[["blacklist"]])
        
        blacklist_BC_matrix <- FeatureMatrix(fragments = fragment, 
                                             features = blacklist_GR, 
                                             process_n = param(.Object)[["n"]])
        
        print("Saving fragments in blacklist informations......")
        frag_in_blacklist <- colSums(blacklist_BC_matrix)
        
        saveRDS(object = frag_in_blacklist, 
                file = output(.Object)[["fragInBlacklist.rds"]])
        
        .Object
    }
)

setMethod(
    f = "genReport",
    signature = "SCQC",
    definition = function(.Object, ...){
        report(.Object)[["nucleosomeQC.rds"]] <- output(.Object)[["nucleosomeQC.rds"]]
        report(.Object)[["nucleosomeQC.pdf"]] <- output(.Object)[["nucleosomeQC.pdf"]]
        report(.Object)[["tssQC.rds"]] <- output(.Object)[["tssQC.rds"]]
        report(.Object)[["tssQC.pdf"]] <- output(.Object)[["tssQC.pdf"]]
        report(.Object)[["peak_bc_matrix.h5"]] <- output(.Object)[["peak_bc_matrix.h5"]]
        report(.Object)[["peak_bc_matrix.rds"]] <- output(.Object)[["peak_bc_matrix.rds"]]
        report(.Object)[["fragInPeaks.rds"]] <- output(.Object)[["fragInPeaks.rds"]]
        report(.Object)[["fragInBlacklist.rds"]] <- output(.Object)[["fragInBlacklist.rds"]]
        report(.Object)[["Violin.pdf"]] <- output(.Object)[["Violin.pdf"]]
        .Object
    }
)




#' @name SCQC
#'
#' @title Get Single Cell Pre-processing Information
#'
#' @description
#' Get scATAC-seq pre-processing results and create Fragment Object.
#'
#' @param atacProcFrag \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' ###################
#' 
#' @param atacProcPeak \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' ###################
#' 
#' @param fragInput scATAC-seq fragment file.
#' 
#' @param csvInput scATAC-seq csv record file.
#' 
#' @param peak scATAC-seq peak file in BED format.
#' 
#' @param gene.annotation scATAC-seq gene.annotation file in BED format.
#' Note: Please using the results from function GetGRangesFromEnsDb or 
#' using the default data.
#' 
#' @param blacklist Genome blacklist file in BED format.
#' 
#' @param n Number of regions to process at a time.
#' Default: 2000.
#' 
#' @param ... Additional arguments, currently unused.
#'
#' @return An invisible \code{\link{ATACProc-class}} object scalar.
#'
#' @author Wei Zhang
#'
#' @examples
#' print(123)
setGeneric("atacSCQC", function(atacProcFrag, 
                                atacProcPeak, 
                                fragInput = NULL,
                                csvInput = NULL,
                                peak = NULL,
                                gene.annotation = NULL,
                                blacklist = NULL, 
                                n = 2000, ...) standardGeneric("atacSCQC"))

#' @rdname SCQC
#' @aliases atacSCQC
#' @export
setMethod(
    f = "atacSCQC",
    signature = "ATACProc",
    definition = function(atacProcFrag, 
                          atacProcPeak, 
                          fragInput = NULL,
                          csvInput = NULL,
                          peak = NULL,
                          gene.annotation = NULL,
                          blacklist = NULL, 
                          n = 2000, ...){
        allpara <- c(list(Class = "SCQC", prevSteps = list(atacProc)),
                     as.list(environment()), list(...))
        step <- do.call(new, allpara)
        invisible(step)
    }
)


#' @rdname SCQC
#' @aliases atacscQC
#' @export
atacscQC <- function(fragInput = NULL,
                     csvInput = NULL,
                     peak = NULL,
                     gene.annotation = NULL,
                     blacklist = NULL, 
                     n = 2000, ...){
    
    allpara <- c(list(Class = "SCQC", prevSteps = list()), as.list(environment()), list(...))
    step <- do.call(new,allpara)
    invisible(step)
}





