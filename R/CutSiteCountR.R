setClass(Class = "CutSiteCountR",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "CutSiteCountR",
    definition = function(.Object, prevSteps = list(), ...){
        allparam <- list(...)
        atacProcCutSite  <- NULL
        atacProcMotifScan <- NULL
        if(length(prevSteps)>0){
            atacProcCutSite <- prevSteps[[1]]
            atacProcMotifScan <- prevSteps[[2]]
        }
        csInput <- allparam[["csInput"]]
        motif_info <- allparam[["motif_info"]]
        chr <- allparam[["chr"]]
        matrixOutput <- allparam[["matrixOutput"]]
        strandLength <- allparam[["strandLength"]]
        FootPrint <- allparam[["FootPrint"]]
        prefix <- allparam[["prefix"]]
        editable <- allparam[["editable"]]

        # necessary parameters
        if(!is.null(atacProcCutSite)){
            param(.Object)[["csfile.dir"]] <- getParam(atacProcCutSite, "csfile.dir");
        }else{
            param(.Object)[["csfile.dir"]] <- csInput
        }
        if(!is.null(atacProcMotifScan)){
            param(.Object)[["motif_info"]] <- readRDS(getParam(atacProcMotifScan, "rdsOutput"))
        }else{
            param(.Object)[["motif_info"]] <- readRDS(motif_info)
        }

        param(.Object)[["chr"]] <- as.list(chr)
        param(.Object)[["strandLength"]] <- strandLength
        param(.Object)[["FootPrint"]] <- FootPrint
        if(is.null(prefix)){
            param(.Object)[["prefix"]] <- "ATAC_CutSite"
            warning("Please specify a prefix, otherwise your file will be overwrite!")
        }else{
            param(.Object)[["prefix"]] <- prefix
        }

        if(is.null(matrixOutput)){
            output(.Object)[["matrixfile.dir"]] <- 
                getStepWorkDir(.Object,
                                paste0("/Footprint_",
                                       param(.Object)[["prefix"]]))
            output(.Object)[["fp_pdf.dir"]] <- paste0(
                output(.Object)[["matrixfile.dir"]],
                "/FP_PDF_",
                param(.Object)[["prefix"]]
            )
            output(.Object)[["footprint.data"]] <- paste0(
                output(.Object)[["matrixfile.dir"]], "/Footprint_",
                param(.Object)[["prefix"]],
                "_data.rds"
            )
        }else{
            output(.Object)[["matrixfile.dir"]] <- matrixOutput
            output(.Object)[["fp_pdf.dir"]] <- paste0(
                output(.Object)[["matrixfile.dir"]],
                "/FP_PDF_",
                param(.Object)[["prefix"]]
            )
            output(.Object)[["footprint.data"]] <- paste0(
                output(.Object)[["matrixfile.dir"]],
                "/Footprint_",
                param(.Object)[["prefix"]],
                "_data.rds"
            )
        }
        dir.create(output(.Object)[["matrixfile.dir"]])
        dir.create(output(.Object)[["fp_pdf.dir"]])
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "CutSiteCountR",
    definition = function(.Object,...){

        motif_num <- nrow(param(.Object)[["motif_info"]])
        # list to save footprint data
        footprint_data <- list()
        for(i in seq(motif_num)){
            motif_name <- param(.Object)[["motif_info"]][i,1]
            motif_file <- param(.Object)[["motif_info"]][i,2]
            motif_length <- param(.Object)[["motif_info"]][i,3]
            matrixsave.dir <- file.path(output(.Object)[["matrixfile.dir"]], motif_name)
            dir.create(matrixsave.dir)
            footprint.path <- file.path(
                output(.Object)[["fp_pdf.dir"]],
                paste(param(.Object)[["prefix"]], "_", motif_name, ".pdf", sep = "")
            )
            # start!
            writeLog(.Object, sprintf("Generating footprint for %s......", motif_name))
            # .Object <- writeLog(.Object, sprintf("Matrix Destination:%s", matrixsave.dir))
            # .Object <- writeLog(.Object, sprintf("Footprint PDF Destination:%s", output(.Object)[["fp_pdf.dir"]]))
            tmp_dir <- paste(tempdir(), "/", Sys.getpid(), sep="")
            # using tmp dir to save temp data
            dir.create(tmp_dir, FALSE, TRUE, "0700")
            motif_file_index <- .chr_separate_call(ReadsIfile = motif_file,
                                                   ReadsOpath = tmp_dir,
                                                   Name = "/Motif")
            motif_tmp <- paste(tmp_dir, "/Motif", sep = "")
            motif_file_index <- normalizePath(motif_file_index)


            chr <- param(.Object)[["chr"]]
            chr_len <- length(chr)
            data <- data.frame()  # save all matrix
            for(i in seq(1:chr_len)){
                # echo_str <- paste("Now, processing chr", chr[[i]], "......", sep = "")
                # print(echo_str)
                CutSiteInput <- paste0(param(.Object)[["csfile.dir"]], "_chr", chr[[i]], ".cs", collapse = "")
                MotifInput <- normalizePath(
                    paste0(motif_tmp, "_chr", chr[[i]], ".bed", collapse = "")
                )
                if(!file.exists(CutSiteInput)){
                    # echo_str <- paste("There is no cut site in chr", chr[[i]], ", skip!", sep = "")
                    # print(echo_str)
                    next
                }
                if(!(MotifInput %in% motif_file_index)){
                    # echo_str <- paste("There is no motif occurance in chr", chr[[i]], ", skip!", sep = "")
                    # print(echo_str)
                    next
                }
                MatrixOutput <- paste0(matrixsave.dir, "/", motif_name , "_chr", chr[[i]], ".matrix", collapse = "")
                # only two file exist, the program will run
                .CutSiteCount(readsfile = CutSiteInput, motiffile = MotifInput, matrixfile = MatrixOutput,
                              motif_len = motif_length, strand_len = param(.Object)[["strandLength"]])
                temp <- try(read.table(MatrixOutput), silent = TRUE)
                if(inherits(temp, "try-error")){
                    temp <- data.frame()
                }
                data <- rbind(data, temp)

                # echo_str <- paste("Now, finishing chr", chr[[i]], "......", sep = "")
                # print(echo_str)
            }
            if(param(.Object)[["FootPrint"]]){
                if(nrow(data) == 0){
                    tmp_len <- motif_length + 2 * param(.Object)[["strandLength"]]
                    fp <- seq(from = 0, to = 0, length.out = tmp_len)
                }else{
                    fp <- apply(data, 2, sum)
                }
                footprint_data[[motif_name]] <- fp
                pdf(file = footprint.path)
                plot(fp, type = "l", col = "blue", lwd = 2, xlab = "Relative Distance From Motif (bp)", ylab = "Cut Site Count", xaxt = "n", yaxt = "n")
                axis(1, at = seq(1, param(.Object)[["strandLength"]], len = 3),
                     labels = -(param(.Object)[["strandLength"]] + 1 - seq(1, param(.Object)[["strandLength"]] + 1, len = 3)),
                     padj = -1.0, tck = -0.01)
                axis(1, at = param(.Object)[["strandLength"]] + motif_length + seq(1, param(.Object)[["strandLength"]], len = 3),
                     labels = seq(0, param(.Object)[["strandLength"]], len = 3),
                     padj = -1.0, tck = -0.01)
                axis(2, padj = 1.0,tck = -0.02)
                abline(v = c(param(.Object)[["strandLength"]], param(.Object)[["strandLength"]] + motif_length + 1),
                       lty = 2)
                dev.off()
            }

        }
        saveRDS(object = footprint_data, file = output(.Object)[["footprint.data"]])
        
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "CutSiteCountR",
    definition = function(.Object, ...){
        footprint_data <- readRDS( output(.Object)[["footprint.data"]])
        report(.Object)[["footprint.data"]] <- footprint_data
        report(.Object)[["pdf.dir"]] <- output(.Object)[["fp_pdf.dir"]]
        .Object
    }
)



#' @name CutSiteCountR
#' @title Count cut site number in given motif region and plot footprint.

#' @description This function is used to count cut site number in given motif
#' regions and plot footprint. Multi-motif is supported.
#' NOTE: The input parameter is a a little bit complex,
#' \code{atacExtractCutSite} and \code{atacMotifScan} is recommended to use which
#' makes the entire procedure easier.
#' @param atacProcCutSite \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacExtractCutSite}}.
#' @param atacProcMotifScan \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacMotifScan}}.
#' @param csInput Your cut site information file(from atacExtractCutSite
#' function, separated by chromatin name and all cut site are sorted) path
#' with prefix. e.g. "/your_cut_site_information_path/prefix".
#' @param motif_info A rds file from function \code{\link{atacMotifScan}}.
#' In the rds file, it saves 3 column information(motif, motif exact position
#' information file path and motif length).
#' @param chr Which chromatin the program will processing. It must be identical
#' with the filename of cut site information files or subset of .
#' Default:c(1:22, "X", "Y").
#' @param matrixOutput The output directory, where to save your cut site count
#' of every motif position. an empty folder would be great.
#' Default:tmpdir/Footprint
#' @param strandLength How many bp(base pair) do you want to count
#' up/downstream of the motif. default:100.
#' @param FootPrint TRUE or FALSE, plot footprint or not.
#' @param prefix prefix for the pdf file.
#' @param ... Additional arguments, currently unused.
#' @details The parameter is simplified because of too many input file.
#' parameter \code{atacProcCutSite} and \code{atacProcMotifScan} contains all
#' input information so function \code{\link{atacExtractCutSite}} and
#' \code{\link{atacMotifScan}} is recommended to use together. For instance,
#' if you want footprint of 3 TFs (transcription factor) of human in
#' chr1-22, X, Y, then you need 24 chromatin cut site files, 3 motif position
#' files as well as 3 integers of the motif. Function \code{atacExtractCutSite} and
#' \code{atacMotifScan} will do all this, you just specify which motif you want.
#' Therefore, \code{\link{atacExtractCutSite}} and \code{\link{atacMotifScan}} is
#' recommended to use together.
#' @return An invisible \code{\link{ATACProc-class}} object scalar.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' ## processing bed file
#' fra_path <- system.file("extdata", "chr20.50000.bed.bz2", package="esATAC")
#' frag <- as.vector(bunzip2(filename = fra_path,
#' destname = file.path(getwd(), "chr20.50000.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' cs.data <- extractcutsite(bedInput = frag, prefix = "ATAC")
#'
#' ## find motif position
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="esATAC")
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN = bzfile, overwrite=TRUE, remove = FALSE))
#' # motif <- readRDS(system.file("extdata", "MotifPFM.rds", package="esATAC"))
#' # motif.data <- motifscan(peak = peak1_path, genome = BSgenome.Hsapiens.UCSC.hg19, motifs = motif)
#'
#' ## plot footprint
#' # atacCutSiteCount(atacProcCutSite = cs.data, atacProcMotifScan = motif.data)
#'
#'
#' @seealso
#' \code{\link{atacExtractCutSite}}
#' \code{\link{atacMotifScan}}
#'

setGeneric("atacCutSiteCount",
           function(atacProcCutSite, atacProcMotifScan = NULL, csInput = NULL,
                    motif_info = NULL, chr = c(1:22, "X", "Y"), matrixOutput = NULL,
                    strandLength = 100, FootPrint = TRUE, prefix = NULL, ...) standardGeneric("atacCutSiteCount"))


#' @rdname CutSiteCountR
#' @aliases atacCutSiteCount
#' @export
setMethod(
    f = "atacCutSiteCount",
    signature = "ATACProc",
    definition = function(atacProcCutSite, atacProcMotifScan = NULL, csInput = NULL,
                          motif_info = NULL, chr = c(1:22, "X", "Y"), matrixOutput = NULL,
                          strandLength = 100, FootPrint = TRUE, prefix = NULL, ...){

        allpara <- c(list(Class = "CutSiteCountR", prevSteps = list(atacProcCutSite)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname CutSiteCountR
#' @aliases cutsitecount
#' @export
cutsitecount <- function(csInput = NULL, motif_info = NULL, chr = c(1:22, "X", "Y"), matrixOutput = NULL,
                         strandLength = 100, FootPrint = TRUE, prefix = NULL, ...){

    allpara <- c(list(Class = "CutSiteCountR", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
