CutSiteCountR <- R6::R6Class(
    classname = "CutSiteCountR",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProcCutSite, atacProcMotifScan, csInput = NULL,
                              motif_info = NULL, chr = c(1:22, "X", "Y"), matrixOutput = NULL,
                              strandLength = NULL, FootPrint = TRUE, prefix = NULL, editable = FALSE){
            super$initialize("CutSiteCountR", editable, list(arg1 = atacProcCutSite, arg2 = atacProcMotifScan))

            # necessary parameters
            if(!is.null(atacProcCutSite)){
                private$paramlist[["csfile.dir"]] <- atacProcCutSite$getParam("csfile.dir");
            }else{
                private$paramlist[["csfile.dir"]] <- csInput
            }
            if(!is.null(atacProcMotifScan)){
                private$paramlist[["motif_info"]] <- readRDS(atacProcMotifScan$getParam("rdsOutput"))
            }else{
                private$paramlist[["motif_info"]] <- read.table(motif_info)
            }
            private$paramlist[["chr"]] <- as.list(chr)
            private$paramlist[["strandLength"]] <- strandLength
            private$paramlist[["FootPrint"]] <- FootPrint
            private$paramlist[["prefix"]] <- prefix

            if(is.null(matrixOutput)){
                private$paramlist[["matrixfile.dir"]] <- paste(
                    .obtainConfigure("tmpdir"),
                    "/Footprint_",
                    private$paramlist[["prefix"]],
                    sep = ""
                )
                private$paramlist[["footprint.data"]] <- paste(
                    private$paramlist[["matrixfile.dir"]],
                    "_data.rds",
                    sep = ""
                )
            }else{
                private$paramlist[["matrixfile.dir"]] <- matrixOutput
                private$paramlist[["footprint.data"]] <- paste(
                    private$paramlist[["matrixfile.dir"]],
                    "/Footprint_",
                    private$paramlist[["prefix"]],
                    "_data.rds",
                    sep = ""
                )
            }
            dir.create(private$paramlist[["matrixfile.dir"]])


            # parameter check
            private$paramValidation()
        } # initialization end

    ), # public end

    private = list(
        processing = function(){
            private$writeLog(paste0("Now, start processing!"))
            motif_num <- nrow(private$paramlist[["motif_info"]])
            # list to save footprint data
            footprint_data <- list()
            for(i in seq(motif_num)){
                motif_name <- private$paramlist[["motif_info"]][i,1]
                motif_file <- private$paramlist[["motif_info"]][i,2]
                motif_length <- private$paramlist[["motif_info"]][i,3]
                matrixsave.dir <- file.path(private$paramlist[["matrixfile.dir"]], motif_name)
                dir.create(matrixsave.dir)
                footprint.path <- file.path(
                    .obtainConfigure("tmpdir"),
                    paste(private$paramlist[["prefix"]], "_", motif_name, ".pdf", sep = "")
                )
                # start!
                private$writeLog(sprintf("Start Processing %s", motif_name))
                private$writeLog(sprintf("Matrix Destination:%s", matrixsave.dir))
                private$writeLog(sprintf("Footprint PDF Destination:%s", matrixsave.dir))
                tmp_dir <- paste(tempdir(), "/", Sys.getpid(), sep="")
                # using tmp dir to save temp data
                dir.create(tmp_dir, FALSE, TRUE, "0700")
                motif_file_index <- .chr_separate_call(ReadsIfile = motif_file,
                                                       ReadsOpath = tmp_dir,
                                                       Name = "/Motif")
                motif_tmp <- paste(tmp_dir, "/Motif", sep = "")
                motif_file_index <- normalizePath(motif_file_index)

                chr <- private$paramlist[["chr"]]
                chr_len <- length(chr)
                for(i in seq(1:chr_len)){
                    echo_str <- paste("Now, processing chr", chr[[i]], "......", sep = "")
                    print(echo_str)
                    CutSiteInput <- paste0(private$paramlist[["csfile.dir"]], "_chr", chr[[i]], ".cs", collapse = "")
                    MotifInput <- normalizePath(
                        paste0(motif_tmp, "_chr", chr[[i]], ".bed", collapse = "")
                    )
                    if(!file.exists(CutSiteInput)){
                        echo_str <- paste("There is no cut site in chr", chr[[i]], ", skip!", sep = "")
                        print(echo_str)
                        next
                    }
                    if(!(MotifInput %in% motif_file_index)){
                        echo_str <- paste("There is no motif occurance in chr", chr[[i]], ", skip!", sep = "")
                        print(echo_str)
                        next
                    }
                    MatrixOutput <- paste0(matrixsave.dir, "/", motif_name , "_chr", chr[[i]], ".matrix", collapse = "")
                    .CutSiteCount(readsfile = CutSiteInput, motiffile = MotifInput, matrixfile = MatrixOutput,
                                  motif_len = motif_length, strand_len = private$paramlist[["strandLength"]])
                    if(i == 1){
                        data <- try(read.table(MatrixOutput), silent = TRUE)
                        if(inherits(data, "try-error")){
                            data <- data.frame()
                        }
                    }else{
                        temp <- try(read.table(MatrixOutput), silent = TRUE)
                        if(inherits(temp, "try-error")){
                            temp <- data.frame()
                        }
                        data <- rbind(data, temp)
                    }
                }
                if(private$paramlist[["FootPrint"]]){
                    fp <- apply(data, 2, sum)
                    footprint_data[[motif_name]] <- fp
                    pdf(file = footprint.path)
                    plot(fp, type = "l", col = "blue", lwd = 2, xlab = "Relative Distance From Motif (bp)", ylab = "Cut Site Count", xaxt = "n", yaxt = "n")
                    axis(1, at = seq(1, private$paramlist[["strandLength"]], len = 3),
                         labels = -(private$paramlist[["strandLength"]] + 1 - seq(1, private$paramlist[["strandLength"]] + 1, len = 3)),
                         padj = -1.0, tck = -0.01)
                    axis(1, at = private$paramlist[["strandLength"]] + motif_length + seq(1, private$paramlist[["strandLength"]], len = 3),
                         labels = seq(0, private$paramlist[["strandLength"]], len = 3),
                         padj = -1.0, tck = -0.01)
                    axis(2, padj = 1.0,tck = -0.02)
                    abline(v = c(private$paramlist[["strandLength"]], private$paramlist[["strandLength"]] + motif_length + 1),
                           lty = 2)
                    dev.off()
                }

            }
            saveRDS(object = footprint_data, file = private$paramlist[["footprint.data"]])

        }, # processing end

        checkRequireParam = function(){
            if(is.null(private$paramlist[["csfile.dir"]])){
                stop("Parameter csInput is required!")
            }
        }, # checkRequireParam end

        checkAllPath = function(){
            private$checkPathExist(private$paramlist[["csfile.dir"]])
        }, # checkAllPath end

        getReportValImp = function(item){
            if(item == "footprint.data"){
                fp <- readRDS(private$paramlist[["footprint.data"]])
                return(fp)
            }else if(item == "pdf.dir"){
                return(.obtainConfigure("tmpdir"))
            }
        },

        getReportItemsImp = function(){
            return(c("footprint.data", "pdf.dir"))
        }
    ) # private end

) # class end


#' @name atacCutSiteCount
#' @aliases atacCutSiteCount
#' @aliases cutsitecount
#' @title Count cut site number in given motif region.
#' @description This function is used to count cut site number in given motif
#' region and plot footprint. Multi-motif is supported.
#' NOTE: The input parameter is a a little bit complex,
#' \code{atacCutSitePre} and \code{atacMotifScan} is recommended to use which
#' makes the entire procedure easier.
#' @param atacProcCutSite \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacCutSitePre}}.
#' @param atacProcMotifScan \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacMotifScan}}.
#' @param csInput Your cut site information file(from atacCutSitePre function,
#' separated by chromatin name and all cut site are sorted) path with prefix.
#' e.g. "/your_cut_site_information_path/prefix"
#' @param motif_info A rds file from function \code{\link{atacMotifScan}}.
#' In the rds file, it saves 3 column information(motif, motif exact position
#' information file path and motif length).
#' @param chr Which chromatin the program will processing. It must be identical
#' with the file name of cut site information file.
#' Default:c(1:22, "X", "Y").
#' @param matrixOutput The output directory, where to save your cut site count
#' of every motif position. an empty folder would be great.
#' Default:tmpdir/Footprint
#' @param strandLength How many bp(base pair) do you want to count
#' up/downstream of the motif. default:100.
#' @param FootPrint TRUE or FALSE, plot footprint or not.
#' @param prefix prefix for the pdf file.
#' @details The parameter is simplified because of too many input file.
#' parameter \code{atacProcCutSite} and \code{atacProcMotifScan} contains all
#' input information so function \code{\link{atacCutSitePre}} and
#' \code{\link{atacMotifScan}} is recommended to use together. For instance,
#' if you want footprint of 3 TFs (transcription factor) of human in
#' chr1-22, X, Y, then you need 24 chromatin cut site files, 3 motif position
#' files as well as 3 integers of the motif. Function \code{atacCutSitePre} and
#' \code{atacMotifScan} will do all this, you just specify which motif you want.
#' Therefore, \code{\link{atacCutSitePre}} and \code{\link{atacMotifScan}} is
#' recommended to use together.
#' @return An invisible \code{\link{ATACProc}} object scalar.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' ## processing bed file
#' fra_path <- system.file("extdata", "chr20.50000.bed.bz2", package="ATACFlow")
#' frag <- as.vector(bunzip2(filename = fra_path,
#' destname = file.path(getwd(), "chr20.50000.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' cs.data <- extractcutsite(bedInput = frag, prefix = "ATAC")
#'
#' ## find motif position
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="ATACFlow")
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN = bzfile, overwrite=TRUE, remove = FALSE))
#' pwm <- readRDS(system.file("extdata", "motifPWM.rds", package="ATACFlow"))
#' motif.data <- motifscan(peak = peak1_path, genome = BSgenome.Hsapiens.UCSC.hg19,
#' motifPWM = pwm)
#'
#' ## plot footprint
#' atacCutSiteCount(atacProcCutSite = cs.data, atacProcMotifScan = motif.data)
#'
#'
#' @seealso
#' \code{\link{atacCutSitePre}}
#' \code{\link{atacMotifScan}}
#'

#' @rdname atacCutSiteCount
#' @export
atacCutSiteCount <- function(atacProcCutSite = NULL, atacProcMotifScan = NULL, csInput = NULL,
                             motif_info = NULL, chr = c(1:22, "X", "Y"), matrixOutput = NULL,
                             strandLength = 100, FootPrint = TRUE, prefix = "Motif"){
    tmp <- CutSiteCountR$new(atacProcCutSite, atacProcMotifScan, csInput,
                             motif_info, chr, matrixOutput, strandLength, FootPrint, prefix)
    tmp$process()
    invisible(tmp)
}

#' @rdname atacCutSiteCount
#' @export
cutsitecount <- function(csInput, motif_info, chr = c(1:22, "X", "Y"),
                         matrixOutput = NULL, strandLength = 100,
                         FootPrint = TRUE, prefix = "Motif"){
    tmp <- CutSiteCountR$new(atacProcCutSite = NULL, atacProcMotifScan = NULL, csInput,
                             motif_info, chr, matrixOutput, strandLength, FootPrint, prefix)
    tmp$process()
    invisible(tmp)
}
