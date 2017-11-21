setClass(Class = "RMotifScan",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "RMotifScan",
    definition = function(.Object, atacProc, ..., peak = NULL, genome = NULL,
                          motifPWM = NULL, min.score = NULL,
                          scanO.dir = NULL, n.cores = NULL,
                          prefix = NULL, editable = FALSE){
        .Object <- init(.Object, "RMotifScan", editable, list(arg1 = atacProc))

        # necessary parameters
        if(!is.null(atacProc)){
            .Object@paramlist[["peak"]] <- getParam(atacProc, "bedOutput");
        }else{
            .Object@paramlist[["peak"]] <- peak
        }
        if(!is.null(genome)){
            .Object@paramlist[["genome"]] <- genome
        }else{
            .Object@paramlist[["genome"]] <- .obtainConfigure("bsgenome")
        }
        .Object@paramlist[["motifPWM"]] <- motifPWM
        .Object@paramlist[["motifPWM.len"]] <- lapply(X = .Object@paramlist[["motifPWM"]], FUN = ncol)
        .Object@paramlist[["min.score"]] <- min.score
        if(is.null(prefix)){
            .Object@paramlist[["prefix"]] <- "motifscan"
        }else{
            .Object@paramlist[["prefix"]] <- prefix
        }
        # unnecessary parameters
        if(is.null(scanO.dir)){
            .Object@paramlist[["scanO.dir"]] <- paste(tools::file_path_sans_ext(.Object@paramlist[["peak"]]),
                                                      "_",
                                                      .Object@paramlist[["prefix"]],
                                                      "_MotifScanOutput",
                                                      sep = "")
            dir.create(.Object@paramlist[["scanO.dir"]])
        }else{
            .Object@paramlist[["scanO.dir"]] <- scanO.dir
        }
        .Object@paramlist[["rdsOutput"]] <- paste(
            .Object@paramlist[["scanO.dir"]],
            "/", .Object@paramlist[["prefix"]], "_",
            "RMotifScan.rds",
            sep = ""
        )

        if(is.null(n.cores)){
            .Object@paramlist[["n.cores"]] <- .obtainConfigure("threads")
        }else{
            .Object@paramlist[["n.cores"]] <- n.cores
        }

        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "RMotifScan",
    definition = function(.Object,...){
        .Object <- writeLog(.Object, paste0("processing file:"))
        .Object <- writeLog(.Object, sprintf("peak file:%s", .Object@paramlist[["peak"]]))
        .Object <- writeLog(.Object, sprintf("Output destination:%s", .Object@paramlist[["scanO.dir"]]))

        # running
        peak <- rtracklayer::import(.Object@paramlist[["peak"]])
        save_info <- data.frame()

        # processing 2*n.core motifs in each turn
        k <- .Object@paramlist[["n.cores"]] * 2
        n_motif <- length(.Object@paramlist[["motifPWM"]])
        motif_in_group <- split(.Object@paramlist[["motifPWM"]],
                                rep(1:ceiling(n_motif/k), each=k)[1:n_motif])
        n_group <- length(motif_in_group)

        # write order(motif index) while writing save_info
        WriteMotifOrder <- 1
        cl <- parallel::makeCluster(.Object@paramlist[["n.cores"]])
        for(i in seq(n_group)){
            thisGroup.motif <- motif_in_group[[i]]
            thisGroup.motifname <- names(thisGroup.motif)
            thisGroup.motifnum <- length(thisGroup.motif)
            thisGroup.motifinfo <- paste("Now, processing the following motif: ",
                                         paste(thisGroup.motifname, collapse = ","),
                                         sep = "")
            print(thisGroup.motifinfo)
            sitesetList <- parLapply(cl = cl,
                                     X = thisGroup.motif,
                                     fun = Biostrings::matchPWM,
                                     subject = .Object@paramlist[["genome"]],
                                     min.score = .Object@paramlist[["min.score"]],
                                     with.score = TRUE)

            for(i in seq(thisGroup.motifnum)){
                motif_name <- names(sitesetList[i])
                output_data <- IRanges::subsetByOverlaps(x = sitesetList[[i]],
                                                         ranges = peak,
                                                         ignore.strand = TRUE)
                output_data <- sort(x = output_data, ignore.strand = TRUE)
                output_data <- as.data.frame(output_data)
                output_data <- within(output_data, rm(width))
                output_path <- paste(.Object@paramlist[["scanO.dir"]],
                                     "/", .Object@paramlist[["prefix"]], "_",
                                     motif_name, sep = "")
                motif_len <- .Object@paramlist[["motifPWM.len"]][[motif_name]]
                save_info[WriteMotifOrder, 1] <- motif_name
                save_info[WriteMotifOrder, 2] <- R.utils::getAbsolutePath(output_path)
                save_info[WriteMotifOrder, 3] <- motif_len
                WriteMotifOrder <- WriteMotifOrder + 1
                write.table(x = output_data, file = output_path, row.names = FALSE,
                            col.names = FALSE, quote = FALSE)
            }
        }
        stopCluster(cl)

        saveRDS(object = save_info, file = .Object@paramlist[["rdsOutput"]])
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "RMotifScan",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["peak"]])){
            stop("Parameter peak is required!")
        }
        if(is.null(.Object@paramlist[["motifPWM"]])){
            stop("Parameter motifPWM is required!")
            if(!is.list(.Object@paramlist[["motifPWM"]])){
                stop("Parameter motifPWM must be a list!")
            }
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "RMotifScan",
    definition = function(.Object,...){
        checkFileExist(.Object, .Object@paramlist[["peak"]]);
        checkPathExist(.Object, .Object@paramlist[["scanO.dir"]]);
    }
)


#' @name RMotifScan
#' @title Search Motif Position in Given Regions
#' @description
#' Search motif position in given genome regions according PWM matrix.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}}.
#' @param peak \code{Character} scalar.
#' Input region path. UCSC bed file is recommented. Other file should be able
#' to import as \link[GenomicRanges]{GRanges} objects through
#' \link[rtracklayer]{import}.
#' @param genome A DNAString object.
#' @param motifPWM \code{list} scalar. Default: from \code{\link{setConfigure}}.
#' Every element in the \code{list} contains a motif PWM matrix.
#' e.g. pwm <- list("CTCF" = CTCF_PWMmatrix)
#' @param min.score The minimum score for counting a match. Can be given as a
#' character string containing a percentage (e.g. "85%") of the highest
#' possible score or as a single number.
#' @param scanO.dir \code{Character} scalar.
#' the output file directory. This function will use the index in motifPWM as
#' the file name to save the motif position information in separate files.
#' @param n.cores How many core to run this function.
#' Default: from \code{\link{setConfigure}}.
#' @param prefix prefix for Output file.
#' @param ... Additional arguments, currently unused.
#' @details This function scan motif position in a given genome regions.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(R.utils)
#' p1bz <- system.file("extdata", "Example_peak1.bed.bz2", package="esATAC")
#' peak1_path <- as.vector(bunzip2(filename = p1bz,
#' destname = file.path(getwd(), "Example_peak1.bed"),
#' ext="bz2", FUN = bzfile, overwrite=TRUE, remove = FALSE))
#' pwm <- readRDS(system.file("extdata", "motifPWM.rds", package="esATAC"))
#' #motifscan(peak = peak1_path, genome = BSgenome.Hsapiens.UCSC.hg19,
#' #motifPWM = pwm, prefix = "test")
#'
#'
#' @seealso
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacCutSiteCount}}
#' \link[Biostrings]{matchPWM}
#' \link[IRanges]{subsetByOverlaps}
#'
#' @importFrom rtracklayer import
#' @importFrom IRanges subsetByOverlaps
#' @importFrom R.utils getAbsolutePath
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#'

setGeneric("atacMotifScan",
           function(atacProc, peak = NULL, genome = NULL,
                    motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                    n.cores = NULL, prefix = NULL, ...) standardGeneric("atacMotifScan"))



#' @rdname RMotifScan
#' @aliases atacMotifScan
#' @export
setMethod(
    f = "atacMotifScan",
    signature = "ATACProc",
    function(atacProc, peak = NULL, genome = NULL,
             motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
             n.cores = NULL, prefix = NULL, ...){
        atacproc <- new(
            "RMotifScan",
            atacProc = atacProc,
            peak = peak,
            genome = genome,
            motifPWM = motifPWM,
            min.score = min.score,
            scanO.dir = scanO.dir,
            n.cores = n.cores,
            prefix = prefix)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname RMotifScan
#' @aliases motifscan
#' @export
motifscan <- function(peak = NULL, genome = NULL,
                      motifPWM = NULL, min.score = "85%", scanO.dir = NULL,
                      n.cores = NULL, prefix = NULL, ...){
    atacproc <- new(
        "RMotifScan",
        atacProc = NULL,
        peak = peak,
        genome = genome,
        motifPWM = motifPWM,
        min.score = min.score,
        scanO.dir = scanO.dir,
        n.cores = n.cores,
        prefix = prefix)
    atacproc <- process(atacproc)
    invisible(atacproc)
}



#' @name getMotifPWM
#' @title Processing PFM or PWM file.
#' @description
#' atacMotifScan and atacMotifScanPair accept PWM in a \code{list}, this
#' function convert a PFM or PWM file(in JASPAR format) to a list in R.
#' @param motif.file PFM or PWM file.
#' @param is.PWM TRUE or FALSE. If TRUE, the input file contains PWM, do not
#' need convert to PWM. If FALSE, the input file contains PFM, need convert
#' to PWM. Default:FALSE.
#' @param JASPARdb TRUE or FALSE. Whether use JASPAR database or not.
#' @param Species Taxonomy ID. For human, it's 9606.
#' @param Name The name of the transcription factor.
#' @param ID The ID of the transcription factor.
#' @details Converting a PFM or PWM file(in JASPAR format) to a list in R.
#' @return A list contains PWM.
#' @author Wei Zhang
#' @importFrom TFBSTools PFMatrix
#' @importFrom TFBSTools as.matrix
#' @examples
#'
#' # from files(user customized)
#' pfm_file <- system.file("extdata", "motif.txt", package="esATAC")
#' pwm_list <- getMotifPWM(motif.file = pfm_file, is.PWM = FALSE)
#'
#' # from JASPAR database
#' pwm_list <- getMotifPWM(JASPARdb = TRUE, Name = "TFAP2A")
#'
#' @export
getMotifPWM <- function(motif.file = NULL, is.PWM = FALSE, JASPARdb = FALSE,
                        Species = NULL, Name = NULL, ID = NULL){
    if(!is.null(motif.file)){
        PWMList <- list()
        name.flag <- FALSE
        A.flag <- FALSE
        C.flag <- FALSE
        G.flag <- FALSE
        T.flag <- FALSE

        con <- file(motif.file, "r")
        line <- readLines(con, n = 1)
        while( length(line) != 0 ){

            if(substring(line, 1, 1) == ">"){
                name_str <- unlist(strsplit(x = sub(pattern = ">", replacement = "", x = line), split = "\\s+"))
                motif_name <- tail(name_str, n = 1)
                name.flag <- TRUE
            }else if(substring(line, 1, 1) == "A"){
                A_str_num <- regmatches(line, gregexpr("[[:digit:]]+", line))
                A_num <- as.numeric(unlist(A_str_num))
                A.flag <- TRUE
            }else if(substring(line, 1, 1) == "C"){
                C_str_num <- regmatches(line, gregexpr("[[:digit:]]+", line))
                C_num <- as.numeric(unlist(C_str_num))
                C.flag <- TRUE
            }else if(substring(line, 1, 1) == "G"){
                G_str_num <- regmatches(line, gregexpr("[[:digit:]]+", line))
                G_num <- as.numeric(unlist(G_str_num))
                G.flag <- TRUE
            }else if(substring(line, 1, 1) == "T"){
                T_str_num <- regmatches(line, gregexpr("[[:digit:]]+", line))
                T_num <- as.numeric(unlist(T_str_num))
                T.flag <- TRUE
            }

            if(name.flag & A.flag & C.flag & G.flag & T.flag){
                p_matrix <- matrix(data = c(A_num, C_num, G_num, T_num), nrow = 4,
                                   byrow = TRUE,  dimnames=list(c("A", "C", "G", "T")))

                if(!is.PWM){
                    p_matrix <- TFBSTools::PFMatrix(profileMatrix = p_matrix)
                    p_matrix <- TFBSTools::as.matrix(TFBSTools::toPWM(p_matrix))
                }

                PWMList[[motif_name]] <- p_matrix

                name.flag <- FALSE
                A.flag <- FALSE
                C.flag <- FALSE
                G.flag <- FALSE
                T.flag <- FALSE
            }

            line <- readLines(con, n = 1)
        }
        close(con)
        return(PWMList)
    }else if(!is.null(JASPARdb)){
        print("Now, JASPAR2016 database is in use!")
        if(!is.null(Species)){
            opts <- list()
            opts[["species"]] <- Species
            pwm <- TFBSTools::getMatrixSet(x = JASPAR2016::JASPAR2016, opts = opts)
            pwm <- TFBSTools::toPWM(pwm)
            names(pwm) <- TFBSTools::name(pwm)
            pwm <- lapply(X = pwm, FUN = TFBSTools::as.matrix)
            names(pwm) <- gsub(pattern = "[^a-zA-Z0-9]", replacement = "", x = names(pwm), perl = TRUE)
            return(pwm)
        }else if(!is.null(Name)){
            pwm <- TFBSTools::getMatrixByName(x = JASPAR2016::JASPAR2016, name = Name)
            pwm <- TFBSTools::toPWM(pwm)
            pwm.name <- gsub(pattern = "[^a-zA-Z0-9]", replacement = "", x = TFBSTools::name(pwm), perl = TRUE)
            PWMList <- list(pwm.name = TFBSTools::as.matrix(pwm))
            names(PWMList) <- pwm.name
            return(PWMList)
        }else if(!is.null(ID)){
            pwm <- TFBSTools::getMatrixByID(x = JASPAR2016::JASPAR2016, ID = ID)
            pwm <- TFBSTools::toPWM(pwm)
            pwm.name <- gsub(pattern = "[^a-zA-Z0-9]", replacement = "", x = TFBSTools::name(pwm), perl = TRUE)
            PWMList <- list(pwm.name = TFBSTools::as.matrix(pwm))
            names(PWMList) <- pwm.name
            return(PWMList)
        }
    }

}

