#' @useDynLib esATAC
NULL

#' isRemote
#'
#' Check if path is remote
#'
#' @param x path/s to check
#' @return path
#'
isRemote <- function(x) {
    return(grepl(pattern = "^http|^ftp", x = x))
}


#' Set a default value if an object is null
#'
#' @param x An object to set if it's null
#' @param y The value to provide if x is null
#'
#' @return Returns y if x is null, otherwise returns x.
#'
SetIfNull <- function(x, y) {
    if (is.null(x = x)) {
        return(y)
    } else {
        return(x)
    }
}


#' String to GRanges
#'
#' Convert a genomic coordinate string to a GRanges object
#'
#' @param regions Vector of genomic region strings
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @param ... Additional arguments passed to
#' \code{\link[GenomicRanges]{makeGRangesFromDataFrame}}
#'
#' @return Returns a GRanges object
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tidyr separate
#'
StringToGRanges <- function(regions, sep = c("-", "-"), ...) {
    ranges.df <- data.frame(ranges = regions)
    ranges.df <- separate(
        data = ranges.df,
        col = "ranges",
        sep = paste0(sep[[1]], "|", sep[[2]]),
        into = c("chr", "start", "end")
    )
    granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
    return(granges)
}


#' Run FeatureMatrix on a single Fragment object
#'
#' @param fragment A list of \code{\link{Fragment}} objects. Note that if
#' setting the \code{cells} parameter, the requested cells should be present in
#' the supplied \code{Fragment} objects. However, if the cells information in
#' the fragment object is not set (\code{Cells(fragments)} is \code{NULL}), then
#' the fragment object will still be searched.
#'
#' @param features A GRanges object containing a set of genomic intervals.
#' These will form the rows of the matrix, with each entry recording the number
#' of unique reads falling in the genomic region for each cell.
#'
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#'
#' @param process_n Number of regions to load into memory at a time, per thread.
#' Processing more regions at once can be faster but uses more memory.
#'
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#'
#' @param verbose Display messages
#'
#' @return SingleFeatureMatrix
#'
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#' @importFrom Matrix sparseMatrix
#' @importMethodsFrom GenomicRanges intersect
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom fastmatch fmatch
#'
SingleFeatureMatrix <- function(
        fragment,
        features,
        cells = NULL,
        process_n = 2000,
        sep = c("-", "-"),
        verbose = TRUE
) {
    fragment.path <- GetFragmentData(object = fragment, slot = "path")
    if (!is.null(cells)) {
        # only look for cells that are in the fragment file
        frag.cells <- GetFragmentData(object = fragment, slot = "cells")
        if (is.null(x = frag.cells)) {
            # cells information not set in fragment object
            names(x = cells) <- cells
        } else {
            # first subset frag.cells
            cell.idx <- fmatch(
                x = names(x = frag.cells),
                table = cells,
                nomatch = 0L
            ) > 0
            cells <- frag.cells[cell.idx]
        }
    }
    tbx <- TabixFile(file = fragment.path)
    features <- keepSeqlevels(
        x = features,
        value = intersect(
            x = seqnames(x = features),
            y = seqnamesTabix(file = tbx)
        ),
        pruning.mode = "coarse"
    )
    if (length(x = features) == 0) {
        stop("No matching chromosomes found in fragment file.")
    }
    
    feature.list <- ChunkGRanges(
        granges = features,
        nchunk = ceiling(x = length(x = features) / process_n)
    )
    if (verbose) {
        message("Extracting reads overlapping genomic regions")
    }
    if (nbrOfWorkers() > 1) {
        matrix.parts <- future_lapply(
            X = feature.list,
            FUN = PartialMatrix,
            tabix = tbx,
            cells = cells,
            sep = sep,
            future.globals = list(),
            future.scheduling = FALSE
        )
    } else {
        mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
        matrix.parts <- mylapply(
            X = feature.list,
            FUN = PartialMatrix,
            tabix = tbx,
            cells = cells,
            sep = sep
        )
    }
    # remove any that are NULL (no fragments for any cells in the region)
    null.parts <- sapply(X = matrix.parts, FUN = is.null)
    matrix.parts <- matrix.parts[!null.parts]
    if (is.null(x = cells)) {
        all.cells <- unique(
            x = unlist(x = lapply(X = matrix.parts, FUN = colnames))
        )
        matrix.parts <- lapply(
            X = matrix.parts,
            FUN = AddMissingCells,
            cells = all.cells
        )
    }
    featmat <- do.call(what = rbind, args = matrix.parts)
    if (!is.null(x = cells)) {
        # cells supplied, rename with cell name from object rather than file
        cell.convert <- names(x = cells)
        names(x = cell.convert) <- cells
        colnames(x = featmat) <- unname(obj = cell.convert[colnames(x = featmat)])
    }
    # reorder features
    feat.str <- GRangesToString(grange = features, sep = sep)
    featmat <- featmat[feat.str, , drop=FALSE]
    return(featmat)
}



#' Split a genomic ranges object into evenly sized chunks
#'
#' @param granges A GRanges object
#' @param nchunk Number of chunks to split into
#'
#' @return Returns a list of GRanges objects
ChunkGRanges <- function(granges, nchunk) {
    if (length(x = granges) < nchunk) {
        nchunk <- length(x = granges)
    }
    chunksize <- as.integer(x = (length(granges) / nchunk))
    range.list <- sapply(X = seq_len(length.out = nchunk), FUN = function(x) {
        chunkupper <- (x * chunksize)
        if (x == 1) {
            chunklower <- 1
        } else {
            chunklower <- ((x - 1) * chunksize) + 1
        }
        if (x == nchunk) {
            chunkupper <- length(x = granges)
        }
        return(granges[chunklower:chunkupper])
    })
    return(range.list)
}


#' GRanges to String
#'
#' Convert GRanges object to a vector of strings
#'
#' @param grange A GRanges object
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @return Returns a character vector
#' @export
#' @concept utilities
#' @importMethodsFrom GenomicRanges start end seqnames
#' @examples
#' print("see https://satijalab.org/signac/reference/grangestostring")
GRangesToString <- function(grange, sep = c("-", "-")) {
    regions <- paste0(
        as.character(x = seqnames(x = grange)),
        sep[[1]],
        start(x = grange),
        sep[[2]],
        end(x = grange)
    )
    return(regions)
}


#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors elementNROWS
PartialMatrix <- function(tabix, regions, sep = c("-", "-"), cells = NULL) {
    # construct sparse matrix for one set of regions
    # names of the cells vector can be ignored here, conversion is handled in
    # the parent functions
    open(con = tabix)
    cells.in.regions <- GetCellsInRegion(
        tabix = tabix,
        region = regions,
        cells = cells
    )
    close(con = tabix)
    gc(verbose = FALSE)
    nrep <- elementNROWS(x = cells.in.regions)
    if (all(nrep == 0) & !is.null(x = cells)) {
        # no fragments
        # zero for all requested cells
        featmat <- sparseMatrix(
            dims = c(length(x = regions), length(x = cells)),
            i = NULL,
            j = NULL
        )
        rownames(x = featmat) <- GRangesToString(grange = regions)
        colnames(x = featmat) <- cells
        featmat <- as(object = featmat, Class = "dgCMatrix")
        return(featmat)
    } else if (all(nrep == 0)) {
        # no fragments, no cells requested
        # create empty matrix
        featmat <- sparseMatrix(
            dims = c(length(x = regions), 0),
            i = NULL,
            j = NULL
        )
        rownames(x = featmat) <- GRangesToString(grange = regions)
        featmat <- as(object = featmat, Class = "dgCMatrix")
        return(featmat)
    } else {
        # fragments detected
        if (is.null(x = cells)) {
            all.cells <- unique(x = unlist(x = cells.in.regions))
            cell.lookup <- seq_along(along.with = all.cells)
            names(x = cell.lookup) <- all.cells
        } else {
            cell.lookup <- seq_along(along.with = cells)
            names(cell.lookup) <- cells
        }
        # convert cell name to integer
        cells.in.regions <- unlist(x = cells.in.regions)
        cells.in.regions <- unname(obj = cell.lookup[cells.in.regions])
        all.features <- GRangesToString(grange = regions, sep = sep)
        feature.vec <- rep(x = seq_along(along.with = all.features), nrep)
        featmat <- sparseMatrix(
            i = feature.vec,
            j = cells.in.regions,
            x = rep(x = 1, length(x = cells.in.regions))
        )
        featmat <- as(Class = "dgCMatrix", object = featmat)
        rownames(x = featmat) <- all.features[seq(max(feature.vec))]
        colnames(x = featmat) <- names(x = cell.lookup)[seq(max(cells.in.regions))]
        # add zero columns for missing cells
        if (!is.null(x = cells)) {
            featmat <- AddMissingCells(x = featmat, cells = cells)
        }
        # add zero rows for missing features
        missing.features <- all.features[!(all.features %in% rownames(x = featmat))]
        if (length(x = missing.features) > 0) {
            null.mat <- sparseMatrix(
                i = c(),
                j = c(),
                dims = c(length(x = missing.features), ncol(x = featmat))
            )
            rownames(x = null.mat) <- missing.features
            null.mat <- as(object = null.mat, Class = "dgCMatrix")
            featmat <- rbind(featmat, null.mat)
        }
        return(featmat)
    }
}

#' Get cells in a region
#'
#' Extract cell names containing reads mapped within a given genomic region
#'
#' @param tabix Tabix object
#' @param region A string giving the region to extract from the fragments file
#' @param cells Vector of cells to include in output. If NULL, include all cells
#'
#' @importFrom Rsamtools scanTabix
#' @importFrom methods is
#' @importFrom fastmatch fmatch
#' @export
#' @concept utilities
#' @return Returns a list
GetCellsInRegion <- function(tabix, region, cells = NULL) {
    if (!is(object = region, class2 = "GRanges")) {
        region <- StringToGRanges(regions = region)
    }
    reads <- scanTabix(file = tabix, param = region)
    reads <- lapply(X = reads, FUN = ExtractCell)
    if (!is.null(x = cells)) {
        reads <- sapply(X = reads, FUN = function(x) {
            x <- x[fmatch(x = x, table = cells, nomatch = 0L) > 0L]
            if (length(x = x) == 0) {
                return(NULL)
            } else {
                return(x)
            }
        }, simplify = FALSE)
    }
    return(reads)
}


#' @importFrom Matrix sparseMatrix
AddMissingCells <- function(x, cells) {
    # add columns with zeros for cells not in matrix
    missing.cells <- setdiff(x = cells, y = colnames(x = x))
    if (!(length(x = missing.cells) == 0)) {
        null.mat <- sparseMatrix(
            i = c(),
            j = c(),
            dims = c(nrow(x = x), length(x = missing.cells))
        )
        rownames(x = null.mat) <- rownames(x = x)
        colnames(x = null.mat) <- missing.cells
        x <- cbind(x, null.mat)
    }
    x <- x[, cells, drop = FALSE]
    return(x)
}


#' Extract cell
#'
#' Extract cell barcode from list of tab delimited character
#' vectors (output of \code{\link{scanTabix}})
#'
#' @param x List of character vectors
#' @return Returns a string
#' @importFrom stringi stri_split_fixed
ExtractCell <- function(x) {
    if (length(x = x) == 0) {
        return(NULL)
    } else {
        x <- stri_split_fixed(str = x, pattern = "\t")
        n <- length(x = x)
        x <- unlist(x = x)
        return(unlist(x = x)[5 * (seq(n)) - 1])
    }
}


#' Extract genomic ranges from EnsDb object
#'
#' Pulls the transcript information for all chromosomes from an EnsDb object.
#' This wraps \code{\link[biovizBase]{crunch}} and applies the extractor
#' function to all chromosomes present in the EnsDb object.
#'
#' @param ensdb An EnsDb object
#' @param standard.chromosomes Keep only standard chromosomes
#' @param biotypes Biotypes to keep
#' @param verbose Display messages
#' @return GRanges
#'
#' @importFrom GenomeInfoDb keepStandardChromosomes seqinfo
#' @importFrom biovizBase crunch
#' @concept utilities
#' @export
#' @examples
#' print("see https://satijalab.org/signac/reference/getgrangesfromensdb")
GetGRangesFromEnsDb <- function(
        ensdb,
        standard.chromosomes = TRUE,
        biotypes = c("protein_coding", "lincRNA", "rRNA", "processed_transcript"),
        verbose = TRUE
) {
    # convert seqinfo to granges
    whole.genome <-  as(object = seqinfo(x = ensdb), Class = "GRanges")
    if (standard.chromosomes) {
        whole.genome <- keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
    }
    
    # extract genes from each chromosome
    if (verbose) {
        tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
            biovizBase::crunch(
                obj = ensdb,
                which = whole.genome[x],
                columns = c("tx_id", "gene_name", "gene_id", "gene_biotype"))
        })
    } else {
        tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
            suppressMessages(expr = biovizBase::crunch(
                obj = ensdb,
                which = whole.genome[x],
                columns = c("tx_id", "gene_name", "gene_id", "gene_biotype")))
        })
    }
    
    # combine
    tx <- do.call(what = c, args = tx)
    tx <- tx[tx$gene_biotype %in% biotypes]
    return(tx)
}


#' Construct a feature x cell matrix from a genomic fragments file
#'
#' @param fragments A list of \code{\link{Fragment}} objects. Note that if
#' setting the \code{cells} parameter, the requested cells should be present in
#' the supplied \code{Fragment} objects. However, if the cells information in
#' the fragment object is not set (\code{Cells(fragments)} is \code{NULL}), then
#' the fragment object will still be searched.
#' @param features A GRanges object containing a set of genomic intervals.
#' These will form the rows of the matrix, with each entry recording the number
#' of unique reads falling in the genomic region for each cell.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param process_n Number of regions to load into memory at a time, per thread.
#' Processing more regions at once can be faster but uses more memory.
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @param verbose Display messages
#'
#' @return Returns a sparse matrix
#'
#' @export
#' @importFrom SeuratObject RowMergeSparseMatrices
#'
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="scUtils")
#' ppath <- system.file("extdata", "peaks.rds", package="scUtils")
#' peaks <- readRDS(ppath)
#' fragments <- CreateFragmentObject(fpath)
#' FeatureMatrix(
#'  fragments = fragments,
#'   features = peaks
#' )
FeatureMatrix <- function(
        fragments,
        features,
        cells = NULL,
        process_n = 2000,
        sep = c("-", "-"),
        verbose = TRUE
) {
    if (!inherits(x = features, what = "GRanges")) {
        if (inherits(x = features, what = "character")) {
            features <- StringToGRanges(
                regions = features,
                sep = sep
            )
        } else {
            stop("features should be a GRanges object")
        }
    }
    if (!inherits(x = fragments, what = "list")) {
        if (inherits(x = fragments, what = "Fragment")) {
            fragments <- list(fragments)
        } else {
            stop("fragments should be a list of Fragment objects")
        }
    }
    # if cells is not NULL, iterate over all fragment objects
    # and find which objects contain cells that are requested
    if (!is.null(x = cells)) {
        obj.use <- c()
        for (i in seq_along(along.with = fragments)) {
            if (is.null(x = Cells(fragments[[i]]))) {
                # cells information not set for fragment object
                obj.use <- c(obj.use, i)
            } else if (any(cells %in% Cells(x = fragments[[i]]))) {
                obj.use <- c(obj.use, i)
            }
        }
    } else {
        obj.use <- seq_along(along.with = fragments)
    }
    # create a matrix from each fragment file
    mat.list <- sapply(
        X = obj.use,
        FUN = function(x) {
            SingleFeatureMatrix(
                fragment = fragments[[x]],
                features = features,
                cells = cells,
                sep = sep,
                verbose = verbose,
                process_n = process_n
            )
        })
    # merge all the matrices
    if (length(x = mat.list) == 1) {
        return(mat.list[[1]])
    } else {
        featmat <- Reduce(f = RowMergeSparseMatrices, x = mat.list)
        return(featmat)
    }
}


#' Run groupCommand for the first n lines, convert the cell barcodes in the file
#' to the cell names that appear in the fragment object, and subset the output to
#' cells present in the fragment object
#'
#' Every cell in the fragment file will be present in the output dataframe. If
#' the cell information is not set, every cell barcode that appears in the first
#' n lines will be present.
#'
#' @param fragments A Fragment object
#' @param n Number of lines to read from the beginning of the fragment file
#' @param verbose Display messages
#'
#' @return Returns a data.frame
ExtractFragments <- function(fragments, n = NULL, verbose = TRUE) {
    fpath <- GetFragmentData(object = fragments, slot = "path")
    if (isRemote(x = fpath)) {
        stop("Remote fragment files not supported")
    }
    fpath <- normalizePath(path = fpath, mustWork = TRUE)
    cells <- GetFragmentData(object = fragments, slot = "cells")
    if (!is.null(x = cells)) {
        cells.use <- as.character(x = cells)
    } else {
        cells.use <- NULL
    }
    verbose <- as.logical(x = verbose)
    n <- SetIfNull(x = n, y = 0)
    n <- as.numeric(x = n)
    n <- round(x = n, digits = 0)
    counts <- groupCommand(
        fragments = fpath,
        some_whitelist_cells = cells.use,
        max_lines = n,
        verbose = verbose
    )
    # convert cell names
    if (!is.null(x = cells)) {
        # every cell will be present in the output, even if 0 counts
        converter <- names(x = cells)
        names(x = converter) <- cells
        counts$CB <- converter[counts$CB]
    }
    return(counts)
}



# Generate cut matrix for many regions
#
# Run CutMatrix on multiple regions and add them together.
# Assumes regions are pre-aligned.
#
# @param object A Fragment object
# @param regions A set of GRanges
# @param group.by Name of grouping variable to use
# @param fragments A list of Fragment objects
# @param assay Name of the assay to use
# @param cells Vector of cells to include
# @param verbose Display messages
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom SeuratObject DefaultAssay
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels
MultiRegionCutMatrix <- function(
        object,
        regions,
        group.by = NULL,
        fragments = NULL,
        assay = NULL,
        cells = NULL,
        verbose = FALSE
) {
    fragments <- SetIfNull(x = fragments, y = object)
    
    if (length(x = fragments) == 0) {
        stop("No fragment files present in parameter 'object' or 'fragments'!")
    }
    
    res <- list()
    
    for (i in seq_along(along.with = fragments)) {
        frag.path <- GetFragmentData(object = fragments[[i]], slot = "path")
        cellmap <- GetFragmentData(object = fragments[[i]], slot = "cells")
        if (is.null(x = cellmap)) {
            cellmap <- as.character(fragments[[i]]@cells)
            names(x = cellmap) <- cellmap
        }
        tabix.file <- Rsamtools::TabixFile(file = frag.path)
        open(con = tabix.file)
        # remove regions that aren't in the fragment file
        common.seqlevels <- intersect(
            x = seqlevels(x = regions),
            y = Rsamtools::seqnamesTabix(file = tabix.file)
        )
        regions <- GenomeInfoDb::keepSeqlevels(
            x = regions,
            value = common.seqlevels,
            pruning.mode = "coarse"
        )
        cm <- SingleFileCutMatrix(
            cellmap = cellmap,
            tabix.file = tabix.file,
            region = regions,
            verbose = verbose
        )
        close(con = tabix.file)
        res[[i]] <- cm
    }
    # each matrix contains data for different cells at same positions
    # bind all matrices together
    res <- do.call(what = rbind, args = res)
    return(res)
}


#' Find transcriptional start sites
#'
#' Get the TSS positions from a set of genomic ranges containing gene positions.
#' Ranges can contain exons, introns, UTRs, etc, rather than the whole
#' transcript. Only protein coding gene biotypes are included in output.
#'
#' @param ranges A GRanges object containing gene annotations.
#' @param biotypes Gene biotypes to include. If NULL, use all biotypes in the
#' supplied gene annotation.
#'
#' @return transcriptional start sites
#'
#' @export
#' @importFrom GenomicRanges resize
#' @importFrom S4Vectors mcols
GetTSSPositions <- function(ranges, biotypes = "protein_coding") {
    if (!("gene_biotype" %in% colnames(x = mcols(x = ranges)))) {
        stop("Gene annotation does not contain gene_biotype information")
    }
    if (!is.null(x = biotypes)){
        ranges <- ranges[ranges$gene_biotype == "protein_coding"]
    }
    gene.ranges <- CollapseToLongestTranscript(ranges = ranges)
    # shrink to TSS position
    tss <- resize(gene.ranges, width = 1, fix = 'start')
    return(tss)
}


#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
CollapseToLongestTranscript <- function(ranges) {
    range.df <- as.data.table(x = ranges)
    range.df$strand <- as.character(x = range.df$strand)
    range.df$strand <- ifelse(
        test = range.df$strand == "*",
        yes = "+",
        no = range.df$strand
    )
    collapsed <- range.df[
        , .(unique(seqnames),
            min(start),
            max(end),
            strand[[1]],
            gene_biotype[[1]],
            gene_name[[1]]),
        "gene_id"
    ]
    colnames(x = collapsed) <- c(
        "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
    )
    collapsed$gene_name <- make.unique(names = collapsed$gene_name)
    gene.ranges <- makeGRangesFromDataFrame(
        df = collapsed,
        keep.extra.columns = TRUE
    )
    return(gene.ranges)
}



#' Extend
#'
#' Resize GenomicRanges upstream and or downstream.
#' From \url{https://support.bioconductor.org/p/78652/}
#'
#' @param x A range
#' @param upstream Length to extend upstream
#' @param downstream Length to extend downstream
#' @param from.midpoint Count bases from region midpoint,
#' rather than the 5' or 3' end for upstream and downstream
#' respectively.
#'
#' @importFrom GenomicRanges trim
#' @importFrom BiocGenerics start strand end width
#' @importMethodsFrom GenomicRanges strand start end width
#' @importFrom IRanges ranges IRanges "ranges<-"
#' @export
#' @concept utilities
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
Extend <- function(
        x,
        upstream = 0,
        downstream = 0,
        from.midpoint = FALSE
) {
    if (any(strand(x = x) == "*")) {
        warning("'*' ranges were treated as '+'")
    }
    on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
    if (from.midpoint) {
        midpoints <- start(x = x) + (width(x = x) / 2)
        new_start <- midpoints - ifelse(
            test = on_plus, yes = upstream, no = downstream
        )
        new_end <- midpoints + ifelse(
            test = on_plus, yes = downstream, no = upstream
        )
    } else {
        new_start <- start(x = x) - ifelse(
            test = on_plus, yes = upstream, no = downstream
        )
        new_end <- end(x = x) + ifelse(
            test = on_plus, yes = downstream, no = upstream
        )
    }
    ranges(x = x) <- IRanges(start = new_start, end = new_end)
    x <- trim(x = x)
    return(x)
}


# Generate matrix of integration sites
#
# Generates a cell-by-position matrix of Tn5 integration sites
# centered on a given region (usually a DNA sequence motif). This
# matrix can be used for downstream footprinting analysis.
#
# @param cellmap A mapping of cell names in the fragment file to cell names in
# the Seurat object. Should be a named vector where each element is a cell name
# that appears in the fragment file and the name of each element is the
# name of the cell in the Seurat object.
# @param region A set of GRanges containing the regions of interest
# @param cells Which cells to include in the matrix. If NULL, use all cells in
# the cellmap
# @param tabix.file A \code{\link[Rsamtools]{TabixFile}} object.
# @param verbose Display messages
#' @importFrom Matrix sparseMatrix
#' @importFrom Rsamtools TabixFile
#' @importMethodsFrom GenomicRanges width start end
# @return Returns a sparse matrix
SingleFileCutMatrix <- function(
        cellmap,
        region,
        cells = NULL,
        tabix.file,
        verbose = TRUE
) {
    # if multiple regions supplied, must be the same width
    cells <- SetIfNull(x = cells, y = names(x = cellmap))
    if (length(x = region) == 0) {
        return(NULL)
    }
    fragments <- GetReadsInRegion(
        region = region,
        cellmap = cellmap,
        cells = cells,
        tabix.file = tabix.file,
        verbose = verbose
    )
    start.lookup <- start(x = region)
    names(start.lookup) <- seq_along(region)
    # if there are no reads in the region
    # create an empty matrix of the correct dimension
    if (nrow(x = fragments) == 0) {
        cut.matrix <- sparseMatrix(
            i = NULL,
            j = NULL,
            dims = c(length(x = cells), width(x = region)[[1]])
        )
    } else {
        fragstarts <- start.lookup[fragments$ident] + 1
        cut.df <- data.frame(
            position = c(fragments$start, fragments$end) - fragstarts,
            cell = c(fragments$cell, fragments$cell),
            stringsAsFactors = FALSE
        )
        cut.df <- cut.df[
            (cut.df$position > 0) & (cut.df$position <= width(x = region)[[1]]),
        ]
        cell.vector <- seq_along(along.with = cells)
        names(x = cell.vector) <- cells
        cell.matrix.info <- cell.vector[cut.df$cell]
        cut.matrix <- sparseMatrix(
            i = cell.matrix.info,
            j = cut.df$position,
            x = 1,
            dims = c(length(x = cells), width(x = region)[[1]])
        )
    }
    rownames(x = cut.matrix) <- cells
    colnames(x = cut.matrix) <- seq_len(width(x = region)[[1]])
    return(cut.matrix)
}


#' GetReadsInRegion
#'
#' Extract reads for each cell within a given genomic region or set of regions
#'
#' @param cellmap A mapping of cell names in the fragment file to cell names in
#' the Seurat object. Should be a named vector where each element is a cell name
#' that appears in the fragment file and the name of each element is the
#' name of the cell in the Seurat object.
#' @param region A genomic region, specified as a string in the format
#' 'chr:start-end'. Can be a vector of regions.
#' @param tabix.file A TabixFile object.
#' @param cells Cells to include. Default is all cells present in the object.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#'
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom SeuratObject Idents
#' @importFrom fastmatch fmatch
#'
#' @return Returns a data frame
GetReadsInRegion <- function(
        cellmap,
        region,
        tabix.file,
        cells = NULL,
        verbose = TRUE,
        ...
) {
    file.to.object <- names(x = cellmap)
    names(x = file.to.object) <- cellmap
    
    if (verbose) {
        message("Extracting reads in requested region")
    }
    if (!is(object = region, class2 = "GRanges")) {
        region <- StringToGRanges(regions = region, ...)
    }
    # remove regions that aren't in the fragment file
    common.seqlevels <- intersect(
        x = seqlevels(x = region),
        y = seqnamesTabix(file = tabix.file)
    )
    region <- keepSeqlevels(
        x = region,
        value = common.seqlevels,
        pruning.mode = "coarse"
    )
    reads <- scanTabix(file = tabix.file, param = region)
    reads <- TabixOutputToDataFrame(reads = reads)
    reads <- reads[
        fmatch(x = reads$cell, table = cellmap, nomatch = 0L) > 0,
    ]
    # convert cell names to match names in object
    reads$cell <- file.to.object[reads$cell]
    if (!is.null(x = cells)) {
        reads <- reads[reads$cell %in% cells, ]
    }
    if (nrow(reads) == 0) {
        return(reads)
    }
    reads$length <- reads$end - reads$start
    return(reads)
}


# TabixOutputToDataFrame
#
# Create a single dataframe from list of character vectors
#
# @param reads List of character vectors (the output of \code{\link{scanTabix}})
# @param record.ident Add a column recording which region the reads overlapped
# with
#' @importFrom stringi stri_split_fixed
#' @importFrom S4Vectors elementNROWS
TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
    if (record.ident) {
        nrep <- elementNROWS(x = reads)
    }
    reads <- unlist(x = reads, use.names = FALSE)
    if (length(x = reads) == 0) {
        df <- data.frame(
            "chr" = "",
            "start" = "",
            "end" = "",
            "cell" = "",
            "count" = ""
        )
        df <- df[-1, ]
        return(df)
    }
    reads <- stri_split_fixed(str = reads, pattern = "\t")
    n <- length(x = reads[[1]])
    unlisted <- unlist(x = reads)
    e1 <- unlisted[n * (seq_along(along.with = reads)) - (n - 1)]
    e2 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 2)])
    e3 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 3)])
    e4 <- unlisted[n * (seq_along(along.with = reads)) - (n - 4)]
    e5 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 5)])
    df <- data.frame(
        "chr" = e1,
        "start" = e2,
        "end" = e3,
        "cell" = e4,
        "count" = e5,
        stringsAsFactors = FALSE,
        check.rows = FALSE,
        check.names = FALSE
    )
    if (record.ident) {
        df$ident <- rep(x = seq_along(along.with = nrep), nrep)
    }
    return(df)
}


#' Apply function to integration sites per base per group
#'
#' Perform colSums on a cut matrix with cells in the rows
#' and position in the columns, for each group of cells
#' separately.
#'
#' @param mat A cut matrix. See CutMatrix
#' @param groups A vector of group identities, with the name
#' of each element in the vector set to the cell name.
#' @param fun Function to apply to each group of cells.
#' For example, colSums or colMeans.
#' @param group.scale.factors Scaling factor for each group. Should
#' be computed using the number of cells in the group and the average number of
#' counts in the group.
#' @param normalize Perform sequencing depth and cell count normalization
#' @param scale.factor Scaling factor to use. If NULL (default), will use the
#' median normalization factor for all the groups.
#'
#' @return ApplyMatrixByGroup
#'
#' @importFrom stats median
ApplyMatrixByGroup <- function(
        mat,
        groups,
        fun,
        normalize = TRUE,
        group.scale.factors = NULL,
        scale.factor = NULL
) {
    if (normalize) {
        if (is.null(x = group.scale.factors) | is.null(x = scale.factor)) {
            stop("If normalizing counts, supply group scale factors")
        }
    }
    all.groups <- as.character(x = unique(x = groups))
    if (any(is.na(x = groups))) {
        all.groups <- c(all.groups, NA)
    }
    ngroup <- length(x = all.groups)
    npos <- ncol(x = mat)
    
    group <- unlist(
        x = lapply(X = all.groups, FUN = function(x) rep(x, npos))
    )
    position <- rep(x = as.numeric(x = colnames(x = mat)), ngroup)
    count <- vector(mode = "numeric", length = npos * ngroup)
    
    for (i in seq_along(along.with = all.groups)) {
        grp <- all.groups[[i]]
        if (is.na(x = grp)) {
            pos.cells <- names(x = groups)[is.na(x = groups)]
        } else {
            pos.cells <- names(x = groups)[groups == all.groups[[i]]]
        }
        if (length(x = pos.cells) > 1) {
            totals <- fun(x = mat[pos.cells, ])
        } else {
            totals <- mat[pos.cells, ]
        }
        count[((i - 1) * npos + 1):((i * npos))] <- totals
    }
    
    # construct dataframe
    coverages <- data.frame(
        "group" = group, "position" = position, "count" = count,
        stringsAsFactors = FALSE
    )
    
    if (normalize) {
        scale.factor <- SetIfNull(
            x = scale.factor, y = median(x = group.scale.factors)
        )
        coverages$norm.value <- coverages$count /
            group.scale.factors[coverages$group] * scale.factor
    } else {
        coverages$norm.value <- coverages$count
    }
    return(coverages)
}


#' Compute total fragment counts in genomic regions for update single cell
#' information table.
#'
#' Compute total fragment counts in genomic regions for every cell.
#'
#' @param fragment A fragment object.
#' @param GR Genomic regions saved in Granges.
#' @param csvFile The csv file used for saving single cell information.
#' @param csvOutputFile The updated csv file.
#' @param name The column name of new data.
#' @param process_n Number of regions to load into memory at a time, per thread.
#' Processing more regions at once can be faster but uses more memory.
#' Default: 5000.
#'
#' @return a vector of counts number for every cells.
#'
#' @export
#'
#' @importFrom Matrix colSums
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr full_join
FragInRegions <- function (fragment = NULL,
                           GR = NULL,
                           csvFile = NULL,
                           csvOutputFile = NULL,
                           name = NULL,
                           process_n = 2000) {
    # detect cells
    if (is.null(csvFile)) {
        Stop("csvFile must be specified")
    }else{
        print("csv file is specified, detecting cells in csv file.")
        metadata <- read.csv(
            file = csvFile,
            header = TRUE
        )
        cells <- metadata$barcode
        cellNum <- length(cells)
        print(paste0("Now, processing ", cellNum, " cells......"))
    }
    
    if (is.null(csvOutputFile)) {
        Stop("csvOutputFile must be specified")
    }
    
    # get feature x cell matrix for all cells with signal
    frag_in_GR <- FeatureMatrix(fragments = fragment, features = GR, process_n = process_n)
    
    record <- dim(frag_in_GR)
    print(paste0("In all ", record[1], " genome regions, ", record[2], " cells have signal."))
    
    # compute total fragment counts for every cell
    frag_in_GR <- Matrix::colSums(frag_in_GR)
    
    # match cells
    frag_in_GR_detected <- frag_in_GR[which(names(frag_in_GR) %in% cells)]
    
    print(paste0("Find ", length(frag_in_GR_detected), " cells which match the input csv file."))
    
    # update csvfiles
    df <- tibble::rownames_to_column(data.frame(frag_in_GR_detected))
    colnames(df) <- c("barcode", name)
    metadata <- full_join(x = metadata, y = df, by = "barcode")
    metadata[is.na(metadata)] <- 0
    write.csv(x = metadata, file = csvOutputFile, quote = FALSE)
}



#' Quality Control for Nucleosome Signal per Cell
#'
#' Calculate the strength of the nucleosome signal per cell.
#' Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to
#' fragments < 147 bp (nucleosome-free)
#'
#' @param frags A Fragment object. If the Fragment object contains "cells" slot,
#' the output will be filtered by "cells" slot. if not, all cells are reported.
#' @param n Number of lines to read from the fragment file. If NULL, read all
#' lines. Default: (cellnumber * 5e3) or 4e7.
#' @param verbose Display messages, TRUE or FALSE
#' @param ... Arguments passed to other functions
#'
#' @return Returns a dataframe for mononucleosomal, nucleosome_free,
#' nucleosome_signal and nucleosome_percentile of each cell.
#'
#' @importFrom stats ecdf
#' @importFrom fastmatch fmatch
#'
#' @export
#' @concept QualityControl
nsQC <- function(frags = NULL,
                 n = NULL,
                 verbose = TRUE,
                 ...){
    
    if (!inherits(x = frags, what = "Fragment")) {
        stop("The input frags shoule be a 'Fragment' object!")
    }
    
    if (is.null(n)) {
        if (is.null(frags@cells)) {
            n = 4e7
        }else{
            n = length(frags@cells) * 5e3
        }
    }
    
    verbose <- as.logical(x = verbose)
    
    counts <- ExtractFragments(
        fragments = frags,
        n = n,
        verbose = verbose
    )
    
    if (is.null(fragments@cells)) {
        cells.keep <- fastmatch::fmatch(x = counts$CB,
                                        table = counts$CB,
                                        nomatch = 0L)
    } else {
        cells.keep <- fastmatch::fmatch(x = counts$CB,
                                        table = as.character(fragments@cells),
                                        nomatch = 0L)
    }
    
    rownames(x = counts) <- counts$CB
    
    counts <- counts[
        cells.keep > 0, c("mononucleosomal", "nucleosome_free")
    ]
    af <- counts
    
    af$nucleosome_signal <- af$mononucleosomal / af$nucleosome_free
    
    e.dist <- ecdf(x = af$nucleosome_signal)
    
    af$nucleosome_percentile <- round(
        x = e.dist(af$nucleosome_signal),
        digits = 2
    )
    
    return(af)
}


#' Compute TSS enrichment score per cell
#'
#' Compute the transcription start site (TSS) enrichment score for each cell,
#' as defined by ENCODE:
#' \url{https://www.encodeproject.org/data-standards/terms/}.
#'
#' @param object A Fragment object, must contain slot 'cells'.
#' @param gene.annotation A GRanges object containing the TSS positions.
#' @param n Number of TSS positions to use. This will select the first _n_
#' TSSs from the set. If NULL, use all TSSs (slower).
#' @param cells A vector of cells to include. If NULL (default), use all cells
#' in the object
#' @param process_n Number of regions to process at a time if using \code{fast}
#' option.
#' @param verbose Display messages
#'
#' @importFrom Matrix rowMeans
#' @importFrom methods slot
#' @importFrom stats ecdf
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges start width strand
#' @importFrom SeuratObject DefaultAssay
#'
#' @return Returns two matrix
#' @export
#' @concept QualityControl
tssQC <- function (object = NULL,
                   gene.annotation = NULL,
                   n = NULL,
                   cells = NULL,
                   process_n = 2000,
                   verbose = TRUE) {
    
    # first check that fragments are present
    if (!inherits(x = object, what = "Fragment")) {
        stop("The input object shoule be a 'Fragment' object!")
    }
    
    # convert to a list
    object <- list(object)
    
    # check tss.positions
    if (is.null(x = gene.annotation)) {
        stop("No gene annotations present in assay")
    }
    
    tss.positions <- GetTSSPositions(ranges = gene.annotation)
    
    
    if (!is.null(x = n)) {
        if (n > length(x = tss.positions)) {
            n <- length(x = tss.positions)
        }
        tss.positions <- tss.positions[seq(n), ]
    }
    
    # exclude chrM
    sn <- seqnames(x = tss.positions)
    tss.positions <- tss.positions[!as.character(sn) %in% c("chrM", "Mt", "MT")]
    
    regions <- Extend(
        x = tss.positions,
        upstream = 1000,
        downstream = 1000,
        from.midpoint = TRUE
    )
    
    if (length(x = regions) == 0) {
        stop("No gene supplied")
    }
    
    # split into strands
    on_plus <- GenomicRanges::strand(x = regions) == "+" | GenomicRanges::strand(x = regions) == "*"
    plus.strand <- regions[on_plus, ]
    minus.strand <- regions[!on_plus, ]
    
    if (verbose) {
        message("Finding + strand cut sites")
    }
    cut.matrix.plus <- MultiRegionCutMatrix(
        regions = plus.strand,
        object = object,
        assay = assay,
        cells = cells,
        verbose = verbose
    )
    
    if (verbose) {
        message("Finding - strand cut sites")
    }
    cut.matrix.minus <- MultiRegionCutMatrix(
        regions = minus.strand,
        object = object,
        assay = assay,
        cells = cells,
        verbose = FALSE
    )
    
    # reverse minus strand and add together
    if (is.null(x = cut.matrix.plus)) {
        full.matrix <- cut.matrix.minus[, rev(x = colnames(x = cut.matrix.minus))]
    } else if (is.null(x = cut.matrix.minus)) {
        full.matrix <- cut.matrix.plus
    } else {
        full.matrix <- cut.matrix.plus + cut.matrix.minus[, rev(
            x = colnames(x = cut.matrix.minus)
        )]
    }
    
    # rename so 0 is center
    region.width <- GenomicRanges::width(x = regions)[[1]]
    midpoint <- round(x = (region.width / 2))
    colnames(full.matrix) <- seq_len(length.out = region.width) - midpoint
    
    cutmatrix <- full.matrix
    
    if (verbose) {
        message("Computing mean insertion frequency in flanking regions")
    }
    
    flanking.mean <- Matrix::rowMeans(x = cutmatrix[, c(seq(from = 1, to = 100), seq(from = 1902, to = 2001))])
    flanking.mean[is.na(x = flanking.mean)] <- 0
    flanking.mean[flanking.mean == 0] <- mean(flanking.mean, na.rm = TRUE)
    
    if (verbose) {
        message("Normalizing TSS score")
    }
    norm.matrix <- cutmatrix / flanking.mean
    
    # compute TSS enrichment
    TSS.enrichment <- Matrix::rowMeans(x = norm.matrix[, 500:1500], na.rm = TRUE)
    e.dist <- ecdf(x = TSS.enrichment)
    
    # compute TSS percentile
    TSS.percentile <- round(
        x = e.dist(TSS.enrichment),
        digits = 2
    )
    
    expected.insertions <- rep(1, ncol(x = cutmatrix))
    expected.insertions <- t(x = as.matrix(x = expected.insertions))
    rownames(x = expected.insertions) <- "expected"
    
    motif.vec <- t(x = matrix(
        data = c(
            rep(x = 0, 1000),
            1,
            rep(x = 0, 1000)
        )
    )
    )
    rownames(x = motif.vec) <- "motif"
    
    norm.matrix <- rbind(norm.matrix, expected.insertions)
    norm.matrix <- rbind(norm.matrix, motif.vec)
    
    positionEnrichment <- norm.matrix
    
    out <- list(TSS.enrichment = TSS.enrichment,
                TSS.percentile = TSS.percentile,
                positionEnrichment = positionEnrichment)
    
    return(out)
    
}

#' Plot signal enrichment around TSSs
#'
#' Plot the normalized TSS enrichment score at each position relative to the
#' TSS.
#'
#' @param tssInfo TSS enrichment information from function 'tssQC'
#' @param threshold Threshold to separate low and high quality cells
#'
#' @importFrom Matrix colMeans
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab theme_classic
#' ggtitle facet_wrap theme element_blank
#'
#' @return Returns a \code{\link[ggplot2]{ggplot2}} object
#' @export
#' @concept QualityControl
plotTssQC <- function (tssInfo = NULL,
                       threshold = 3) {
    if (is.null(tssInfo)) {
        stop("tssInfo is NULL!")
    }
    
    groups <- ifelse(tssInfo$TSS.enrichment > threshold, 'High', 'Low')
    
    enrichment.matrix <- tssInfo$positionEnrichment
    
    # remove motif and expected
    enrichment.matrix <- enrichment.matrix[seq((nrow(x = enrichment.matrix) - 2)), ]
    
    # average the signal per group per base
    groupmeans <- ApplyMatrixByGroup(
        mat = enrichment.matrix,
        groups = groups,
        fun = colMeans,
        normalize = FALSE
    )
    
    p <- ggplot(data = groupmeans,
                mapping = aes(x = position, y = norm.value, color = group)) +
        geom_line(stat = "identity", size = 0.2) +
        facet_wrap(facets = ~group) +
        xlab("Distance from TSS (bp)") +
        ylab(label = "Mean TSS enrichment score") +
        theme_classic() +
        theme(
            legend.position = "none",
            strip.background = element_blank()
        ) +
        ggtitle("TSS enrichment")
    
    return(p)
}











