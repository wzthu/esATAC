################################################################################################
###############################This file using the following license############################
################################################################################################
##
##
## MIT License
##
## Copyright (c) 2021 Tim Stuart
##
##     Permission is hereby granted, free of charge, to any person obtaining a copy
##     of this software and associated documentation files (the "Software"), to deal
##     in the Software without restriction, including without limitation the rights
##     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
##     copies of the Software, and to permit persons to whom the Software is
##     furnished to do so, subject to the following conditions:
##
##         The above copyright notice and this permission notice shall be included in all
##         copies or substantial portions of the Software.
##
##         THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
##         IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
##         FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
##         AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
##         LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
##         OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##         SOFTWARE.
##
##################################################################################################

#' @include generics.R
#'
NULL


#' @title The Fragment class
#'
#' @description 
#' The Fragment class is designed to hold information needed for working with
#' fragment files. This class is from Signac.
#'
#' @slot path Path to the fragment file on disk.
#' See \url{https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments}
#' @slot hash A vector of two md5sums: first element is the md5sum of the
#' fragment file, the second element is the md5sum of the index.
#' @slot cells A named vector of cells where each element is the cell barcode
#' as it appears in the fragment file.
#'
#' @name Fragment-class
#' 
#' @exportClass Fragment
#' 
#' @concept fragments
#'
#' @examples
#' frag <- new("Fragment", path = "file path", hash = "file hash")
#'
#' @keywords internal
#'
Fragment <- setClass(
    Class = "Fragment",
    slots = list(
        path = "character",
        hash = "character",
        cells = "ANY"
    )
)


#' @title Get Fragment object data
#'
#' @description 
#' Extract data from a \code{\link{Fragment-class}} object
#'
#' @param object A \code{\link{Fragment}} object
#' @param slot Information to pull from object (path, hash, cells, prefix, suffix)
#' 
#' @return Fragment object
#' @examples
#' frag <- new("Fragment", path = "file path", hash = "file hash")
#' GetFragmentData(frag)
#'
#' @keywords internal
#'
GetFragmentData <- function(object, slot = "path") {
    return(slot(object = object, name = slot))
}


#' @title Return the first rows of a fragment file
#'
#' @describeIn 
#' Returns the first \code{n} rows of a fragment file. This allows the content
#' of a fragment file to be inspected.
#'
#' @param x a \code{Fragment} object
#' @param n an integer specifying the number of rows to return from the fragment
#' file
#' @param ... additional arguments passed to \code{\link[utils]{read.table}}
#'
#' @return The first \code{n} rows of a fragment file as a \code{data.frame}
#' with the following columns: chrom, start, end, barcode, readCount.
#' 
#' @importFrom utils head read.table
#' 
#' @method head Fragment
#'
#' @keywords internal
#'
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="scUtils")
#' frag <- CreateFragmentObject(fpath)
#' head(frag)
#' 
head.Fragment <- function(x, n = 6L, ...) {
    fpath <- GetFragmentData(object = x, slot = "path")
    df <- read.table(file = fpath, nrows = n, ...)
    if (ncol(x = df) == 5) {
        colnames(x = df) <- c("chrom", "start", "end", "barcode", "readCount")
    }
    return(df)
}




#' @title Count fragments
#'
#' @description 
#' Count total fragments per cell barcode present in a fragment file.
#'
#' @param fragments Path to a fragment file. If a list of fragment files is
#' provided, the total fragments for each cell barcode across all files will be
#' returned
#' @param cells Cells to include. If NULL, include all cells
#' @param max_lines Maximum number of lines to read from the fragment file. If
#' NULL, read all lines in the file.
#' @param verbose Display messages
#'
#' @return Returns a data.frame with the following columns:
#'
#' \itemize{
#'   \item{CB: the cell barcode}
#'   \item{frequency_count: total number of fragments sequenced for the cell}
#'   \item{mononucleosome: total number of fragments with length between 147 bp and 294 bp}
#'   \item{nucleosome_free: total number of fragments with length <147 bp}
#'   \item{reads_count: total number of reads sequenced for the cell}
#' }
#'
#' @keywords internal
#'
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="scUtils")
#' counts <- CountFragments(fragments = fpath)
#' 
CountFragments <- function(
        fragments,
        cells = NULL,
        max_lines = NULL,
        verbose = TRUE
) {
    if (!inherits(x = fragments, what = "list")) {
        fragments <- list(fragments)
    }
    for (i in seq_along(along.with = fragments)) {
        if (isRemote(x = i)) {
            stop("Remote fragment files not supported")
        }
        fragments[[i]] <- normalizePath(path = fragments[[i]], mustWork = TRUE)
        max_lines <- SetIfNull(x = max_lines, y = 0)
        verbose = as.logical(x = verbose)
        if (!is.null(x = cells)) {
            cells <- unique(x = cells)
        }
        counts <- groupCommand(
            fragments = fragments[[i]],
            some_whitelist_cells = cells,
            max_lines = max_lines,
            verbose = verbose
        )
        rownames(x = counts) <- counts$CB
        counts$CB <- NULL
        if (i == 1) {
            # first file
            allcounts <- counts
        } else {
            # merge
            common <- intersect(
                x = rownames(x = allcounts), y = rownames(x = counts)
            )
            allcounts[common, ] <- allcounts[common, ] + counts[common, ]
            missing_cells <- setdiff(x = rownames(x = counts), y = common)
            allcounts <- rbind(allcounts, counts[missing_cells, ])
        }
    }
    # reformat for backwards compatibility
    allcounts$CB <- rownames(x = allcounts)
    rownames(x = allcounts) <- NULL
    allcounts <- allcounts[, c(5, 1, 2, 3, 4)]
    return(allcounts)
}


#' @title Filter cells from fragment file
#'
#' @description 
#' Remove all fragments that are not from an allowed set of cell barcodes from
#' the fragment file. This will create a new file on disk that only contains
#' fragments from cells specified in the \code{cells} argument. The output file
#' is block gzip-compressed and indexed.
#'
#' @param fragments Path to a fragment file
#' @param cells A vector of cells to keep
#' @param outfile Name for output file
#' @param buffer_length Size of buffer to be read from the fragment file. This
#' must be longer than the longest line in the file.
#' @param verbose Display messages
#'
#' @return Filter cells from fragment file
#'
#' @importFrom Rsamtools bgzip indexTabix
#'
#' @keywords internal
#'
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="scUtils")
#' cells <- readLines(system.file("extdata", "cell_name.txt", package="scUtils"))
#' tmpf <- tempfile(fileext = ".gz")
#' FilterCells(
#'   fragments = fpath,
#'   cells = cells,
#'   outfile = tmpf
#' )
#' file.remove(tmpf)
#'
FilterCells <- function(
        fragments,
        cells,
        outfile = NULL,
        buffer_length = 256L,
        verbose = TRUE
) {
    if (isRemote(x = fragments)) {
        stop("Remote fragment files not supported")
    }
    fragments <- normalizePath(path = fragments, mustWork = TRUE)
    if (is.null(x = outfile)) {
        outfile <- paste0(fragments, ".filtered")
    }
    if (file.exists(outfile)) {
        warning("Output file already exists, file will be overwritten",
                immediate. = TRUE)
    }
    verbose <- as.logical(x = verbose)
    tf <- tempfile(pattern = "filtercells", tmpdir = tempdir(), fileext = "")
    buffer_length <- as.integer(x = buffer_length)
    filtered <- filterCells(
        fragments = fragments,
        keep_cells = cells,
        outfile = tf,
        buffer_length = buffer_length,
        verbose = verbose
    )
    if (filtered == 1) {
        stop("Error: cannot open requested file")
    }
    # bgzip and index output
    if (verbose) {
        message("\nCompressing filtered file")
    }
    bgzip(file = tf, dest = outfile, overwrite = TRUE)
    file.remove(tf)
    if (verbose) {
        message("Indexing fragment file")
    }
    idx <- indexTabix(file = outfile, format = "bed")
}


#' @title 
#' Split fragment file by cell identities
#' 
#' @description 
#' Splits a fragment file into separate files for each group of cells. If
#' splitting multiple fragment files containing common cell types, fragments
#' originating from different files will be appended to the same file for one
#' group of cell identities.
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param group.by Name of grouping variable to group cells by
#' @param idents List of identities to include
#' @param buffer_length Size of buffer to be read from the fragment file. This
#' must be longer than the longest line in the file.
#' @param outdir Directory to write output files
#' @param file.suffix Suffix to add to all file names (before file extension).
#' If splitting multiple fragment files without the \code{append} option set to
#' TRUE, an additional numeric suffix will be added to each file (eg, .1, .2).
#' @param append If splitting multiple fragment files, append cells from the
#' same group (eg cluster) to the same file. Note that this can cause the output
#' file to be unsorted.
#' @param verbose Display messages
#'
#' @return SplitFragments
#'
#' @keywords internal
#'
SplitFragments <- function(
        object,
        assay = NULL,
        group.by = NULL,
        idents = NULL,
        outdir = getwd(),
        file.suffix = "",
        append = TRUE,
        buffer_length = 256L,
        verbose = TRUE
) {
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    frags <- Fragments(object = object[[assay]])
    groups <- GetGroups(
        object = object,
        group.by = group.by,
        idents = idents
    )
    # replace space with underscore
    cells <- names(x = groups)
    groups <- gsub(pattern = " ", replacement = "_", x = groups)
    groups <- gsub(pattern = .Platform$file.sep, replacement = "_", x = groups)
    names(x = groups) <- cells
    buffer_length <- as.integer(x = buffer_length)
    file.suffix <- as.character(x = file.suffix)
    idents <- as.character(x = unname(obj = groups))
    unique_idents <- unique(x = idents)
    outdir <- normalizePath(path = outdir, mustWork = TRUE)

    # split cells from each fragment file
    # append to existing file when more than one fragment file used
    for (i in seq_along(along.with = frags)) {
        if (!append & (length(x = frags) > 1)) {
            suffix.use <- paste0(file.suffix, ".", i)
        } else {
            suffix.use <- file.suffix
        }
        fragpath <- GetFragmentData(object = frags[[i]], slot = "path")
        if (isRemote(x = fragpath)) {
            message("Remote fragment files not supported, skipping fragment file")
        } else {
            # convert cell names
            cellmap <- GetFragmentData(object = frags[[i]], slot = "cells")
            cell.in.frag <- cells %in% names(x = cellmap)
            cells.use <- cellmap[cells[cell.in.frag]]
            idents.use <- as.character(x = unname(obj = groups[names(x = cells.use)]))
            frag.cell.name <- as.character(x = unname(obj = cells.use))

            if (verbose) {
                message("Processing file ", fragpath)
            }
            splitfiles <- splitFragments(
                fragments = fragpath,
                outdir = paste0(outdir, .Platform$file.sep),
                suffix = suffix.use,
                append = append,
                cells = frag.cell.name,
                idents = idents.use,
                unique_idents = unique_idents,
                buffer_length = buffer_length,
                verbose = verbose
            )
            if (verbose) {
                message("\n")
            }
            if (splitfiles == 1) {
                stop("Error: cannot open requested file")
            }
        }
    }
}


#' @title Create a Fragment object
#' 
#' @description 
#' Create a \code{Fragment} object to store fragment file information.
#' This object stores a 32-bit MD5 hash of the fragment file and the fragment
#' file index so that any changes to the files on-disk can be detected. A check
#' is also performed to ensure that the expected cells are present in the
#' fragment file.
#'
#' @param path A path to the fragment file. The file should contain a tabix
#' index in the same directory.
#' @param cells A named character vector containing cell barcodes contained in
#' the fragment file. This does not need to be all cells in the fragment file,
#' but there should be no cells in the vector that are not present in the
#' fragment file. A search of the file will be performed until at least one
#' fragment from each cell is found. If NULL, don't check for expected cells.
#'
#' Each element of the vector should be a cell barcode that appears in the
#' fragment file, and the name of each element should be the corresponding cell
#' name in the object.
#' @param validate.fragments Check that expected cells are present in the
#' fragment file.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{ValidateCells}
#'
#' @return Fragment object
#'
#' @importFrom tools md5sum file_ext
#'
#' @keywords internal
#'
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="SignacSlim")
#' fragments <- CreateFragmentObject(fpath)
#' 
CreateFragmentObject <- function(
        path,
        cells = NULL,
        validate.fragments = TRUE,
        verbose = TRUE,
        ...
) {
    # check that file exists and is indexed
    # don't check if supplying remote file
    is.remote <- isRemote(x = path)
    if (is.remote) {
        validate.fragments <- FALSE
    }
    if (!file.exists(path) & !is.remote) {
        stop("Fragment file does not exist.")
    }
    index.file <- paste0(path, ".tbi")
    if (!file.exists(index.file) & !is.remote) {
        stop("Fragment file is not indexed.")
    }
    if (is.remote) {
        con <- gzcon(con = url(description = path))
    } else {
        con <- path
    }
    df <- readLines(con = con, n = 10000)
    for (i in df) {
        if (grepl(pattern = '^#', x = i)) {
            next
        } else {
            if (length(x = strsplit(x = i, split = "\t")[[1]]) != 5) {
                stop("Incorrect number of columns found in fragment file")
            } else {
                break
            }
        }
    }
    if (!is.null(x = cells)) {
        if (is.null(names(x = cells))) {
            # assume cells are as they appear in the assay
            names(x = cells) <- cells
        }
    }
    # compute hash of the file and index
    if (verbose) {
        message("Computing hash......")
    }
    if (!is.remote) {
        path <- normalizePath(path = path, mustWork = TRUE)
    }
    # will be NA if file remote
    hashes <- md5sum(files = c(path, index.file))
    # create object
    frags <- new(
        Class = "Fragment",
        path = path,
        hash = unname(obj = hashes),
        cells = cells
    )
    # validate cells
    if (!is.null(x = cells) & validate.fragments) {
        if (ValidateCells(object = frags, verbose = verbose, ...)) {
            return(frags)
        } else {
            stop("Not all cells requested could be found in the fragment file.")
        }
    } else {
        return(frags)
    }
}

#' @title Validate cells present in fragment file
#'
#' @description 
#' Search for a fragment from each cell that should exist in the fragment file.
#' Will iterate through chunks of the fragment file until at least one fragment
#' from each cell barcode requested is found.
#'
#' @param object A \code{\link{Fragment}} object
#' @param cells A character vector containing cell barcodes to search for.
#' If NULL, use the cells stored in the Fragment object.
#' @param tolerance Fraction of input cells that can be unseen before returning
#' TRUE. For example, \code{tolerance = 0.01} will return TRUE when 99% of cells
#' have observed fragments in the file. This can be useful if there are cells
#' present that have much fewer total counts, and would require extensive
#' searching before a fragment from those cells are found.
#' @param max.lines Maximum number of lines to read in without finding the
#' required number of cells before returning FALSE. Setting this value avoids
#' having to search the whole file if it becomes clear that the expected cells
#' are not present. Setting this value to NULL will enable an exhaustive search
#' of the entire file.
#' @param verbose Display messages
#'
#' @return ValidateCells
#'
#' @keywords internal
#'
#' @examples
#' print("see https://satijalab.org/signac/reference/validatecells")
ValidateCells <- function(
        object,
        cells = NULL,
        tolerance = 0.5,
        max.lines = 5e7,
        verbose = TRUE
) {
    cell_barcodes <- GetFragmentData(object = object, slot = "cells")
    cells <- SetIfNull(x = cells, y = cell_barcodes)
    if (is.null(x = cells)) {
        warning("No cells stored in object")
        return(TRUE)
    }
    max.lines <- SetIfNull(x = max.lines, y = 0)
    filepath <- GetFragmentData(object = object, slot = "path")
    filepath <- normalizePath(path = filepath, mustWork = TRUE)
    is.remote <- isRemote(x = filepath)
    find_n <- as.integer(x = length(x = cells) * (1 - tolerance))
    # if remote, return TRUE
    if (is.remote) {
        return(TRUE)
    }
    valid <- validateCells(
        fragments = filepath,
        cells = cells,
        find_n = find_n,
        max_lines = max.lines,
        verbose = verbose
    )
    return(valid)
}

#' @title Validate hashes for Fragment object
#'
#' @param object A \code{\link{Fragment}} object
#' @param verbose Display messages
#'
#' @return ValidateHash
#'
#'
#' @keywords internal
#'
#' @importFrom tools md5sum
#' 
ValidateHash <- function(object, verbose = TRUE) {
    path <- GetFragmentData(object = object, slot = "path")
    index.file <- paste0(path, ".tbi")
    is.remote <- isRemote(x = path)
    # if remote, return TRUE
    if (is.remote) {
        return(TRUE)
    }
    if (!all(file.exists(path, index.file))) {
        return(FALSE)
    }
    oldhash <- GetFragmentData(object = object, slot = "hash")
    if (verbose) {
        message("Computing hash")
    }
    newhash <- md5sum(files = c(path, index.file))
    return(all(newhash == oldhash))
}

#' @title Validate Fragment object
#'
#' @description 
#' Verify that the cells listed in the object exist in the fragment file
#' and that the fragment file or index have not changed since creating the
#' fragment object.
#'
#' @param object A \code{\link{Fragment}} object
#' @param verbose Display messages
#' @param ... Additional parameters passed to \code{\link{ValidateCells}}
#'
#' @return ValidateFragments
#'
#' @keywords internal
#'
#' @concept fragments
#' 
ValidateFragments <- function(
        object,
        verbose = TRUE,
        ...
) {
    valid.cells <- ValidateCells(object = object, verbose = verbose, ...)
    valid.hash <- ValidateHash(object = object, verbose = verbose, ...)
    return(valid.cells & valid.hash)
}

#' @title Set and get cell barcode information for a \code{\link{Fragment}} object
#' 
#' @describeIn 
#' This returns the names of cells in the object that are contained in the
#' fragment file. These cell barcodes may not match the barcodes present in the
#' fragment file. The \code{\link{Fragment}} object contains an internal mapping
#' of the cell names in the ChromatinAssay object to the cell
#' names in the fragment file, so that cell names can be changed in the
#' assay without needing to change the cell names on disk.
#' 
#' To access the cell names that are stored in the fragment file itself, use
#' \code{GetFragmentData(object = x, name = "cells")}.
#' 
#' @param x A Fragment object
#' @param ... Arguments passed to other methods
#' @return cell barcode information
#' 
#' @method Cells Fragment
#' @importFrom SeuratObject Cells
#'
#' @keywords internal
#'
Cells.Fragment <- function(x, ...) {
    cells <- slot(object = x, name = "cells")
    return(names(x = cells))
}


#' @title Cells
#' 
#' @param value A vector of cell names to store in the \code{\link{Fragment}}
#' object
#' 
#' @method Cells<- Fragment
#'
#' @keywords internal
#'
"Cells<-.Fragment" <- function(x, ..., value) {
    if (is.null(x = names(x = value))) {
        stop("Cells must be a named vector")
    }
    slot(object = x, name = "cells") <- value
    if (!ValidateCells(object = x, verbose = FALSE, ...)) {
        stop("Cells not present in fragment file")
    } else {
        return(x)
    }
}


#' @title Update the file path for a Fragment object
#'
#' @description 
#' Change the path to a fragment file store in a \code{\link{Fragment}}
#' object. Path must be to the same file that was used to create the fragment
#' object. An MD5 hash will be computed using the new path and compared to the
#' hash stored in the Fragment object to verify that the files are the same.
#'
#' @param object A \code{\link{Fragment}} object
#' @param new.path Path to the fragment file
#' @param verbose Display messages
#'
#' @return UpdatePath
#' 
#' @keywords internal
#'
#' @examples
#' print("see https://satijalab.org/signac/reference/updatepath")
UpdatePath <- function(object, new.path, verbose = TRUE) {
    new.is.remote <- isRemote(x = new.path)
    if (!new.is.remote) {
        new.path <- normalizePath(path = new.path, mustWork = TRUE)
        index.file <- paste0(new.path, ".tbi")
        if (!file.exists(new.path)) {
            stop("Fragment file not found")
        } else if (!file.exists(index.file)) {
            stop("Fragment file not indexed")
        }
    }
    old.path <- GetFragmentData(object = object, slot = "path")
    old.is.remote <- isRemote(x = old.path)
    if (identical(x = old.path, y = new.path)) {
        return(object)
    }
    if (!old.is.remote & new.is.remote) {
        warning("Replacing local file path with a remote file")
    }
    slot(object = object, name = "path") <- new.path
    if (ValidateHash(object = object, verbose = verbose)) {
        return(object)
    } else {
        stop("MD5 sum does not match previously computed sum")
    }
}

#' @title AssignFragCellnames
#' 
#' @description  If Fragment object does not contain cell names, assign to given list of names
#' Used in CreateChromatinAssay when adding a Fragment object with no cell names
#' 
#' @importFrom methods slot slot<-
#' 
#' @keywords internal
#'
AssignFragCellnames <- function(fragments, cellnames) {
    if (is.null(x = Cells(x = fragments))) {
        cells <- cellnames
        names(x = cells) <- cells
        slot(object = fragments, name = "cells") <- cells
    }
    return(fragments)
}
