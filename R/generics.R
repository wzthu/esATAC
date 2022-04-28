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

#' Get the Fragment objects
#'
#' @param ... Arguments passed to other methods
#' @return Returns a list of \code{\link{Fragment}} objects. If there are
#' no Fragment objects present, returns an empty list.
#' @rdname Fragments
#' @export Fragments
Fragments <- function(object, ...) {
    UseMethod(generic = "Fragments", object = object)
}


#' @param value A \code{\link{Fragment}} object or list of Fragment objects
#'
#' @rdname Fragments
#' @export Fragments<-
#'
"Fragments<-" <- function(object, ..., value) {
    UseMethod(generic = 'Fragments<-', object = object)
}

#' Binarize counts
#'
#' Set counts >1 to 1 in a count matrix
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname BinarizeCounts
#' @export BinarizeCounts
#' @examples
#' print("see https://satijalab.org/signac/reference/binarizecounts")
BinarizeCounts <- function(object, ...) {
    UseMethod(generic = "BinarizeCounts", object = object)
}


#' Set and get cell barcode information for a Fragment object
#'
#' @param x A Seurat object
#' @param value A character vector of cell barcodes
#' @param ... Arguments passed to other methods
#' @return cell barcode information
#' @export Cells<-
#' @examples
#' print("see https://satijalab.org/signac/reference/allelefreq")
"Cells<-" <- function(x, ..., value) {
    UseMethod(generic = "Cells<-", object = x)
}



#' Get or set links information
#'
#' Get or set the genomic link information for a Seurat object
#'
#' @param ... Arguments passed to other methods
#' @return Links
#'
#' @rdname Links
#' @export Links
#' @examples
#' print("see https://satijalab.org/signac/articles/data_structures.html")
Links <- function(object, ...) {
    UseMethod(generic = "Links", object = object)
}

#' @param value A \code{\link[GenomicRanges]{GRanges}} object
#' @rdname Links
#' @export Links<-
"Links<-" <- function(object, ..., value) {
    UseMethod(generic = "Links<-", object = object)
}

#' Region enrichment analysis
#'
#' Count fragments within a set of regions for different groups of
#' cells.
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname RegionMatrix
#' @export RegionMatrix
#' @examples
#' print("see https://satijalab.org/signac/reference/regionheatmap")
RegionMatrix <- function(object, ...) {
    UseMethod(generic = "RegionMatrix", object = object)
}

#' Compute base composition information for genomic ranges
#'
#' Compute the GC content, region lengths, and dinucleotide base frequencies
#' for regions in the assay and add to the feature metadata.
#'
#' @param object A Seurat object, Assay object, or set of genomic ranges
#' @param ... Arguments passed to other methods
#' @return Returns a dataframe
#' @rdname RegionStats
#' @export RegionStats
#' @examples
#' print("see https://satijalab.org/signac/reference/regionstats")
RegionStats <- function(object, ...) {
    UseMethod(generic = "RegionStats", object = object)
}




