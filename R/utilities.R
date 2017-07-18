#' Using Rbowtie to build genome reference.
#'
#' @param ref_dir A full path of the fa file(containing fa file).
#' @param out_dir A path you will save the ref.
#' @param file_prefix A prefix for your ref file, default "index".
#' @return to be continued!!!
#' @export
Rbowtie_ref_build <- function(ref_dir, out_dir, file_prefix) {
  Rbowtie::bowtie_build(references = ref_dir, outdir = out_dir, prefix = file_prefix,
                        force = TRUE)
}


#' Using GenomicRanges function to find overlap of 2 bed files.
#'
#' Find the overlap peak of the Input1 file compared with Input2 file.
#' The two input file must have the same column, because the function will processing 3 and 6 column file with
#' difference methods according to the strand column.
#' The function using 0-based coordinate as bed file.
#' @param Input1 The peak bed file(3-6 columns) that you want to find difference compared with the reference bed file.
#' @param Input2 The reference peak bed file.
#' @param bed_output The output bed file.
#' @param n.col How many columns of the input bed file.
#' @export
BedIntersect <- function(Input1, Input2, bed_output, n.col){
  bed_1 <- read.table(Input1, header = FALSE)
  bed_2 <- read.table(Input2, header = FALSE)
  if(n.col == 3){
    colnames(bed_1) <- c('chr', 'start', 'end')
    colnames(bed_2) <- c('chr', 'start', 'end')
    data1 <- with(bed_1, GRanges(chr, IRanges(start, end)))
    data2 <- with(bed_2, GRanges(chr, IRanges(start, end)))
    output_data <- subsetByOverlaps(data1, data2, ignore.strand = TRUE)
    write.table(output_data, file = bed_output,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }else if(n.col == 6){
    colnames(bed_1) <- c('chr', 'start', 'end', 'id', 'score', 'strand')
    colnames(bed_2) <- c('chr', 'start', 'end', 'id', 'score', 'strand')
    data1 <- with(bed_1, GRanges(chr, IRanges(start, end), strand, score, id = id))
    data2 <- with(bed_2, GRanges(chr, IRanges(start, end), strand, score, id = id))
    output_data <- subsetByOverlaps(data1, data2, ignore.strand = FALSE)
    write.table(output_data, file = bed_output,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}























