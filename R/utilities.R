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


#' Cutting sequence according a bed file and save these sequence as fastq or fasta.
#'
#' In this program, the strand infomation will not be used.
#' @param ref_path The reference fasta file.
#' @param save_path Where you want save these sequences, only in fastq or fasta format.
#' @param bed_path bed file.
#' @param save_format Fastq or fasta.
Sequence_Cut <- function(ref_path, save_path, bed_path, save_format){
  ref <- Biostrings::readDNAStringSet(ref_path)
  gr_a <-read.table(bed_path, header = FALSE)
  gr_a <- gr_a[, 1:3]
  colnames(gr_a) <- c('chr', 'start', 'end')
  gr_a$start <- gr_a$start + 1
  gr_a$end <- gr_a$end
  chr_index <- names(ref)
  line_num <- nrow(gr_a)
  DNA_string <- Biostrings::DNAStringSet()
  for(i in seq(line_num)){
    index_num <- which(gr_a$chr[i] == chr_index)
    if(any(index_num)){
      DNA_string <- append(DNA_string, Biostrings::subseq(ref[index_num], gr_a$start[i], gr_a$end[i]),
                           after=length(DNA_string))
    }else{
      stop("An unexpected chromatin detected!")
    }
  }
  names(DNA_string) <- paste(gr_a$chr, ":", gr_a$start, "-", gr_a$end, sep = "")
  Biostrings::writeXStringSet(DNA_string, filepath = save_path,
                              append = FALSE, compress = FALSE,
                              compression_level = NA, format = save_format)
}

















