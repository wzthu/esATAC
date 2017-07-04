#' using Rbowtie to build genome reference
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
