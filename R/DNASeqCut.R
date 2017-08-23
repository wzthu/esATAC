DNASeqCut <- R6::R6Class(
  classname = "DNASeqCut",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, seqRef = NULL, seqOutput = NULL,
                          bedInput = NULL, format = NULL,
                          editable = FALSE){
      super$initialize("DNASeqCut", editable, list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){ # get parameter from class PeakCallingFseq
        private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
      }else{
        private$paramlist[["bedInput"]] <- bedInput
      }
      private$paramlist[["seqRef"]] <- seqRef
      private$paramlist[["seqOutput"]] <- seqOutput
      private$paramlist[["format"]] <- format
      # parameter check
      private$paramValidation()

    } # initialization end

  ), # public end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("Bed file:%s", private$paramlist[["bedInput"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["seqOutput"]]))
      private$writeLog(sprintf("Format:%s", private$paramlist[["format"]]))
      Sequence_Cut(private$paramlist[["seqRef"]], private$paramlist[["seqOutput"]],
                   private$paramlist[["bedInput"]], private$paramlist[["format"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["seqRef"]])){
        stop("Parameter seqRef is required!")
      }
      if(is.null(private$paramlist[["seqOutput"]])){
        stop("Parameter seqOutput is required!")
      }
      if(is.null(private$paramlist[["bedInput"]])){
        stop("Parameter bedInput is required!")
      }
      if(is.null(private$paramlist[["format"]])){
        stop("Parameter format is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["seqRef"]])
      private$checkFileExist(private$paramlist[["bedInput"]])
      private$checkPathExist(private$paramlist[["seqOutput"]])
    } # checkAllPath end

  ) # private end

) # class end


# Cutting sequence according a bed file and save these sequence as fastq or fasta.
#
# In this program, the strand infomation will not be used.
# ref_path The reference fasta file.
# save_path Where you want save these sequences, only in fastq or fasta format.
# bed_path bed file.
# save_format Fastq or fasta.
Sequence_Cut <- function(ref_path, save_path, bed_path, save_format){
  ref <- Biostrings::readDNAStringSet(ref_path)
  gr_a <- ChIPseeker::readPeakFile(peakfile = bed_path, header = FALSE)
  chr_index <- names(ref)
  line_num <- length(gr_a)
  DNA_string <- Biostrings::DNAStringSet()
  for(i in seq(line_num)){
    index_num <- which(as.vector(GenomicRanges::seqnames(gr_a[i])) == chr_index)
    if(any(index_num)){
      DNA_string <- append(DNA_string, Biostrings::subseq(ref[index_num], GenomicRanges::start(gr_a[i]),
                                                          GenomicRanges::end(gr_a[i])), after=length(DNA_string))
    }else{
      stop("An unexpected chromatin detected!")
    }
  }
  names(DNA_string) <- paste(GenomicRanges::seqnames(gr_a), ":", GenomicRanges::start(gr_a),
                             "-", end(gr_a), sep = "")
  Biostrings::writeXStringSet(DNA_string, filepath = save_path,
                              append = FALSE, compress = FALSE,
                              compression_level = NA, format = save_format)
}





#' Cutting sequence according a bed file and save these sequence as fastq or fasta.
#'
#' In this program, the strand infomation will not be used.
#' @param seqRef The reference fasta file.
#' @param seqOutput Where you want save these sequences, only in fastq or fasta format.
#' @param bedInput Input bed file.
#' @param format Fastq or fasta.
DnaSeqCut <- function(atacProc = NULL, seqRef = NULL, seqOutput = NULL,
                      bedInput = NULL, format = NULL){
  tmp <- DNASeqCut$new(atacProc, seqRef, seqOutput, bedInput, format)
  tmp$process()
  return(tmp)
}
