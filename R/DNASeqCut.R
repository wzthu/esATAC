DNASeqCut <- R6::R6Class(
  classname = "DNASeqCut",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, ref_path = NULL, save_path = NULL,
                          bed_path = NULL, save_format = NULL,
                          editable = FALSE){
      super$initialize("DNASeqCut", editable, list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){
        print("Parameter atacProc is not using now! We will add more functions in the future!")
      }
      private$paramlist[["ref_path"]] <- ref_path
      private$paramlist[["save_path"]] <- save_path
      private$paramlist[["bed_path"]] <- bed_path
      private$paramlist[["save_format"]] <- save_format
      # parameter check
      private$paramValidation()

    } # initialization end

  ), # public end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("Bed file:%s", private$paramlist[["bed_path"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["save_path"]]))
      private$writeLog(sprintf("Format:%s", private$paramlist[["save_format"]]))
      Sequence_Cut(private$paramlist[["ref_path"]], private$paramlist[["save_path"]],
                   private$paramlist[["bed_path"]], private$paramlist[["save_format"]])
    }, # processing end

    checkRequireParam = function(){
      if(is.null(private$paramlist[["ref_path"]])){
        stop("Parameter ref_path is required!")
      }
      if(is.null(private$paramlist[["save_path"]])){
        stop("Parameter save_path is required!")
      }
      if(is.null(private$paramlist[["bed_path"]])){
        stop("Parameter bed_path is required!")
      }
      if(is.null(private$paramlist[["save_format"]])){
        stop("Parameter save_format is required!")
      }
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["ref_path"]])
      private$checkFileExist(private$paramlist[["bed_path"]])
      private$checkPathExist(private$paramlist[["save_path"]])
    } # checkAllPath end

  ) # private end

) # class end


#' Cutting sequence according a bed file and save these sequence as fastq or fasta.
#'
#' In this program, the strand infomation will not be used.
#' @param ref_path The reference fasta file.
#' @param save_path Where you want save these sequences, only in fastq or fasta format.
#' @param bed_path bed file.
#' @param save_format Fastq or fasta.
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
#' @param ref_path The reference fasta file.
#' @param save_path Where you want save these sequences, only in fastq or fasta format.
#' @param bed_path bed file.
#' @param save_format Fastq or fasta.
DnaSeqCut <- function(atacProc = NULL, ref_path = NULL, save_path = NULL,
                      bed_path = NULL, save_format = NULL){
  tmp <- DNASeqCut$new(atacProc, ref_path, save_path, bed_path, save_format)
  tmp$process()
  return(tmp)
}
