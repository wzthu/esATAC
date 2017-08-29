DNASeqCut <- R6::R6Class(
  classname = "DNASeqCut",
  inherit = BaseProc,

  public = list(
    initialize = function(atacProc, seqRef = NULL, seqOutput = NULL,
                          bedInput = NULL, editable = FALSE){
      super$initialize("DNASeqCut", editable, list(arg1=atacProc))

      # necessary parameters
      if(!is.null(atacProc)){ # get parameter from class PeakCallingFseq
        private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
      }else{
        private$paramlist[["bedInput"]] <- bedInput
      }
      private$paramlist[["seqRef"]] <- seqRef
      private$paramlist[["seqOutput"]] <- seqOutput
      # parameter check
      private$paramValidation()

    } # initialization end

  ), # public end


  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("Bed file:%s", private$paramlist[["bedInput"]]))
      private$writeLog(sprintf("destination:%s", private$paramlist[["seqOutput"]]))
      genome <- GenomeInfoDb::Seqinfo(genome = NA_character_)
      gr_a <- BiocGenerics::unstrand(rtracklayer::import(private$paramlist[["bedInput"]],
                                                         format = "bed", genome = genome))
      DNA_string <- Biostrings::getSeq(private$paramlist[["seqRef"]], gr_a)
      names(DNA_string) <- paste(GenomicRanges::seqnames(gr_a), ":", GenomicRanges::start(gr_a),
                                 "-", GenomicRanges::end(gr_a), sep = "")
      Biostrings::writeXStringSet(DNA_string, filepath = private$paramlist[["seqOutput"]],
                                  append = FALSE, compress = FALSE,
                                  compression_level = NA, format = "fasta")
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
    }, # checkRequireParam end

    checkAllPath = function(){
      private$checkFileExist(private$paramlist[["bedInput"]])
      private$checkPathExist(private$paramlist[["seqOutput"]])
    } # checkAllPath end

  ) # private end

) # class end


#' Cutting sequence according a bed file and save these sequence as fastq or fasta.
#'
#' In this program, the strand infomation will not be used.
#' @param seqRef The reference fasta file.
#' @param seqOutput Where you want save these sequences, only in fastq or fasta format.
#' @param bedInput Input bed file.
DnaSeqCut <- function(atacProc = NULL, seqRef = NULL, seqOutput = NULL,
                      bedInput = NULL){
  tmp <- DNASeqCut$new(atacProc, seqRef, seqOutput, bedInput)
  tmp$process()
  return(tmp)
}
