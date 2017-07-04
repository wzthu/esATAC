
atacPipe <- function(fastqInput1,fastqInput2=NULL,tmpdir=NULL,saveTmpFiles=FALSE,renamer=TRUE,removeAdapter=TRUE){
  atacInputFile(fastqInput1 = fastqInput1, fastqInput2 = fastqInput2)
  gph<-GraphMng$new()
  vlist<-gph$getProcList();
  if(!renamer){
    vlist$renamer<-NULL
  }
  if(!removeAdapter){
    vlist$removeAdapter<-NULL
  }
  return(gph$getSubGraphTopo(as.numeric(vlist),atacInputFile,list(Renamer=Renamer,RemoveAdapter=RemoveAdapter)))
}


atacInputFile <- function(fastqInput1,fastqInput2=NULL){
  seqFile <- SeqFile$new(fastqInput1,fastqInput2);
  seqFile$processing();
  return(seqFile);
}

atacRenamer <- function(atacProc, fastqOutput1=NULL, fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL){
  renamer <- Renamer$new(atacProc, fastqOutput1, fastqOutput2,fastqInput1,fastqInput2,FALSE)
  renamer$processing()
  return(renamer)
}

atacRemoveAdapter <- function(atacProc,adapter1=NULL,adapter2=NULL,fastqOutput1=NULL,reportPrefix=NULL,
                              fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL){
  removeAdapter <- RemoveAdapter$new(atacProc,adapter1,adapter2,fastqOutput1,reportPrefix,
                                     fastqOutput2,fastqInput1, fastqInput2,editable=FALSE)
  removeAdapter$processing()
  return(removeAdapter)
}



#' Mapping reads to the reference using Rbowtie, if output file do not be specified, the output will be named mapping_result.sam
#' @param seq_file A full path of the fa file(containing fa file). For single end, using a list; for paired end, using a list(length = 2).
#' @param ref_file Character scalar. The path to the bowtie index and prefix to align against, in the form </path/to/index>/<prefix>.
#' @param out_file path and output file name, 'E:\\RATAC_test\\output\\output.sam'
#' @param seq_type sequence type, "single", "paired"
#' @export
atacMapping <- function(atacProc = NULL, fileInput = NULL, Reference = NULL, fileOutput = NULL, In_type = NULL){
  tmp <- Mapping$new(atacProc, fileInput, Reference, fileOutput, In_type)
  tmp$processing()
  return(tmp)
}


#' convert sam to bed
#' @param ATAC_obj obj returned from ATAC_mapping
#' @param samfile sam file dir
#' @param bedfile bed file dir
#' @param readlen reads length
#' @export
atacSam2Bed <- function(atacProc = NULL, samfile = NULL, bedfile = NULL, readlen = NULL){
  tmp <- SamToBed$new(atacProc, samfile, bedfile, readlen)
  tmp$processing()
  return(tmp)
}

#' Quality control using Quasr::qQCreport
#' @param input_file a c()
#' @export
atacQCReport <- function(atacProc = NULL, input_file = NULL, output_file = NULL){
  tmp <- QCreporter$new(atacProc, input_file, output_file)
  tmp$processing()
  return(tmp)
}

#' sam2bam using Rbowtie::asBam
#' @export
atacSam2Bam <- function(atacProc = NULL, samfile = NULL, bamfile = NULL){
  tmp <- SamToBam$new(atacProc, samfile, bamfile)
  tmp$processing()
  return(tmp)
}


#' sortbam using Rbowtie::sortBam
#' the output bam file do not have bam header, can not use QCreport function
#' @export
atacBamSort <- function(atacProc = NULL, inputbam = NULL, outputbam = NULL){
  tmp <- Rsortbam$new(atacProc, inputbam, outputbam)
  tmp$processing()
  return(tmp)
}




















atacRenamerResult <- function(fastqOutput1=NULL, fasqOutput2=NULL){
  renamer <- Renamer$new(NULL, fastqOutput1, fasqOutput2,NULL, NULL,editable=TRUE)
  return(renamer)
}

atacRemoveResult <- function(fastqOutput1=NULL,fastqOutput2=NULL){
  removeAdapter <- RemoveAdapter$new(NULL,NULL,NULL,fastqOutput1,NULL,fastqOutput2,NULL, NULL,editable=FALSE)
  return(removeAdapter)
}





#atacRenamerObj<-atacRenamer("A/fastq1","B/fastq2","C/output")

#atacRemoveAdapterObj<-atacRemoveAdapter(atacRenamerObj)

#atacMappingObj<-Mapping(atacRemoveAdapterObj)
