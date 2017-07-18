
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

atacMappingBt2 <- function(atacProc,bowtie2Index=NULL,samOutput=NULL,
                           fastqInput1=NULL, fastqInput2=NULL){
    bt2Mapping<-Bowtie2Mapping$new(atacProc,bowtie2Index,samOutput,NULL,
                       fastqInput1, fastqInput2)
    bt2Mapping$processing()
    return(bt2Mapping)
}

atacPeakCalling <- function(atacProc,bedInput=NULL,background=NULL,genomicReadsCount=NULL,
                            fragmentSize=NULL,featureLength=NULL,bedOutput=NULL,
                            outputFormat=c("bed","wig","npf"), ploidyDir=NULL,
                            wiggleTrackStep=NULL,threshold=NULL,verbose=NULL,
                            wgThresholdSet=NULL){
    peakcalling <- PeakCallingFseq$new(atacProc,bedInput,background,genomicReadsCount,
                        fragmentSize,featureLength,bedOutput,outputFormat, ploidyDir,
                        wiggleTrackStep,threshold,verbose,wgThresholdSet)
    peakcalling$processing();
    return(peakcalling)
}

#' Mapping reads to the reference using Rbowtie, if output file do not be specified, the output will be named mapping_result.sam
#' @param seq_file A full path of the fa file(containing fa file). For single end, using a list; for paired end, using a list(length = 2).
#' @param ref_file Character scalar. The path to the bowtie index and prefix to align against, in the form </path/to/index>/<prefix>.
#' @param out_file path and output file name, 'E:\\RATAC_test\\output\\output.sam'
#' @param seq_type sequence type, "single", "paired"
#' @export
atacMappingBt <- function(atacProc = NULL, fileInput = NULL, Reference = NULL, fileOutput = NULL, In_type = NULL){
  tmp <- BowtieMapping$new(atacProc, fileInput, Reference, fileOutput, In_type)
  tmp$processing()
  return(tmp)
}


#' convert sam to bed
#' @param ATAC_obj obj returned from ATAC_mapping
#' @param samfile sam file dir
#' @param bedfile bed file dir
#' @param readlen reads length
#' @export
atacSam2Bed <- function(atacProc, merge = TRUE, posOffset = +4, negOffset= -5, chrFilterList= NULL,
                        samInput = NULL, bedOutput = NULL){
  tmp <- SamToBed$new(atacProc, merge, posOffset, negOffset, chrFilterList, samInput , bedOutput)
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

#' bam2bed using rtracklayer package
#' @export
atacBam2Bed <- function(atacProc = NULL, bamfile = NULL, bedfile = NULL){
  tmp <- BamToBed$new(atacProc, bamfile, bedfile)
  tmp$processing()
  return(tmp)
}

#' separate genome information by chromatin name.
#' @param atacProc Not using now, we will use it in the future.
#' @param ReadsIfile Input bed file path, the first column is chromatin name.
#' @param ReadsOpath The output path, an empty folder would be great, please using "/" even in windows OS.
#' @param prefix the prefix of the output name, format:prefix_chr*.bed, default:output.
#' @param sort TRUE or FALSE, sort every output file by a column.
#' @param sort_col Which column you want to sort, if sort = TRUE, this parameter must be specified.
#' @export
atacChrDivi <- function(atacProc = NULL, ReadsIfile = NULL, ReadsOpath = NULL,
                        prefix = NULL, sort = NULL, sort_col = null){
  tmp <- ChrDivi$new(atacProc, ReadsIfile, ReadsOpath, prefix, sort, sort_col)
  tmp$processing()
  return(tmp)
}

#' extract cut site information
#' @param atacProc Do not use this parameter, we will add nore functions in the future!
#' @param InputFile Input file path, the No.1-3 column is chromatin name, start cut site, end cut site.
#' @param OutputFile The output path, an empty folder would be great.
#' @param prefix Output file name prefix, e.g. prefix_chr*.bed, default "output".
#' @export
atacCutSitePre <- function(atacProc = NULL, InputFile = NULL, OutputFile = NULL, prefix = NULL){
  tmp <- CutSitePre$new(atacProc, InputFile, OutputFile, prefix)
  tmp$processing()
  return(tmp)
}

#' Counting cut site around motif.
#' @param atacProc Do not use this parameter, we will add nore functions in the future!
#' @param CutSiteFile Your cut site infoemation file(from atacCutSitePre function) path with prefix.
#' e.g. "/your_cut_site_information_path/prefix"
#' @param MotifFile Your cut site infoemation file(from atacChrDivi function) path with prefix.
#' e.g. "/your_motif_information_path/prefix"
#' @param MatrixPath The output path with a prefix, an empty folder would be great.
#' e.g. "/where_you_want_to_save_output/prefix"
#' @param motif_length Motif length.
#' @param strand_length How many bp(base pair) do you want to count up/downstream of the motif.
#' @export
atacCutSiteCount <- function(atacProc = NULL, CutSiteFile = NULL, MotifFile = NULL, MatrixPath = NULL,
                             motif_length = NULL, strand_length = NULL, FootPrint = FALSE){
  tmp <- CutSiteCountR$new(atacProc, CutSiteFile, MotifFile, MatrixPath,
                        motif_length, strand_length, FootPrint)
  tmp$processing()
  return(tmp)
}


#' Find overlaps of the two bed files
#' @param atacProc Do not use this parameter, we will add nore functions in the future!
#' @param InputFile Input file path, the No.1-3 column is chromatin name, start cut site, end cut site.
#' @param OutputFile The output path, an empty folder would be great.
#' @param prefix Output file name prefix, e.g. prefix_chr*.bed, default "output".
#' @export
atacBedOverlaps <- function(atacProc = NULL, BedInput1 = NULL, BedInput2 = NULL, Output = NULL,
                            n.col = NULL){
  tmp <- BedOverlaps$new(atacProc, BedInput1, BedInput2, Output, n.col)
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
