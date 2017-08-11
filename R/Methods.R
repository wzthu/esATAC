
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

atacPrintMap <-function(atacProc=NULL,preProc=FALSE,nextProc=TRUE,curProc=TRUE){
    if(is.null(atacProc)){
        GraphMng$new()$printMap()
    }else if(class(atacProc)=="character"){
        GraphMng$new()$printMap(atacProc,preProc,nextProc,curProc)
    }else{
        atacProc$printMap(preProc,nextProc,curProc)
    }
}

atacPrintNextList<-function(atacProc){
    if(class(atacProc)=="character"){
        GraphMng$new()$getNextProcs(atacProc)
    }else{
        GraphMng$new()$getNextProcs(atacProc$getProcName())
    }
}

atacPrintPrevList<-function(atacProc){
    if(class(atacProc)=="character"){
        GraphMng$new()$getPrevProcs(atacProc)
    }else{
        GraphMng$new()$getPrevProcs(atacProc$getProcName())
    }
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
                            fragmentSize=0,featureLength=NULL,bedOutput=NULL,
                            outputFormat=c("bed","wig","npf"), ploidyDir=NULL,
                            wiggleTrackStep=NULL,threshold=NULL,verbose=TRUE,
                            wgThresholdSet=NULL){
    peakcalling <- PeakCallingFseq$new(atacProc,bedInput,background,genomicReadsCount,
                        fragmentSize,featureLength,bedOutput,outputFormat, ploidyDir,
                        wiggleTrackStep,threshold,verbose,wgThresholdSet)
    peakcalling$processing();
    return(peakcalling)
}


atacReadsLenDistr<-function(atacProc,reportPrefix=NULL,bedInput=NULL){
    distr<-ReadsLenDistribute$new(atacProc,reportPrefix,bedInput)
    distr$processing()
    return(distr)
}

atacLibComplexQC<-function(atacProc,reportPrefix=NULL,samInput=NULL,paired = FALSE,subsample=TRUE,subsampleSize=4*10e6){
    libqc<-LibComplexQC$new(atacProc,reportPrefix=reportPrefix,samInput=samInput,paired = paired,
                            subsample=subsample,subsampleSize=subsampleSize,editable=FALSE)
    libqc$processing()
    return(libqc)
}

atacTSSQC<-function(atacProc, txdb.knownGene = NULL,reportPrefix=NULL,bedInput = NULL,fregLenRange=c(0,2000),tssUpdownstream=1000){
    tssQC<-TSSQC$new(atacProc=atacProc, txdb.knownGene=txdb.knownGene,reportPrefix=reportPrefix,bedInput=bedInput,fregLenRange=fregLenRange,tssUpdownstream=tssUpdownstream,editable=FALSE)
    tssQC$processing()
    return(tssQC)
}

atacFripQC<-function(atacProcReads,atacProcPeak,reportPrefix=NULL,readsBedInput=NULL,peakBedInput){
    fripQC<-FRiPQC$new(atacProcReads=atacProcReads,atacProcPeak=atacProcPeak,reportPrefix=reportPrefix,readsBedInput=readsBedInput,peakBedInput=peakBedInput,editable=FALSE)
    fripQC$processing()
    return(fripQC)
}

atacDHSQC<-function(atacProc, reportPrefix=NULL,bedDHS = NULL,bedInput = NULL){
    dhsQC<-DHSQC$new(atacProc, reportPrefix=reportPrefix,bedDHS = bedDHS,bedInput = bedInput)
    dhsQC$processing()
    return(dhsQC)
}

atacBlacklistQC<-function(atacProc, reportPrefix=NULL,bedBlacklist = NULL,bedInput = NULL){
    blacklistQC<-BlacklistQC$new(atacProc, reportPrefix=reportPrefix,bedBlacklist = bedBlacklist,bedInput = bedInput,editable=FALSE)
    blacklistQC$processing()
    return(blacklistQC)
}

atacGenicQC<-function(atacProc, txdb.knownGene = NULL,reportPrefix=NULL,bedInput = NULL,promoterRange=c(-2000,2000)){
    genicQC<-GenicQC(atacProc, txdb.knownGene = txdb.knownGene,reportPrefix=reportPrefix,bedInput = bedInput,promoterRange=promoterRange)
    genicQC$processing()
    return(genicQC)
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
atacSam2Bed <- function(atacProc, merge = TRUE, posOffset = +4, negOffset= -5, chrFilterList= "chrUn.*|chrM|.*random.*",
                        samInput = NULL, bedOutput = NULL, sortBed = TRUE, minFregLen = 0,maxFregLen = 100,
                        saveExtLen = FALSE,uniqueBed = TRUE){
  tmp <- SamToBed$new(atacProc, merge, posOffset, negOffset, chrFilterList, samInput, bedOutput, sortBed, uniqueBed, minFregLen, maxFregLen, saveExtLen)
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
  tmp$processing()
  return(tmp)
}

#' Finding motif binding sites in DNA sequences.
#'
#' @param Seq_file DNA sequence file, fasta or fastq.
#' @param Motif Motif PFM matrix generated by TFBSTools. Must be a SiteSet Class or SiteSetList Class.
#' @param output_file Output file path and name, default: Seq_file_path/output.
MotifScan <- function(atacProc = NULL, Seq_file = NULL,
                      Motif = NULL, output_file = NULL){
  tmp <- RMotifScan$new(atacProc, Seq_file, Motif, output_file)
  tmp$processing()
  return(tmp)
}


#' Using ChIPseeker to annotate peak file.
#'
#' Just for test, other parameters will be add.
#' @param Input peak file.
#' @param Output Output file, default Input_path/output.
#'
PeakAnno <- function(atacProc = NULL, Input = NULL,
                     Output = NULL){
  tmp <- RPeakAnno$new(atacProc, Input, Output)
  tmp$processing()
  return(tmp)
}

#' Using clusterProfiler to do GO analysis(David method).
#'
#'
#'
GODavid <- function(atacProc = NULL, gene = NULL, idType = "ENTREZ_GENE_ID", listType = "Gene",
                    annotation = "GOTERM_BP_FAT",pvalueCutoff = 0.05, pAdjustMethod = "BH",
                    qvalueCutoff = 0.2,species = NA, david.user = NULL, output = NULL){
  tmp <- RGoDavid$new(atacProc, gene, idType, listType,annotation, pvalueCutoff,
                      pAdjustMethod,qvalueCutoff, species, david.user, output)
  tmp$processing()
  return(tmp)
}

#' Using clusterProfiler to do GO analysis.
#'
#'
#'
GOAnalysis <- function(atacProc = NULL, gene = NULL, OrgDb = NULL, keytype = "ENTREZID", ont = "MF",
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = NULL, qvalueCutoff = 0.2,
                       readable = FALSE, pool = FALSE, output = NULL){
  tmp <- RGo$new(atacProc, gene, OrgDb, keytype, ont, pvalueCutoff,
                      pAdjustMethod , universe, qvalueCutoff, readable, pool, output)
  tmp$processing()
  return(tmp)
}


#' Using FunciSNP to do SNP analysis.
#'
#' @param snp.regions.file path: Location of the regions file.
#' @param bio.features.loc path: Location of the biological features folder.
#'
SNPana <- function(atacProc = NULL, snp.regions.file = NULL, bio.features.loc = NULL,
                   built.in.biofeatures = TRUE,
                   par.threads = parallel::detectCores()/2,
                   verbose = par.threads < 2, method.p = "BH",
                   search.window = 200000, output = NULL){
  tmp <- RSNPs$new(atacProc, snp.regions.file, bio.features.loc,
                   built.in.biofeatures, par.threads,
                   verbose, method.p, search.window, output)
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
