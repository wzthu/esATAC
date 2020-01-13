#' @importFrom Rbowtie2 bowtie2_build
#' @importFrom Biostrings masks DNAStringSet injectHardMask
#' @importFrom BSgenome available.genomes
#' @importFrom BSgenome installed.genomes
#' @importFrom BiocManager install



#' @importFrom Rcpp  evalCpp
#' @importFrom igraph  graph
#' @importFrom igraph vertex.attributes
#' @importFrom igraph vertex.attributes<-
#' @importFrom igraph are.connected
#' @importFrom Rcpp sourceCpp
#' @importFrom rJava .jpackage
#' @importFrom rJava .jnew
#' @importFrom rJava .jcall
#' @importFrom rtracklayer import
#' @importFrom rtracklayer export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_ribbon
# @importFrom DiagrammeR render_graph
# @importFrom DiagrammeR to_igraph
# @importFrom DiagrammeR select_nodes
# @importFrom DiagrammeR trav_in
# @importFrom DiagrammeR trav_out
# @importFrom DiagrammeR set_node_attrs_ws
# @importFrom DiagrammeR clear_selection
# @importFrom DiagrammeR create_node_df
# @importFrom DiagrammeR create_edge_df
# @importFrom DiagrammeR create_graph
# @importFrom DiagrammeR add_global_graph_attrs
# @importFrom DiagrammeR export_graph
#' @importFrom magrittr %>%
#' @importFrom digest digest
#' @importFrom BSgenome getBSgenome
#' @importFrom Biostrings writeXStringSet
#' @importFrom GenomeInfoDb seqnames
#' @importFrom AnnotationDbi saveDb
#' @importFrom AnnotationDbi loadDb
#' @importFrom GenomicFeatures makeTxDbFromUCSC
#' @importFrom R.utils isGzipped
#' @importFrom R.utils gunzip
#' @importFrom R.utils isBzipped
#' @importFrom R.utils bunzip2
#' @importFrom GenomicRanges coverage
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocGenerics subset
#' @importFrom rmarkdown render
#' @importFrom knitr knit
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics plot
#' @importFrom stats binom.test
#' @importFrom stats fft
#' @importFrom stats runif
#' @importFrom utils download.file
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @import IRanges
#' @import S4Vectors
#' @import tools
#' @import Rsamtools
#' @import methods
# @importFrom markdown markdownToHTML
#' @useDynLib esATAC



checkAndInstallBSgenomeTestgenome <- function(refFilePath){
    genome <- getGenome()
    if(genome == "testgenome"){
        genome <- "hg19"
    }
    checkAndInstallBSgenome(refFilePath, genome)
}

checkAndInstallGenomeSize <- function(refFilePath){
    genome <- getGenome()
    if(genome == "testgenome"){
        genome <- "hg19"
    }
    bsgenomobj <- getBSgenome(genome)
    lens <-seqlengths(bsgenomobj)
    dt<-data.frame(chrom=names(lens), size = as.integer(lens))
    write.table(dt,file=refFilePath,sep='\t',row.names = FALSE, col.names = FALSE,quote = FALSE)
}


checkAndInstallBt2Idx <- function(refFilePath){
    threads <- getThreads()
    message(paste("--threads",as.character(threads)))
    fastaFilePath <- getRefFiles("fasta")
    fileprefix <- substring(refFilePath[1],1,nchar(refFilePath[1]) -6)
    bowtie2_build(fastaFilePath,fileprefix,"-q","--threads",as.character(threads),overwrite=TRUE)
    file.path(getGenome(),"genome")
}


downloadAndGunzip <- function(urlplaceholder,refFilePath){
    genome <- getGenome()
    download.file(url = sprintf(urlplaceholder,genome),
                  destfile = paste0(refFilePath,'.gz'))
    gunzip(paste0(refFilePath,'.gz'),remove = TRUE)
}

downloadDHS <- function(refFilePath){
    downloadAndGunzip("https://wzthu.github.io/esATAC/refdata/%s.DHS.bed.gz", refFilePath) 
}

downloadBlacklist <- function(refFilePath){
    downloadAndGunzip("https://wzthu.github.io/esATAC/refdata/%s.blacklist.bed.gz", refFilePath)        
}

downloadSNP <- function(refFilePath){
    downloadAndGunzip("https://wzthu.github.io/esATAC/refdata/%s.snp.txt.gz", refFilePath)        
}

checkAndInstall <- function(check = TRUE, ...){
    runWithFinishCheck(func = checkAndInstallBSgenome,refName = "bsgenome")
    runWithFinishCheck(func = checkAndInstallGenomeSize, refName = "chromSize",refFilePath = "chrom.size.txt")
    runWithFinishCheck(func = checkAndInstallTxDb,refName = "knownGene")
    runWithFinishCheck(func = checkAndInstallOrgDb,refName = "annoDb")
    runWithFinishCheck(func = checkAndInstallGenomeFa,refName = "fasta", refFilePath = "genome.fa")
    runWithFinishCheck(func = checkAndInstallBt2Idx,refName = "bt2Idx", refFilePath = c("genome.1.bt2",
                                                                                        "genome.2.bt2",
                                                                                        "genome.3.bt2",
                                                                                        "genome.4.bt2",
                                                                                        "genome.rev.1.bt2",
                                                                                        "genome.rev.2.bt2"))
    runWithFinishCheck(func = downloadDHS,refName = "DHS", refFilePath = "DHS.bed", genome = c("hg19","hg38","mm9","mm10"))
    runWithFinishCheck(func = downloadBlacklist,refName = "blacklist", refFilePath = "blacklist.bed", genome = c("hg19","hg38","mm9","mm10"))
    runWithFinishCheck(func = downloadSNP,refName = "SNP", refFilePath = "snp.txt", genome = "hg19")
}




