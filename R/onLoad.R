.onLoad <- function(libname, pkgname) {
    .jpackage(pkgname, lib.loc = libname)
    initPipeFrame(availableGenome = c("hg19",
                                      "hg38",
                                      "mm9",
                                      "mm10",
                                      "danRer10",
                                      "galGal5",
                                      "galGal4",
                                      "rheMac3",
                                      "rheMac8",
                                      "panTro4",
                                      "rn5",
                                      "rn6",
                                      "sacCer2",
                                      "sacCer3",
                                      "susScr3"),
                  defaultJobName = paste0(pkgname,"-pipeline"),
                  defaultCheckAndInstallFunc = checkAndInstall
    )


    addEdges(edges = c(
                  "PeakCallingFseq","FRiPQC",
                  "PeakCallingFseq","RPeakComp",
                  "FindAdapter", "RemoveAdapter"
                  
              ),
             argOrder = 2)
    addEdges(edges = c(
                  "UnzipAndMerge", "Renamer",
                  "UnzipAndMerge", "FastQC",
                  "UnzipAndMerge", "RemoveAdapter",
                  "UnzipAndMerge", "FindAdapter",
                  "Renamer", "RemoveAdapter",
                  "Renamer", "FastQC",
                  "Renamer", "FindAdapter",
                  "SamToBam", "Rsortbam",
                  "SamToBam", "BamToBed",
                  "Rsortbam", "BamToBed",
                  "RemoveAdapter", "Bowtie2Mapping",
                  "Bowtie2Mapping", "SamToBed",
                  "Bowtie2Mapping", "SamToBam",
                  "Bowtie2Mapping", "LibComplexQC",
                  "SamToBed", "PeakCallingFseq",
                  "PeakCallingFseq", "RPeakComp",
                  "SamToBed", "FragLenDistr",
                  "SamToBed", "TSSQC",
                  "SamToBed", "FRiPQC",
                  "SamToBed", "BedToBigWig",
                  #"SamToBed", "PeakQC",
                  "BamToBed", "BedUtils",
                  #"BamToBed", "PeakQC",
                  "BedUtils", "BedToBigWig",
                  "BedUtils", "BedUtils",
                  "BedUtils", "PeakCallingFseq",
                  "BedUtils", "FragLenDistr",
                  "BedUtils", "TSSQC",
                  "BedUtils", "FRiPQC",
                  "BedUtils", "CutSitePre",
                  "BedUtils", "PeakQC",
                  "SamToBed", "BedUtils",
                  "PeakCallingFseq", "PeakQC",
                  "PeakCallingFseq", "RMotifScan",
                  "PeakCallingFseq", "RPeakAnno",
                  "PeakCallingFseq", "RSNPs",
                  "SamToBed", "CutSitePre",
                  "CutSitePre", "CutSiteCountR",
                  "RMotifScan", "CutSiteCountR",
                  "RPeakAnno", "RGo",
                  "RMotifScan", "RSNPs",
                  "RPeakComp", "RMotifScanPair",
                  "RMotifScanPair",  "CutSiteCountR"
              ),
             argOrder = 1)

}
