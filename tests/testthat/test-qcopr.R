context("SAM, BAM and BED file quality control operation and fseq")



test_that("test fragments length distribute",{
    library(R.utils)
    td <- tempdir()
    dir.create(file.path(td,"fraglen"))
    td<-file.path(td,"fraglen")
    options(atacConf=setConfigure("tmpdir",td))
    
    bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
    bedfile <- file.path(td,"chr20.50000.bed")
    bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
    fragLenDistr(bedfile)
    #'
    dir(td)
    expect_true(file.exists(file.path(td,"chr20.50000.bed")))
    expect_true(file.exists(file.path(td,"chr20.50000.FragLenDistr.dnagroove.pdf")))
    expect_true(file.exists(file.path(td,"chr20.50000.FragLenDistr.histone.pdf")))
    expect_true(file.exists(file.path(td,"chr20.50000.FragLenDistr.lendistr.pdf")))
    expect_true(file.exists(file.path(td,"chr20.50000.FragLenDistr.lendistr.txt")))
    FragLenDistrLog <- dir(td)[grepl(pattern = "FragLenDistr\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,FragLenDistrLog)))
    
})


test_that("test TSS location enrichment ",{
    
    library(R.utils)
    td <- tempdir()
    dir.create(file.path(td,"tss"))
    td<-file.path(td,"tss")
    options(atacConf=setConfigure("tmpdir",td))
    #'
    bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
    bedfile <- file.path(td,"chr20.50000.bed")
    bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(BSgenome.Hsapiens.UCSC.hg19)
    tssQC(bedfile,TxDb.Hsapiens.UCSC.hg19.knownGene,BSgenome.Hsapiens.UCSC.hg19,fragLenRange=c(180,247),tssUpdownstream = 100)
    
    
    expect_true(file.exists(file.path(td,"chr20.50000.bed")))
    expect_true(file.exists(file.path(td,"chr20.50000.TSSQC.pdf")))
    expect_true(file.exists(file.path(td,"chr20.50000.TSSQC.report.txt")))
    expect_true(file.exists(file.path(td,"chr20.50000.TSSQC.txt")))
    TSSQCLog <- dir(td)[grepl(pattern = "TSSQC\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,TSSQCLog)))
    
})


test_that("test FRiP and peak enrichment quality control and peak calling ",{
    library(R.utils)
    library(magrittr)
    library(rJava)
    td <- tempdir()
    dir.create(file.path(td,"frip"))
    td<-file.path(td,"frip")
    options(atacConf=setConfigure("tmpdir",td))
    #'
    bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
    bedfile <- file.path(td,"chr20.50000.bed")
    bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
    #'
    readsProc<-bedUtils(bedInput = bedfile,maxFragLen = 100, chrFilterList = NULL)
    peaksProc<- readsProc %>% atacPeakCalling
    
    library(BSgenome.Hsapiens.UCSC.hg19)
    
    atacFripQC(readsProc,peaksProc,bsgenome=BSgenome.Hsapiens.UCSC.hg19)
    
    expect_true(file.exists(file.path(td,"chr20.50000.bed")))
    expect_true(file.exists(file.path(td,"chr20.50000.BedUtils.bed")))
    expect_true(file.exists(file.path(td,"chr20.50000.BedUtils.report")))
    expect_true(file.exists(file.path(td,"chr20.50000.PeakCallingFseq.bed")))
    expect_true(file.exists(file.path(td,"chr20.50000.PeakCallingFseq.FRiPQC.report.txt")))
    BedUtilsLog <- dir(td)[grepl(pattern = "BedUtils\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,BedUtilsLog)))
    FRiPQCLog <- dir(td)[grepl(pattern = "FRiPQC\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,FRiPQCLog)))
    PeakCallingFseqLog <- dir(td)[grepl(pattern = "PeakCallingFseq\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,PeakCallingFseqLog)))
    
    blacklistfile <- system.file(package="esATAC", "extdata", "hg19.blacklist.bed")
    peaksProc %>% atacPeakQC(qcbedInput = blacklistfile, bsgenome = BSgenome.Hsapiens.UCSC.hg19)
    
    expect_true(file.exists(file.path(td,"chr20.50000.PeakQC.report.txt")))
    PeakQCLog <- dir(td)[grepl(pattern = "PeakQC\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,PeakQCLog)))
    
})




