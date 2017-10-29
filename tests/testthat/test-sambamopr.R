context("SAM and BAM file operation")



test_that("SAM To BAM ",{
    library(R.utils)
    td <- tempdir()
    dir.create(file.path(td,"sam2bam"))
    td<-file.path(td,"sam2bam")
    options(atacConf=setConfigure("tmpdir",td))
    sam_bz <- system.file("extdata", "Example.sam.bz2", package="esATAC")
    sam_path <- as.vector(bunzip2(filename = sam_bz,
    destname = file.path(td, "Example.sam"),
    ext="bz2", FUN=bzfile, remove = FALSE))
    sam2bam(samInput = sam_path)
    expect_true(file.exists(file.path(td,"Example.SamToBam.bam")))
    expect_true(file.exists(file.path(td,"Example.SamToBam.bam.bai")))
    SamToBamLog <- dir(td)[grepl(pattern = "SamToBam\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,SamToBamLog)))
    
})


test_that("SAM to BED ",{
    library(R.utils)
    library(magrittr)
    td <- tempdir()
    dir.create(file.path(td,"sam2bed"))
    td<-file.path(td,"sam2bed")
    options(atacConf=setConfigure("tmpdir",td))
    #'
    sambzfile <- system.file(package="esATAC", "extdata", "Example.sam.bz2")
    samfile <- file.path(td,"Example.sam")
    bunzip2(sambzfile,destname=samfile,overwrite=TRUE,remove=FALSE)
    atacproc<-samToBed(samInput = samfile) %>%
    atacBedUtils(maxFragLen = 100, chrFilterList = NULL)
    expect_true(file.exists(file.path(td,"Example.BedUtils.bed")))
    expect_true(file.exists(file.path(td,"Example.BedUtils.report")))
    expect_true(file.exists(file.path(td,"Example.sam")))
    expect_true(file.exists(file.path(td,"Example.SamToBed.bed")))
    expect_true(file.exists(file.path(td,"Example.SamToBed.report")))
    BedUtilsLog <- dir(td)[grepl(pattern = "BedUtils\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,BedUtilsLog)))
    SamToBedLog <- dir(td)[grepl(pattern = "SamToBed\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,SamToBedLog)))
    
})


test_that("sort BAM file and to BED",{
    library(Rsamtools)
    library(magrittr)
    td <- tempdir()
    dir.create(file.path(td,"bam2bed"))
    td<-file.path(td,"bam2bed")
    options(atacConf=setConfigure("tmpdir",td))
    ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    bamsort(bamInput = ex1_file) %>% atacBam2Bed
    expect_true(file.exists(file.path(td,"ex1.Rsortbam.bam")))
    expect_true(file.exists(file.path(td,"ex1.Rsortbam.bam.bai")))
    expect_true(file.exists(file.path(td,"ex1.BamToBed.bed")))
    RsortbamLog <- dir(td)[grepl(pattern = "Rsortbam\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,RsortbamLog)))
    BamToBedLog <- dir(td)[grepl(pattern = "BamToBed\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,BamToBedLog)))
    
})



test_that("bed to bigwig",{
    library(R.utils)
    td <- tempdir()
    dir.create(file.path(td,"bed2bw"))
    td<-file.path(td,"bed2bw")
    options(atacConf=setConfigure("tmpdir",td))
    
    bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
    bedfile <- file.path(td,"chr20.50000.bed")
    bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
    
    library(BSgenome.Hsapiens.UCSC.hg19)
    bedToBigWig(bedfile, BSgenome.Hsapiens.UCSC.hg19)
    expect_true(file.exists(file.path(td,"chr20.50000.bed")))
    expect_true(file.exists(file.path(td,"chr20.50000.BedToBigWig.bw")))
    BedToBigWigLog <- dir(td)[grepl(pattern = "BedToBigWig\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,BedToBigWigLog)))
})
