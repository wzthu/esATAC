context("SAM and BAM file operate")


test_that("SAM To BAM ",{
    #' library(R.utils)
    #' sam_bz <- system.file("extdata", "Example.sam.bz2", package="ATACpipe")
    #' sam_path <- as.vector(bunzip2(filename = sam_bz,
    #' destname = file.path(getwd(), "Example.sam"),
    #' ext="bz2", FUN=bzfile, remove = FALSE))
    #' sam2bam(samInput = sam_path)
})


test_that("sam to bed ",{
    #' library(R.utils)
    #' library(magrittr)
    #' td <- tempdir()
    #' options(atacConf=setConfigure("tmpdir",td))
    #'
    #' sambzfile <- system.file(package="ATACpipe", "extdata", "Example.sam.bz2")
    #' samfile <- file.path(td,"Example.sam")
    #' bunzip2(sambzfile,destname=samfile,overwrite=TRUE,remove=FALSE)
    #' atacproc<-samToBed(samInput = samfile) %>%
    #' atacBedUtils(maxFregLen = 100, chrFilterList = NULL)
})


test_that("sort BAM file ",{
    #' library(Rsamtools)
    #' ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    #' bamsort(bamInput = ex1_file)
})


test_that("bam to bed ",{
    #' library(Rsamtools)
    #' ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
    #' bam2bed(bamInput = ex1_file)
})


test_that("bed to bigwig",{
    #' library(R.utils)
    #' td <- tempdir()
    #' options(atacConf=setConfigure("tmpdir",td))
    #'
    #' bedbzfile <- system.file(package="ATACpipe", "extdata", "chr20.50000.bed.bz2")
    #' bedfile <- file.path(td,"chr20.50000.bed")
    #' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
    #'
    #' library(BSgenome.Hsapiens.UCSC.hg19)
    #' bedToBigWig(bedfile, BSgenome.Hsapiens.UCSC.hg19)
    #'
    #' dir(td)
})
