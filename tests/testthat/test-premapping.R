context("test the subpipeline from raw fastq file to mapped reads in SAM file")

test_that("Unzip and merge FASTQ files, remove adapter",{
    library(magrittr)
    td <- tempdir()
    dir.create(file.path(td,"umr"))
    td<-file.path(td,"umr")
    options(atacConf=setConfigure("tmpdir",td))
    #'
    # Identify adapters
    prefix<-system.file(package="ATACpipe", "extdata", "uzmg")
    (reads_1 <-file.path(prefix,"m1",dir(file.path(prefix,"m1"))))
    (reads_2 <-file.path(prefix,"m2",dir(file.path(prefix,"m2"))))
    #'
    expect_equal(prod(file.exists(reads_1)),1)
    expect_equal(prod(file.exists(reads_2)),1)
    atacproc <-
    atacUnzipAndMerge(fastqInput1 = reads_1,fastqInput2 = reads_2) %>%
    atacRenamer %>% atacRemoveAdapter
    
    expect_true(file.exists(file.path(td,"reads_1.fastq.1.UnzipAndMerge.fq")))
    expect_true(file.exists(file.path(td,"reads_1.RemoveAdapter.fq")))
    expect_true(file.exists(file.path(td,"reads_1.RemoveAdapter.report.adapter1")))
    expect_true(file.exists(file.path(td,"reads_1.RemoveAdapter.report.adapter2")))
    expect_true(file.exists(file.path(td,"reads_1.RemoveAdapter.report.discarded")))
    expect_true(file.exists(file.path(td,"reads_1.RemoveAdapter.report.settings")))
    expect_true(file.exists(file.path(td,"reads_1.RemoveAdapter.report.singleton.truncated")))
    expect_true(file.exists(file.path(td,"reads_1.Renamer.fq")))
    expect_true(file.exists(file.path(td,"reads_2.fastq.1.UnzipAndMerge.fq")))
    expect_true(file.exists(file.path(td,"reads_2.RemoveAdapter.fq")))
    expect_true(file.exists(file.path(td,"reads_2.Renamer.fq")))
    UnzipAndMergeLog <- dir(td)[grepl(pattern = "UnzipAndMerge\\..*\\.log",dir(td))]
    RemoveAdapterLog <- dir(td)[grepl(pattern = "RemoveAdapter\\..*\\.log",dir(td))]
    RenamerLog <- dir(td)[grepl(pattern = "Renamer\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,UnzipAndMergeLog)))
    expect_true(file.exists(file.path(td,RemoveAdapterLog)))
    expect_true(file.exists(file.path(td,RenamerLog)))

})


test_that("test bowtie2 mapping",{
    td <- tempdir()
    dir.create(file.path(td,"map"))
    td<-file.path(td,"map")
    options(atacConf=setConfigure("tmpdir",td))
    
    ## Building a bowtie2 index
    library("Rbowtie2")
    refs <- dir(system.file(package="ATACpipe", "extdata", "bt2","refs"),
    full=TRUE)
    bowtie2_build(references=refs, bt2Index=file.path(td, "lambda_virus"),
    "--threads 4 --quiet",overwrite=TRUE)
    ## Alignments
    reads_1 <- system.file(package="ATACpipe", "extdata", "bt2", "reads",
    "reads_1.fastq")
    reads_2 <- system.file(package="ATACpipe", "extdata", "bt2", "reads",
    "reads_2.fastq")
    expect_true(file.exists(reads_1))
    expect_true(file.exists(reads_2))
    
    expect_true(file.exists(file.path(td, "lambda_virus.1.bt2")))
    expect_true(file.exists(file.path(td, "lambda_virus.2.bt2")))
    expect_true(file.exists(file.path(td, "lambda_virus.3.bt2")))
    expect_true(file.exists(file.path(td, "lambda_virus.4.bt2")))
    expect_true(file.exists(file.path(td, "lambda_virus.rev.1.bt2")))
    expect_true(file.exists(file.path(td, "lambda_virus.rev.2.bt2")))
    
    bowtie2Mapping(NULL,bt2Idx = file.path(td, "lambda_virus"),
       samOutput = file.path(td, "result.sam"),
       fastqInput1=reads_1,fastqInput2=reads_2,threads=3)
    
    expect_true(file.exists(file.path(td, "result.sam")))
    expect_true(file.exists(file.path(td, "reads_1.Bowtie2Mapping.report")))
    
    Bowtie2MappingLog <- dir(td)[grepl(pattern = "Bowtie2Mapping\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,Bowtie2MappingLog)))
    
})


test_that("test fastqc",{
    td <- tempdir()
    dir.create(file.path(td,"fastqc"))
    td<-file.path(td,"fastqc")
    options(atacConf=setConfigure("tmpdir",td))
    library(R.utils)
    fra_path <- system.file("extdata", "chr20_1.2.fq.bz2", package="ATACpipe")
    fq1 <- as.vector(bunzip2(filename = fra_path,
    destname = file.path(getwd(), "chr20_1.fq"),
    ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
    fra_path <- system.file("extdata", "chr20_2.2.fq.bz2", package="ATACpipe")
    fq2 <- as.vector(bunzip2(filename = fra_path,
    destname = file.path(getwd(), "chr20_2.fq"),
    ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
    qcreport(input_file = c(fq1, fq2))
    expect_true(file.exists(file.path(td, "chr20_1_FastQC.pdf")))
    
    FastQCLog <- dir(td)[grepl(pattern = "FastQC\\..*\\.log",dir(td))]
    expect_true(file.exists(file.path(td,FastQCLog)))
})
         