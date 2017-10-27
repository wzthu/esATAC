context("test the subpipeline from raw fastq file to mapped reads in SAM file")

test_that("Unzip and merge FASTQ files, remove adapter and mapping",{
    library(magrittr)
    td <- tempdir()
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
    UnzipAndMergeLog <- dir(td)[grepl(pattern = "UnzipAndMerge.*log",dir(td))]
    RemoveAdapterLog <- dir(td)[grepl(pattern = "RemoveAdapter.*log",dir(td))]
    RenamerLog <- dir(td)[grepl(pattern = "Renamer.*log",dir(td))]

})
         