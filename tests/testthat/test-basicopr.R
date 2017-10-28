context("test basic operation of ATACProc object")



test_that("test fragments length distribute",{
    library(magrittr)
    td <- tempdir()
    options(atacConf=setConfigure("tmpdir",td))
    #'
    # Identify adapters
    prefix<-system.file(package="esATAC", "extdata", "uzmg")
    (reads_1 <-file.path(prefix,"m1",dir(file.path(prefix,"m1"))))
    (reads_2 <-file.path(prefix,"m2",dir(file.path(prefix,"m2"))))
    #'
    reads_merged_1 <- file.path(td,"reads1.fastq")
    reads_merged_2 <- file.path(td,"reads2.fastq")
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
    
    expect_equal(getProcName(atacproc),"RemoveAdapter")
   
    (pitems<-getParamItems(atacproc))
    expect_equivalent(pitems,c("fastqInput1","fastqInput2","interleave","fastqOutput1",
                               "fastqOutput2","reportPrefix","paramList","findParamList"))
    expect_false(getParam(atacproc,pitems[3]))
    #'
    #'
    expect_true(isReady(atacproc))
    expect_false(isSingleEnd(atacproc))
    (ritems<-getReportItems(atacproc))
    expect_equivalent(ritems,c("adapters","adapterslist","adapter1","adapter2",
                               "settings","statistics","settingslist","statisticslist","distribution"))
    expect_equal(getReportVal(atacproc,ritems[3]),
                 "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACCTAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAA")
    expect_equal(getReportVal(atacproc,ritems[4]),
                 "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    clearProcCache(atacproc)
    expect_false(file.exists(file.path(td,RemoveAdapterLog)))
    process(atacproc)
    expect_true(file.exists(file.path(td,RemoveAdapterLog)))
    
})