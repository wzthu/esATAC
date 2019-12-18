context("test basic operation of ATACProc object")



test_that("test basic operation",{
    td<-tempdir()
    setTmpDir(td)
    
    # Identify adapters
    prefix<-system.file(package="esATAC", "extdata", "uzmg")
    (reads_1 <-file.path(prefix,"m1",dir(file.path(prefix,"m1"))))
    (reads_2 <-file.path(prefix,"m2",dir(file.path(prefix,"m2"))))
    
    reads_merged_1 <- file.path(td,"reads_1.fq")
    reads_merged_2 <- file.path(td,"reads_2.fq")
    atacproc <- atacUnzipAndMerge(fastqInput1 = reads_1,fastqInput2 = reads_2)
    dir(td)
    
    expect_true(file.exists(file.path(td,"esATAC-pipeline","Step_00_pipe_UnzipAndMerge","pipeFrame.obj.log")))
    expect_true(file.exists(file.path(td,"esATAC-pipeline","Step_00_pipe_UnzipAndMerge","reads_1.fq")))
    expect_true(file.exists(file.path(td,"esATAC-pipeline","Step_00_pipe_UnzipAndMerge","reads_2.fq")))
    
    
})