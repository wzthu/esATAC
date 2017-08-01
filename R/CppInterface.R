.renamer_call<-function(inputFile,outputFile,fileType="fq"){
  argv<-list(inputFile=inputFile,outputFile=outputFile,fileType=fileType);
  return(renamer(argv));
}

.sam2bed_merge_call <- function(samfile, bedfile,posOffset,negOffset,sortBed,uniqueBed,filterList,minFregLen,maxFregLen,saveExtLen)
{
    argv <- list(samfile = samfile, bedfile = bedfile ,posOffset = posOffset,negOffset = negOffset,
                 sort = sortBed,unique = uniqueBed, minFregLen = minFregLen, maxFregLen = maxFregLen, saveExtLen = saveExtLen)
    print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    return(R_sam2bed_merge_wrapper(argv,filterList))
}

.sam2bed_call <- function(samfile, bedfile){
  argv <- list(samfile = samfile, bedfile = bedfile)
  print(argv)
  return(R_sam2bed_wrapper(argv))
}

# only chr1-chrY will be saved, chrM and others will be removed.
.chr_separate_call <- function(ReadsIfile, ReadsOpath, Name){
  argv <- list(readsIfile = ReadsIfile, readsOpath = ReadsOpath, name = Name)
  print(argv)
  return(ChrDivi_wrapper(argv))
}

.CutSite_call <- function(InputFile, OutputFile){
  argv <- list(readsIfile = InputFile, readsOpath = OutputFile)
  print(argv)
  return(CutCountPre_wrapper(argv))
}


.CutSiteCount <- function(readsfile, motiffile, matrixfile, motif_len, strand_len){
  argv <- list(readsfile = readsfile, motiffile = motiffile, matrixfile = matrixfile,
               motif_len = motif_len, strand_len = strand_len)
  print(argv)
  return(CutSiteCount_wrapper(argv))
}




.identify_adapters_call <- function(inputFile1,inputFile2,threads=1){
  argv<-c("AdapterRemoval","--identify-adapters","--file1",
          inputFile1,"--file2",inputFile2,"--threads",threads);
  print(argv)
  removeAdapter(argv);
  adapter1tb<-read.table(paste(inputFile1,".adapter",sep = ""));
  adapter2tb<-read.table(paste(inputFile2,".adapter",sep = ""));
  adapter<-list(adapter1=as.character(adapter1tb[1,1]),adapter2=as.character(adapter2tb[1,1]));
  return(adapter)
}

.remove_adapters_call <- function(inputFile1,adapter1,outputFile1,basename,
                                 inputFile2=NULL,adapter2=NULL,outputFile2=NULL,threads=1){

  argv<-c("AdapterRemoval","--file1",inputFile1, "--adapter1",adapter1,
          "--output1",outputFile1,"--basename", basename ,"--threads",threads);
  if(!is.null(inputFile2)){
    argv<-c(argv,"--file2",inputFile2, "--adapter2",adapter2,"--output2",outputFile2);
  }
  return(removeAdapter(argv));
}

.bowtie2_paired_call <-function(bowtie2Index,samOutput,
                         fastqInput1, fastqInput2,threads=1){
    argv<-c("bowtie2-align-s","-x",bowtie2Index,"--no-discordant","--no-unal","--no-mixed","-X","2000",
            "-p",as.character(threads),"-1",fastqInput1,"-2",fastqInput2,"-S",samOutput)
    return(bowtie2Mapping(argv))
}
