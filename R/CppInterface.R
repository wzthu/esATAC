.renamer_call<-function(inputFile,outputFile,fileType="fq",interleave = FALSE){
  argv<-list(inputFile=inputFile,outputFile=outputFile,fileType=fileType,interleave=interleave);
  return(renamer(argv));
}

.sam2bed_merge_call <- function(samfile, bedfile,posOffset,negOffset,sortBed,uniqueBed,filterList,minFregLen,maxFregLen,saveExtLen,downSample = 2e9)
{
    argv <- list(samfile = samfile, bedfile = bedfile ,posOffset = posOffset,negOffset = negOffset,
                 sort = sortBed,unique = uniqueBed, minFregLen = minFregLen, maxFregLen = maxFregLen, saveExtLen = saveExtLen, memSize = 8 ,downSample = downSample,removeXS = TRUE)
    print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    return(R_sam2bed_merge_wrapper(argv,filterList))
}

.sam2bed_call <- function(samfile, bedfile,posOffset,negOffset,sortBed,uniqueBed,filterList,downSample = 2e9){
    argv <- list(samfile = samfile, bedfile = bedfile ,posOffset = posOffset,negOffset = negOffset,
                 sort = sortBed,unique = uniqueBed, memSize = 8 ,downSample = downSample,removeXS = TRUE)
    print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    return(R_sam2bed_wrapper(argv,filterList))
}

.lib_complex_qc_call <- function(bedfile, sortedBed, max_reads){

    argv <- list(bedfile = bedfile ,sortedBed = sortedBed,max_reads = max_reads)
    print(argv)
    rs<-lib_complex_qc(argv)
    if(rs["PBC2"]<0){
        rs["PBC2"]=NA
    }
    print(rs)
    return(rs)
}

# only chr1-chrY will be saved, chrM and others will be removed.
.chr_separate_call <- function(ReadsIfile, ReadsOpath, Name){
  argv <- list(readsIfile = ReadsIfile, readsOpath = ReadsOpath, name = Name)
  ChrDivi_wrapper(argv)
  return(TRUE)
}

.CutSite_call <- function(InputFile, OutputFile){
  argv <- list(readsIfile = InputFile, readsOpath = OutputFile)
  print(argv)
  return(CutCountPre_wrapper(argv))
}


.CutSiteCount <- function(readsfile, motiffile, matrixfile, motif_len, strand_len){
  argv <- list(readsfile = readsfile, motiffile = motiffile, matrixfile = matrixfile,
               motif_len = motif_len, strand_len = strand_len)
  return(CutSiteCount_wrapper(argv))
}




#.identify_adapters_call <- function(inputFile1,inputFile2,findParamList=NULL){
#  argv<-c("AdapterRemoval","--identify-adapters","--file1",
#          inputFile1,"--file2",inputFile2,findParamList);
  #print(argv)
#  removeAdapter(argv);
#  adapter1tb<-read.table(paste(inputFile1,".adapter",sep = ""));
#  adapter2tb<-read.table(paste(inputFile2,".adapter",sep = ""));
#  adapter<-list(adapter1=as.character(adapter1tb[1,1]),adapter2=as.character(adapter2tb[1,1]));
#  return(adapter)
#}
identify_adapters <- function(file1,file2,...,basename = NULL,
                              overwrite = FALSE){
    file1<-trimws(as.character(file1))
    if(!is.null(file2)){
        file2<-trimws(as.character(file2))
    }
    if(!is.null(basename)){
        basename<-trimws(as.character(basename))
    }
    checkFileExist(file1,"file1")
    checkFileExist(file2,"file2")
    checkFileCreatable(paste0(basename,".adapter1"),"file1",overwrite)
    checkFileCreatable(paste0(basename,".adapter2"),"file2",overwrite)
    
    paramArray<-
        checkAddArgus("--identify-adapters|--file1|--file2|--basename",...)
    if(is.null(file2)){
        argvs<-c("AdapterRemoval","--identify-adapters","--file1",
                 file1,"--interleaved",paramArray);
    }else{
        argvs<-c("AdapterRemoval","--identify-adapters","--file1",
                 file1,"--file2",file2,paramArray);
    }
    if(!is.null(basename)){
        argvs<-c(argvs,"--basename",basename)
    }
    
    print(argvs)
    removeAdapter(argvs);
    if(is.null(basename)){
        basename<-"your_output"
    }
    adapter1tb<-readLines(paste0(basename,".adapter1"));
    adapter2tb<-readLines(paste0(basename,".adapter2"));
    return(c(adapter1tb,adapter2tb))
}


remove_adapters <- function(file1,...,adapter1 = NULL,output1 = NULL,
                            file2 = NULL,adapter2 = NULL,output2 = NULL,
                            basename = NULL,interleaved = FALSE,
                            overwrite = FALSE){
    file1<-trimws(as.character(file1))
    if(!is.null(adapter1)){
        adapter1<-trimws(as.character(adapter1))
    }
    if(!is.null(output1)){
        output1<-trimws(as.character(output1))
    }
    if(!is.null(file2)){
        if(interleaved){
            stop("Argumnet `seq2` has to be NULL when interleaved=TRUE")
        }else{
            file2<-trimws(as.character(file2))
            if(length(file1)!=length(file2)){
                stop(paste0("The lengths of arguments `file1` ",
                            "and `file2` should be the same length"))
            }
        }
    }
    if(!is.null(adapter2)){
        adapter2<-trimws(as.character(adapter2))
    }
    if(!is.null(output2)){
        output2<-trimws(as.character(output2))
    }
    if(!is.null(basename)){
        basename<-trimws(as.character(basename))
    }
    
    paramArray<-
        checkAddArgus(paste0("--file1|--file2|--adapter1|--adapter2|",
                             "--output1|--output2|--basename|--interleaved"),...)
    checkFileExist(file1,"file1")
    checkFileExist(file2,"file2")
    checkFileCreatable(output1,"output1",overwrite)
    checkFileCreatable(output2,"output2",overwrite)
    checkPathExist(basename,"basename")
    
    argvs<-c("AdapterRemoval","--file1",file1)
    
    if(!is.null(adapter1)){
        argvs<-c(argvs,"--adapter1",adapter1)
    }
    if(!is.null(output1)){
        argvs<-c(argvs,"--output1",output1)
    }
    if(!is.null(file2)){
        argvs<-c(argvs,"--file2",file2)
    }
    if(!is.null(adapter2)){
        argvs<-c(argvs,"--adapter2",adapter2)
    }
    if(!is.null(output2)){
        argvs<-c(argvs,"--output2",output2)
    }
    if(!is.null(basename)){
        argvs<-c(argvs,"--basename",basename)
    }
    if(interleaved){
        argvs<-c(argvs,"--interleaved")
    }
    argvs<-c(argvs,paramArray)
    invisible(removeAdapter(argvs))
    
    
}
.remove_adapters_call <- function(inputFile1,adapter1,outputFile1,basename,
                                  inputFile2=NULL,adapter2=NULL,outputFile2=NULL,paramlist=NULL){

    argv<-c("AdapterRemoval","--file1",inputFile1, "--adapter1",adapter1,
            "--output1",outputFile1,"--basename", basename );
    if(!is.null(inputFile2)){
        argv<-c(argv,"--file2",inputFile2, "--adapter2",adapter2,"--output2",outputFile2);
    }
    argv<-c(argv,paramlist)
    print(argv)
    return(removeAdapter(argv));
}

# .bowtie2_paired_call <-function(bowtie2Index,samOutput,
#                          fastqInput1, fastqInput2,threads=1){
#     argv<-c("bowtie2-align-s","-x",bowtie2Index,"--no-discordant","--no-unal","--no-mixed","-X","2000",
#             "-p",as.character(threads),"-1",fastqInput1,"-2",fastqInput2,"-S",samOutput)
#     return(bowtie2Mapping(argv))
# }
.bowtie2_call<- function(bowtie2Index,samOutput, fastqInput1, fastqInput2=NULL,paramlist=NULL){
    if(is.null(fastqInput2)){
        argv<-c("bowtie2-align-s",paramlist,"-x",bowtie2Index,"-U",fastqInput1,"-S",samOutput)
    }else{
        argv<-c("bowtie2-align-s",paramlist,"-x",bowtie2Index,"-1",fastqInput1,"-2",fastqInput2,"-S",samOutput)
    }
    return(bowtie2Mapping(argv))
}

.bowtie2_build_call<- function(genomeFastaInput, bt2_base,paramlist=NULL){
    argv<-c("bowtie2-build-s",paramlist,genomeFastaInput,bt2_base)

    return(bowtie2Build(argv))
}

bowtie2_build <- function(references,bt2Index,...,overwrite=FALSE){
    if(R.Version()$os=="mingw32"){
        return("bowtie2 is not support 32bit, please use 64bit R instead")
    }
    references<- trimws(as.character(references))
    bt2Index <- trimws(as.character(bt2Index))
    
    paramArray<-checkAddArgus("noinvalid",...)
    
    checkFileExist(references,"references")
    checkPathExist(bt2Index,"bt2Index")
    checkFileCreatable(paste0(bt2Index,".1.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".2.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".3.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".4.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".rev.1.bt2"),"bt2Index",overwrite)
    checkFileCreatable(paste0(bt2Index,".rev.2.bt2"),"bt2Index",overwrite)
    
    references<-paste0(references,collapse = ",")
    argvs <- c("bowtie2-build-s",paramArray,references,bt2Index)
    
    
    invisible(bowtie2Build(argvs = argvs))
    
}

bowtie2 <- function(bt2Index,samOutput,seq1,...,seq2=NULL,interleaved=FALSE,
                    overwrite=FALSE){
    if(R.Version()$os=="mingw32"){
        return("bowtie2 is not support 32bit, please use 64bit R instead")
    }
    bt2Index <-trimws(as.character(bt2Index))
    samOutput<-trimws(as.character(samOutput))
    
    seq1<-trimws(as.character(seq1))
    
    
    if(!is.null(seq2)){
        seq2<-trimws(as.character(seq2))
        if(length(seq1)!=length(seq2)){
            stop(paste0("The lengths of arguments ",
                        "`seq1` and `seq2` should be the same length"))
        }
    }
    paramArray<-checkAddArgus("-x|--interleaved|-U|-1|-2|-S",...)
    
    
    if(interleaved){
        if(length(seq1)>1){
            stop(paste0("Argumnet `seq1` has to be a SINGLE file",
                        " path rather than a vector of paths"))
        }
        if(!is.null(seq2)){
            stop("Argumnet `seq2` has to be NULL when interleaved=TRUE")
        }
    }
    
    
    checkFileExist(seq1,"seq1")
    checkFileExist(seq2,"seq2")
    checkPathExist(bt2Index,"bt2Index")
    checkFileExist(paste0(bt2Index,".1.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".2.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".3.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".4.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".rev.1.bt2"),"bt2Index")
    checkFileExist(paste0(bt2Index,".rev.2.bt2"),"bt2Index")
    checkFileCreatable(samOutput,"samOutput",overwrite)
    
    argvs = c("-x",bt2Index)
    seq1<-paste0(seq1,collapse = ",")
    if(is.null(seq2)){
        if(interleaved){
            argvs <- c(argvs,"--interleaved",seq1)
        }else{
            argvs <- c(argvs,"-U",seq1)
        }
    }else{
        seq2<-paste0(seq2,collapse = ",")
        argvs <- c(argvs,"-1",seq1,"-2",seq2)
    }
    
    argvs <- c("bowtie2-align-s",paramArray,argvs,"-S",samOutput)
    
    invisible(bowtie2Mapping(argvs = argvs))
    
}

.bedOprUtils_call<-function(ibedfile, obedfile,reportPrefix,
                            mergePair,downSample ,posOffset,negOffset,
                            sortBed,uniqueBed,minFregLen,maxFregLen,filterList,select){
    argv <- list(ibedfile=ibedfile, obedfile=obedfile,reportPrefix=reportPrefix,
                 memSize=8,mergePair=mergePair,downSample = downSample,
                 posOffset=posOffset,negOffset=negOffset,
                 sortBed=sortBed,uniqueBed=uniqueBed,
                 minFregLen=minFregLen,maxFregLen=maxFregLen,
                 filterList=filterList,select=select)
    print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    return(bedOprUtils(argv,filterList))

}


checkFileExist <- function(filePath,argname){
    if(!is.null(filePath)){
        for(i in 1:length(filePath)){
            if(!file.exists(filePath[i])){
                stop(sprintf("For argument `%s`, file does not exist: `%s`",
                             argname,filePath[i]))
            }
        }
    }
}
checkPathExist <- function(filePath,argname){
    if(!is.null(filePath)){
        if(!dir.exists(dirname(filePath))){
            stop(sprintf("For argument `%s`,path does not exist: `%s`",
                         argname,filePath))
        }
    }
}
checkFileCreatable <- function(filePath,argname,overwrite){
    if(!is.null(filePath)){
        if(file.exists(filePath)){
            if(overwrite){
                warning(sprintf(paste0("For argument `%s`, file exist:%s. ",
                                       "It will be overwrited"),
                                argname,filePath));
            }else{
                stop(sprintf(paste0("For argument `%s`,file exist: %s. ",
                                    "Use 'overwrite=TRUE' to overwrite"),
                             argname,filePath));
            }
        }else if(!file.create(filePath)){
            stop(sprintf(paste0("For argument `%s`, cannot create file `%s`.",
                                "\nNo such directory or permission denied."),
                         argname,filePath));
            stop("")
        }else{
            unlink(filePath)
        }
    }
}

checkAddArgus<- function(pattern,...){
    paramlist<-trimws(as.character(list(...)))
    paramArray<-c()
    if(length(paramlist)>0){
        for(i in 1:length(paramlist)){
            paramArray<-c(paramArray,strsplit(paramlist[i],"\\s+")[[1]])
        }
    }
    fixed<-grepl(pattern,paramArray)
    if(sum(fixed)>0){
        invalidp<-paste0(paramArray[fixed],collapse = " ")
        stop(sprintf(paste0("Argument(s) `%s` are invalid for additional ",
                            "argument. Please set them as fixed arguments."),
                     invalidp))
    }
    return(paramArray)
}



