.renamer_call<-function(inputFile,outputFile,fileType="fq",interleave = FALSE){
    argv<-list(inputFile=inputFile,outputFile=outputFile,fileType=fileType,interleave=interleave);
    invisible(fastxrenamer(argv));
}

.sam2bed_merge_call <- function(samfile, bedfile,posOffset,negOffset,sortBed,uniqueBed,filterList,minFragLen,maxFragLen,saveExtLen,downSample = 2e9)
{
    argv <- list(samfile = samfile, bedfile = bedfile ,posOffset = posOffset,negOffset = negOffset,
                 sort = sortBed,unique = uniqueBed, minFregLen = minFragLen, maxFregLen = maxFragLen, saveExtLen = saveExtLen, memSize = 8 ,downSample = downSample,removeXS = TRUE)
    #print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    invisible(R_sam2bed_merge_wrapper(argv,filterList))
}

.sam2bed_call <- function(samfile, bedfile,posOffset,negOffset,sortBed,uniqueBed,filterList,downSample = 2e9){
    argv <- list(samfile = samfile, bedfile = bedfile ,posOffset = posOffset,negOffset = negOffset,
                 sort = sortBed,unique = uniqueBed, memSize = 8 ,downSample = downSample,removeXS = TRUE)
    #print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    invisible(R_sam2bed_wrapper(argv,filterList))
}

.lib_complex_qc_call <- function(bedfile, sortedBed, max_reads){

    argv <- list(bedfile = bedfile ,sortedBed = sortedBed,max_reads = max_reads)
    #print(argv)
    rs<-lib_complex_qc(argv)
    if(rs["PBC2"]<0){
        rs["PBC2"]=NA
    }
    #print(rs)
    invisible(rs)
}


.chr_separate_call <- function(ReadsIfile, ReadsOpath, Name){
    argv <- list(readsIfile = ReadsIfile, readsOpath = ReadsOpath, name = Name)
    ChrDivi_wrapper(argv)
}

.CutSite_call <- function(InputFile, OutputFile){
    argv <- list(readsIfile = InputFile, readsOpath = OutputFile)
    CutCountPre_wrapper(argv)
}


.CutSiteCount <- function(readsfile, motiffile, matrixfile, motif_len, strand_len){
    argv <- list(readsfile = readsfile, motiffile = motiffile, matrixfile = matrixfile,
                 motif_len = motif_len, strand_len = strand_len)
    CutSiteCount_wrapper(argv)
}




.bedOprUtils_call<-function(ibedfile, obedfile,reportPrefix,
                            mergePair,downSample ,posOffset,negOffset,
                            sortBed,uniqueBed,minFragLen,maxFragLen,filterList,select){
    argv <- list(ibedfile=ibedfile, obedfile=obedfile,reportPrefix=reportPrefix,
                 memSize=8,mergePair=mergePair,downSample = downSample,
                 posOffset=posOffset,negOffset=negOffset,
                 sortBed=sortBed,uniqueBed=uniqueBed,
                 minFregLen=minFragLen,maxFregLen=maxFragLen,
                 filterList=filterList,select=select)
    #print(argv)
    if(is.null(filterList)){
        filterList = c("NULL");
    }
    invisible(bedOprUtils(argv,filterList))

}





