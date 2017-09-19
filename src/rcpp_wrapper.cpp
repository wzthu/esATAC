

#include <string>
#include <iostream>
#include "renamer.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "sam2bed.h"
#include "ChrDivi.h"
#include "CutCountPre.h"
#include "CutSiteCount.h"
#include "LibComplexQC.h"
#include "BedUtils.h"
#include "RcoutRcerr.h"
#include <Rcpp.h>

// [[Rcpp::export]]
int fastxrenamer(Rcpp::List argvs) {
        std::string ipath = Rcpp::as<std::string>(argvs["inputFile"]);
        std::string opath = Rcpp::as<std::string>(argvs["outputFile"]);
        std::string filetype = Rcpp::as<std::string>(argvs["fileType"]);
        bool interleave= Rcpp::as<bool>(argvs["interleave"]);
        Renamer rn((char*)ipath.c_str(),(char*)opath.c_str());
        if(interleave){
            rn.renameInterleaveFastq();
        }else if(filetype=="fq"){
                rn.renameFastq();
        }else if(filetype=="fa"){
                rn.renameFasta();
        }

        return 0;
}



// [[Rcpp::export]]
void mergeFile(Rcpp::CharacterVector destFile,Rcpp::CharacterVector fileList){
    int fileSize=fileList.size();
    std::string dest = Rcpp::as<std::string>(destFile[0]);
    std::ofstream of_dest(dest.c_str(), std::ios_base::binary);
    for(int i=0;i<fileSize;i++){
        std::string filePath = Rcpp::as<std::string>(fileList[i]);
        std::ifstream if_file(filePath.c_str(), std::ios_base::binary);
        of_dest << if_file.rdbuf();
        if_file.close();
    }
    of_dest.close();
}


// [[Rcpp::export]]
Rcpp::List R_sam2bed_wrapper(Rcpp::List argvs,Rcpp::CharacterVector filterList)
{
    std::string ipath = Rcpp::as<std::string>(argvs["samfile"]);
    std::string opath = Rcpp::as<std::string>(argvs["bedfile"]);
    int pos_offset = Rcpp::as<int>(argvs["posOffset"]);
    int neg_offset = Rcpp::as<int>(argvs["negOffset"]);
    bool sort = Rcpp::as<bool>(argvs["sort"]);
    bool unique = Rcpp::as<bool>(argvs["unique"]);

    int mem_size = Rcpp::as<int>(argvs["memSize"]);
    int down_sample = Rcpp::as<int>(argvs["downSample"]);
    bool remove_xs = Rcpp::as<bool>(argvs["removeXS"]);



    SamToBed SB((char*)ipath.c_str(), (char*)opath.c_str(),mem_size,down_sample,remove_xs);

    int filterSize=filterList.size();


    char **filters = new char* [filterSize];
    if(filterSize == 1){
        int len = filterList[0].size();
        filters[0] = new char[len+1];
        strcpy(filters[0],(char *)(Rcpp::as<std::string>(filterList[0])).c_str());

        if(!strcmp(filters[0],"NULL")){
            delete[] filters[0];
            delete[] filters;

            cout.flush();
            filters = NULL;
            filterSize = 0;


        }
    }else{
        for(int i=0;i<filterSize;i++){
            int len = filterList[i].size();
            filters[i] = new char[len+1];
            strcpy(filters[i],(char *)(Rcpp::as<std::string>(filterList[i])).c_str());
        }
    }
    SB.sam2bed(pos_offset,neg_offset,filters,filterSize,sort,unique);
    if(filters){
        for(int i=0;i<filterSize;i++){
            delete[] filters[i];
        }
        delete[] filters;
    }

    Rcpp::List rs = Rcpp::List::create(Rcpp::Named("total")=SB.getTotalCounter(),
                                       Rcpp::Named("save")=SB.getSaveCounter(),
                                       Rcpp::Named("filted")=SB.getFiltedCounter(),
                                       Rcpp::Named("extlen")=SB.getExtLenCOunter(),
                                       Rcpp::Named("unique")=SB.getUniqueCounter(),
                                       Rcpp::Named("multimap")=SB.getXsCounter());
    
    return rs;


}



// [[Rcpp::export]]
Rcpp::List R_sam2bed_merge_wrapper(Rcpp::List argvs,Rcpp::CharacterVector filterList)
{
    std::string ipath = Rcpp::as<std::string>(argvs["samfile"]);
    std::string opath = Rcpp::as<std::string>(argvs["bedfile"]);
    int pos_offset = Rcpp::as<int>(argvs["posOffset"]);
    int neg_offset = Rcpp::as<int>(argvs["negOffset"]);
    bool sort = Rcpp::as<bool>(argvs["sort"]);
    bool unique = Rcpp::as<bool>(argvs["unique"]);
    int min_freg_len = Rcpp::as<int>(argvs["minFregLen"]);
    int max_freg_len = Rcpp::as<int>(argvs["maxFregLen"]);
    bool save_ext_len = Rcpp::as<bool>(argvs["saveExtLen"]);
    int mem_size = Rcpp::as<int>(argvs["memSize"]);
    int down_sample = Rcpp::as<int>(argvs["downSample"]);
    bool remove_xs = Rcpp::as<bool>(argvs["removeXS"]);

    SamToBed SB((char*)ipath.c_str(), (char*)opath.c_str(),mem_size,down_sample,remove_xs);

    int filterSize=filterList.size();


    char **filters = new char* [filterSize];
    if(filterSize == 1){
        int len = filterList[0].size();
        filters[0] = new char[len+1];
        strcpy(filters[0],(char *)(Rcpp::as<std::string>(filterList[0])).c_str());

        if(!strcmp(filters[0],"NULL")){
            delete[] filters[0];
            delete[] filters;

            cout.flush();
            filters = NULL;
            filterSize = 0;


        }
    }else{
        for(int i=0;i<filterSize;i++){
            int len = filterList[i].size();
            filters[i] = new char[len+1];
            strcpy(filters[i],(char *)(Rcpp::as<std::string>(filterList[i])).c_str());
        }
    }

    SB.sam2bed_merge(pos_offset,neg_offset,filters,filterSize,sort,unique, min_freg_len, max_freg_len, save_ext_len);

    if(filters){
        for(int i=0;i<filterSize;i++){
            delete[] filters[i];
        }
        delete[] filters;
    }

    Rcpp::List rs = Rcpp::List::create(Rcpp::Named("total")=SB.getTotalCounter(),
                                       Rcpp::Named("save")=SB.getSaveCounter(),
                                       Rcpp::Named("filted")=SB.getFiltedCounter(),
                                       Rcpp::Named("extlen")=SB.getExtLenCOunter(),
                                       Rcpp::Named("unique")=SB.getUniqueCounter(),
                                       Rcpp::Named("multimap")=SB.getXsCounter());
    
    return rs;

}

// [[Rcpp::export]]
Rcpp::List bedOprUtils(Rcpp::List argvs,Rcpp::CharacterVector filterList)
{
    std::string ipath = Rcpp::as<std::string>(argvs["ibedfile"]);//
    std::string opath = Rcpp::as<std::string>(argvs["obedfile"]);//
    std::string rpath = Rcpp::as<std::string>(argvs["reportPrefix"]);//
    int mem_size = Rcpp::as<int>(argvs["memSize"]);//
    bool mergePair = Rcpp::as<bool>(argvs["mergePair"]);
    int down_sample = Rcpp::as<int>(argvs["downSample"]);//
    int pos_offset = Rcpp::as<int>(argvs["posOffset"]);
    int neg_offset = Rcpp::as<int>(argvs["negOffset"]);
    bool isSortBed = Rcpp::as<bool>(argvs["sortBed"]);
    bool unique = Rcpp::as<bool>(argvs["uniqueBed"]);
    int min_freg_len = Rcpp::as<int>(argvs["minFregLen"]);
    int max_freg_len = Rcpp::as<int>(argvs["maxFregLen"]);
    //bool report = Rcpp::as<bool>(argvs["report"]);
    bool report = (rpath.size()!=0);
    bool select = Rcpp::as<bool>(argvs["select"]);


    int filterSize=filterList.size();


    char **filters = new char* [filterSize];
    if(filterSize == 1){
        int len = filterList[0].size();
        filters[0] = new char[len+1];
        strcpy(filters[0],(char *)(Rcpp::as<std::string>(filterList[0])).c_str());

        if(!strcmp(filters[0],"NULL")){
            delete[] filters[0];
            delete[] filters;

            cout.flush();
            filters = NULL;
            filterSize = 0;


        }
    }else{
        for(int i=0;i<filterSize;i++){
            int len = filterList[i].size();
            filters[i] = new char[len+1];
            strcpy(filters[i],(char *)(Rcpp::as<std::string>(filterList[i])).c_str());
        }
    }
    BedUtils bedUtils(ipath.c_str(),opath.c_str(),rpath.c_str(), mem_size,
                      mergePair,
                      down_sample,
                      pos_offset,neg_offset,
                      filters ,filterSize ,select,
                      isSortBed ,
                      unique ,
                      min_freg_len, max_freg_len,
                      report);

    bedUtils.bedToBed();

    if(filters){
        for(int i=0;i<filterSize;i++){
            delete[] filters[i];
        }
        delete[] filters;
    }
    
    Rcpp::List rs = Rcpp::List::create(Rcpp::Named("total")=bedUtils.getTotalCounter(),
                                       Rcpp::Named("save")=bedUtils.getSaveCounter(),
                                       Rcpp::Named("filted")=bedUtils.getFiltedCounter(),
                                       Rcpp::Named("extlen")=bedUtils.getExtLenCOunter(),
                                       Rcpp::Named("unique")=bedUtils.getUniqueCounter());
                                     
    
    return rs;
}



// [[Rcpp::export]]
Rcpp::List lib_complex_qc(Rcpp::List argvs)
{

    std::string bedfile = Rcpp::as<std::string>(argvs["bedfile"]);
    bool sortedBed = Rcpp::as<bool>(argvs["sortedBed"]);
    int max_reads = Rcpp::as<int>(argvs["max_reads"]);

    LibComplexQC * qc = NULL;
    cout<<bedfile<<sortedBed<<max_reads<<std::endl;
    cout.flush();
    if(sortedBed){
        cout<<2<<bedfile<<sortedBed<<max_reads<<std::endl;
        cout.flush();
        qc = new LibComplexQC(bedfile);
        qc->calValSorted();
    }else{
        cout<<1<<bedfile<<sortedBed<<max_reads<<std::endl;
        cout.flush();
        qc = new LibComplexQC(bedfile,max_reads);
        qc->calValUnSorted();
    }


    Rcpp::List rs = Rcpp::List::create(Rcpp::Named("NRF")=qc->getNRF(),
        Rcpp::Named("PBC1")=qc->getPBC1(),
        Rcpp::Named("PBC2")=qc->getPBC2(),
        Rcpp::Named("one")=qc->getOne(),
        Rcpp::Named("two")=qc->getTwo(),
        Rcpp::Named("total")=qc->getTotal(),
        Rcpp::Named("reads")=qc->getReads());

    delete qc;

    return rs;

}


// [[Rcpp::export]]
int ChrDivi_wrapper(Rcpp::List argvs)
{
  std::string RIfile = Rcpp::as<std::string>(argvs["readsIfile"]);
  std::string ROfile = Rcpp::as<std::string>(argvs["readsOpath"]);
  std::string Oname = Rcpp::as<std::string>(argvs["name"]);

  // CHR_DIVIDE is the class name
  ChrInfoDivi CHR_DIVIDE(RIfile, ROfile, Oname);
  return(CHR_DIVIDE.DoDivi());
}

// [[Rcpp::export]]
int CutCountPre_wrapper(Rcpp::List argvs)
{
  std::string RIfile = Rcpp::as<std::string>(argvs["readsIfile"]);
  std::string ROfile = Rcpp::as<std::string>(argvs["readsOpath"]);

  // CutCount is the class name
  CutCountPre CutCount(RIfile, ROfile);
  return(CutCount.EXCutCount());
}


// [[Rcpp::export]]
int CutSiteCount_wrapper(Rcpp::List argvs)
{
  std::string readsfile = Rcpp::as<std::string>(argvs["readsfile"]);
  std::string motiffile = Rcpp::as<std::string>(argvs["motiffile"]);
  std::string matrixfile = Rcpp::as<std::string>(argvs["matrixfile"]);
  int motif_len = Rcpp::as<int>(argvs["motif_len"]);
  int strand_len = Rcpp::as<int>(argvs["strand_len"]);

  // CutSite is the class name
  CutSiteCount CutSite(readsfile, motiffile, matrixfile, motif_len, strand_len);
  return(CutSite.DoCutSiteCount());
}

