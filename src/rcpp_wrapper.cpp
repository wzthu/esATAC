
#include <Rcpp.h>
#include <string>
#include <iostream>
#include "adapterremoval/adrm_interface.h"
#include "renamer.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "sam2bed.h"
#include "bowtie2/bowtie2_interface.h"
#include "ChrDivi.h"
#include "CutCountPre.h"
#include "CutSiteCount.h"

// [[Rcpp::export]]
int removeAdapter(Rcpp::CharacterVector argvs) {
	int argc=argvs.size();
	char **argv = new char* [argc];
	for(int i=0;i<argc;i++){
		int len = argvs[i].size();
		argv[i] = new char[len+1];
		strcpy(argv[i],(char *)(Rcpp::as<std::string>(argvs[i])).c_str());
    }

    return interface_adapterremoval_main(argc,argv);
}

// [[Rcpp::export]]
int renamer(Rcpp::List argvs) {
        std::string ipath = Rcpp::as<std::string>(argvs["inputFile"]);
        std::string opath = Rcpp::as<std::string>(argvs["outputFile"]);
        std::string filetype = Rcpp::as<std::string>(argvs["fileType"]);
        Renamer rn((char*)ipath.c_str(),(char*)opath.c_str());
        if(filetype=="fq"){
                return rn.renameFastq();
        }else if(filetype=="fa"){
                return rn.renameFasta();
        }

        return 1;
}

// [[Rcpp::export]]
int bowtie2Mapping(Rcpp::CharacterVector argvs) {
  int argc=argvs.size();
  char **argv = new char* [argc];
  for(int i=0;i<argc;i++){
    int len = argvs[i].size();
    argv[i] = new char[len+1];
    strcpy(argv[i],(char *)(Rcpp::as<std::string>(argvs[i])).c_str());

  }
  return interface_bowtie_main(argc, (const char **)argv);
}

// [[Rcpp::export]]
int bowtie2Build(Rcpp::CharacterVector argvs) {
  int argc=argvs.size();
  char **argv = new char* [argc];
  for(int i=0;i<argc;i++){
    int len = argvs[i].size();
    argv[i] = new char[len+1];
    strcpy(argv[i],(char *)(Rcpp::as<std::string>(argvs[i])).c_str());

  }
  return interface_bowtie_build_main(argc, (const char **)argv);
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
int R_sam2bed_wrapper(Rcpp::List argvs)
{
  std::string ipath = Rcpp::as<std::string>(argvs["samfile"]);
  std::string opath = Rcpp::as<std::string>(argvs["bedfile"]);
  SamToBed SB((char*)ipath.c_str(), (char*)opath.c_str());
  return(SB.sam2bed());
}

// [[Rcpp::export]]
int R_sam2bed_merge_wrapper(Rcpp::List argvs,Rcpp::CharacterVector filterList)
{
  std::cout << "11111" << std::endl;
  std::string ipath = Rcpp::as<std::string>(argvs["samfile"]);
  std::string opath = Rcpp::as<std::string>(argvs["bedfile"]);
  int pos_offset = Rcpp::as<int>(argvs["posOffset"]);
  int neg_offset = Rcpp::as<int>(argvs["negOffset"]);
  bool sort = Rcpp::as<bool>(argvs["sort"]);
  bool unique = Rcpp::as<bool>(argvs["unique"]);
  int min_freg_len = Rcpp::as<int>(argvs["minFregLen"]);
  int max_freg_len = Rcpp::as<int>(argvs["maxFregLen"]);
  bool save_ext_len = Rcpp::as<bool>(argvs["saveExtLen"]);

  std::cout << "22222" << std::endl;
  SamToBed SB((char*)ipath.c_str(), (char*)opath.c_str());

  std::cout << "222221" << std::endl;
  int filterSize=filterList.size();


  char **filters = new char* [filterSize];
  if(filterSize == 1){
      int len = filterList[0].size();
      filters[0] = new char[len+1];
      strcpy(filters[0],(char *)(Rcpp::as<std::string>(filterList[0])).c_str());

      if(!strcmp(filters[0],"NULL")){
          delete[] filters[0];
          delete[] filters;

          std::cout.flush();
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

  int t = SB.sam2bed_merge(pos_offset,neg_offset,filters,filterSize,sort,unique, min_freg_len, max_freg_len, save_ext_len);

  if(filters){
      for(int i=0;i<filterSize;i++){
        delete[] filters[i];
      }
      delete[] filters;
  }
  return(t);

}

// [[Rcpp::export]]
int ChrDivi_wrapper(Rcpp::List argvs)
{
  std::cout << "Your input file will be segmented according to the chrmatin name!" << std::endl;
  std::string RIfile = Rcpp::as<std::string>(argvs["readsIfile"]);
  std::string ROfile = Rcpp::as<std::string>(argvs["readsOpath"]);
  std::string Oname = Rcpp::as<std::string>(argvs["name"]);

  // CHR_DIVIDE is the class name
  ChrInfoDivi CHR_DIVIDE(RIfile, ROfile, Oname);
  std::cout << "segmentation finished! Your output file path is:" << std::endl;
  std::cout << ROfile << std::endl;
  return(CHR_DIVIDE.DoDivi());
}

// [[Rcpp::export]]
int CutCountPre_wrapper(Rcpp::List argvs)
{
  std::cout << "Your input file will be segmented according to the chrmatin name!" << std::endl;
  std::string RIfile = Rcpp::as<std::string>(argvs["readsIfile"]);
  std::string ROfile = Rcpp::as<std::string>(argvs["readsOpath"]);

  // CutCount is the class name
  CutCountPre CutCount(RIfile, ROfile);
  std::cout << "segmentation finished! Your output file path is:" << std::endl;
  std::cout << ROfile << std::endl;
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

