
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

// [[Rcpp::export]]
int removeAdapter(Rcpp::CharacterVector argvs) {
	int argc=argvs.size();
	char **argv = new char* [argc];
	for(int i=0;i<argc;i++){
		int len = argvs[i].size();
		argv[i] = new char[len+1];
		strcpy(argv[i],(char *)(Rcpp::as<std::string>(argvs[i])).c_str());
//(char *)(Rcpp::as<std::string>(argvs[i])).c_str();
		//std::cout<<argv[i]<<std::endl;
	}
	//std::cout<<argc<<std::endl;

/*
        int argc=17;
        char **argv = new char* [argc];
        argv[0]="AdapterRemoval";
        argv[1]="--file1";//"--help";
        argv[2]="./data/SRR891274_1.fastq";
        argv[3]="--file2";
        argv[4]="./data/SRR891274_2.fastq";
        argv[5]="--output1";
        argv[6]="./data/SRR891274_1.fastq.clipped";
        argv[7]="--output2";
        argv[8]="./data/SRR891274_2.fastq.clipped";
        argv[9]="--basename";
        argv[10]="./data/SRR891274.clipped";
        argv[11]="--adapter1";
        argv[12]="CTGTCTCTTATACACATCTCCGAGCCCACGAGACGGACT";
        argv[13]="--adapter2";
        argv[14]="CTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAG";
        argv[15]="--threads";
        argv[16]="28";
*/
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

/*
extern "C" int interface_bwa_main_func(int argc, char *argv[]);

// [[Rcpp::export]]
int bwaMapping(Rcpp::CharacterVector argvs) {
	int argc=argvs.size();
	char **argv = new char* [argc];
	for(int i=0;i<argc;i++){
		int len = argvs[i].size();
		argv[i] = new char[len+1];
		strcpy(argv[i],(char *)(Rcpp::as<std::string>(argvs[i])).c_str());

	}
    return interface_bwa_main_func(argc,argv);
}
*/

// [[Rcpp::export]]
int R_sam2bed_wrapper(Rcpp::List argvs)
{
  std::cout << "11111" << std::endl;
  std::string ipath = Rcpp::as<std::string>(argvs["samfile"]);
  std::string opath = Rcpp::as<std::string>(argvs["bedfile"]);
  int readlen = Rcpp::as<int>(argvs["readlen"]);
  std::cout << "22222" << std::endl;
  SamToBed SB((char*)ipath.c_str(), (char*)opath.c_str(), readlen);

  return(SB.sam2bed());

}
