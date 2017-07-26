#include "SortBed.h"
#include <queue>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iostream>


SortBed::SortBed(const char * outputpath,int  maxline)
{
    tmp_count = 0;
    this->output_path = outputpath;
    this->max_line = maxline;


    this->tmp_prefix = outputpath;

}

SortBed::SortBed(const char * outputpath,int  maxline,const char * tmpprefix)
{
    tmp_count = 0;
    this->output_path = outputpath;
    this->max_line = maxline;


    this->tmp_prefix = tmpprefix ;

}

SortBed::SortBed(const char * outputpath, const char * inputpath,int maxline)
{
    tmp_count = 0;
    this->output_path = outputpath;
    this->max_line = maxline;


    this->tmp_prefix = outputpath;


    this->input_path = inputpath;
    std::ifstream ifs(inputpath);
    std::string line;


    while(std::getline(ifs,line)){

        insertBedLine(new BedLine(line));

    }
    mergeBed();
}

SortBed::SortBed(const char * outputpath, const char * inputpath,int maxline,const char * tmpprefix)
{
    tmp_count = 0;
    this->output_path = outputpath;
    this->max_line = maxline;


    this->tmp_prefix = tmpprefix ;



    this->input_path = inputpath;
    std::ifstream ifs(inputpath);
    std::string line;

   // int count=0;
    while(std::getline(ifs,line)){

        insertBedLine(new BedLine(line));
     //   std::cout<<count++<<std::endl;
     //   std::cout.flush();

    }
    mergeBed();
}

SortBed::~SortBed(void)
{
}


void SortBed::insertBedLine(BedLine * bedLine){

    bed_buf.push(bedLine);
    //std::cout<<bedLine->chr<<bedLine->start<<bedLine->end<<bedLine->extend<<std::endl;
    //std::cout.flush();
    if(bed_buf.size()>=max_line){
        flush_bed_buf();
    }
}




void SortBed::mergeBed(){
    if(bed_buf.size()>0){
        flush_bed_buf();
    }
    std::cout <<"merge start"<<std::endl;
    //std::cout<<"rserserser"<<std::endl;
    if(tmp_count==1){
        std::string str_count;
        std::stringstream ss;
        ss << 0;
        ss>> str_count;
        std::string tmp_file = tmp_prefix + "." + str_count;
        std::remove(output_path.c_str());
        std::rename(tmp_file.c_str(),output_path.c_str());
        return;
    }
    //std::cout<<"rserserser"<<std::endl;
    std::map<int,std::ifstream *>ifss;


    std::ofstream ofs(output_path.c_str());
    std::string line;
    for(int i = 0; i < tmp_count; i++){
        std::string str_count;
        std::stringstream ss;
        ss << i;
        ss>> str_count;
        std::string tmp_file = tmp_prefix + "." + str_count;
        ifss[i] = new std::ifstream(tmp_file.c_str());
        std::getline(*ifss[i],line);
        bed_buf.push(new BedLine(line,i));
    }
    BedLine * bedLine = NULL;
    while(!bed_buf.empty()){
        bedLine = bed_buf.top();
        ofs << bedLine->chr << '\t';
        ofs << bedLine->start << '\t';
        ofs << bedLine->end;
        ofs << bedLine->extend << std::endl;
       // std::cout<<bedLine->extend<<std::endl;
        bed_buf.pop();
        if(std::getline(*ifss[bedLine->tag],line)){
            bed_buf.push(new BedLine(line,bedLine->tag));
        }else{
            ifss[bedLine->tag]->close();
            delete ifss[bedLine->tag];
            ifss.erase(bedLine->tag);
        }
        delete bedLine;
    }
    ofs.flush();
    ofs.close();
    for(int i = 0; i < tmp_count; i++){
        std::string str_count;
        std::stringstream ss;
        ss << i;
        ss>> str_count;
        std::string tmp_file = tmp_prefix + "." + str_count;
        std::remove(tmp_file.c_str());
    }
    std::cout <<"merge finish"<<std::endl;
}



void SortBed::flush_bed_buf(){
    int size = bed_buf.size();
    std::string str_count;
    std::stringstream ss;
    ss << tmp_count;
    ss>> str_count;
    std::string tmp_file = tmp_prefix + "." + str_count;
    tmp_count++;
    std::ofstream ofs(tmp_file.c_str());

    BedLine * bedLine = NULL;
    for(int i = 0; i < size; i++){
        bedLine = bed_buf.top();
      //  std::cout<<bedLine->chr<<bedLine->start<<bedLine->end<<bedLine->extend<<std::endl;
    //    std::cout.flush();
        ofs << bedLine->chr << '\t';
        ofs << bedLine->start << '\t';
        ofs << bedLine->end;
        ofs << bedLine->extend << std::endl;
        bed_buf.pop();
        delete bedLine;
    }
    ofs.flush();
    ofs.close();
    std::cout <<"finish temporary output:"<<tmp_file.c_str()<<std::endl;
}


