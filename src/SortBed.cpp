#include "SortBed.h"
#include <queue>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iostream>
#include "RcoutRcerr.h"




SortBed::SortBed(const char * output_path, bool unique,int  max_line)
{
    this->unique = unique;
    tmp_count = 0;
    this->output_path = output_path;
    this->max_line = max_line;


    this->tmp_prefix = output_path;

}

SortBed::SortBed(const char * output_path, bool unique,int  max_line,const char * tmp_prefix)
{
    this->unique = unique;
    tmp_count = 0;
    this->output_path = output_path;
    this->max_line = max_line;


    this->tmp_prefix = tmp_prefix ;

}

SortBed::SortBed(const char * output_path, bool unique, const char * input_path,int max_line)
{
    this->unique = unique;
    tmp_count = 0;
    this->output_path = output_path;
    this->max_line = max_line;


    this->tmp_prefix = output_path;


    this->input_path = input_path;
    std::ifstream ifs(input_path);
    std::string line;


    while(std::getline(ifs,line)){

        insertBedLine(new BedLine(line));

    }
    mergeBed();
}

SortBed::SortBed(const char * output_path, bool unique, const char * input_path,int max_line,const char * tmp_prefix)
{
    this->unique = unique;
    tmp_count = 0;
    this->output_path = output_path;
    this->max_line = max_line;


    this->tmp_prefix = tmp_prefix ;



    this->input_path = input_path;
    std::ifstream ifs(input_path);
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

    bed_buf.push(*bedLine);
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
    cout <<"merge start"<<std::endl;
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
        bed_buf.push(BedLine(line,i));
    }
    const BedLine * bedLine = NULL;

    if(unique){
        BedLine preLine = BedLine("chrStart",-1,-1,"");
        const BedLine * oldLine = & preLine;
        while(!bed_buf.empty()){
            bedLine = &bed_buf.top();
            if((*bedLine)!=(*oldLine)){
                ofs << bedLine->chr << '\t';
                ofs << bedLine->start << '\t';
                ofs << bedLine->end;
                ofs << bedLine->extend << std::endl;
            }
            oldLine = bedLine;
            // std::cout<<bedLine->extend<<std::endl;
            bed_buf.pop();
            if(std::getline(*ifss[bedLine->tag],line)){
                bed_buf.push(BedLine(line,bedLine->tag));
            }else{
                ifss[bedLine->tag]->close();
                delete ifss[bedLine->tag];
                ifss.erase(bedLine->tag);
            }

        }
    }else{
        while(!bed_buf.empty()){
            bedLine =& bed_buf.top();
            ofs << bedLine->chr << '\t';
            ofs << bedLine->start << '\t';
            ofs << bedLine->end;
            ofs << bedLine->extend << std::endl;
            // std::cout<<bedLine->extend<<std::endl;
            bed_buf.pop();
            if(std::getline(*ifss[bedLine->tag],line)){
                bed_buf.push(BedLine(line,bedLine->tag));
            }else{
                ifss[bedLine->tag]->close();
                delete ifss[bedLine->tag];
                ifss.erase(bedLine->tag);
            }

        }
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
    cout <<"merge finish"<<std::endl;
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

    const BedLine * bedLine = NULL;
    for(int i = 0; i < size; i++){
        bedLine = &bed_buf.top();
      //  std::cout<<bedLine->chr<<bedLine->start<<bedLine->end<<bedLine->extend<<std::endl;
    //    std::cout.flush();
        ofs << bedLine->chr << '\t';
        ofs << bedLine->start << '\t';
        ofs << bedLine->end;
        ofs << bedLine->extend << std::endl;
        bed_buf.pop();

    }
    ofs.flush();
    ofs.close();
    cout <<"finish temporary output:"<<tmp_file.c_str()<<std::endl;
}


