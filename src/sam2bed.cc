#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef PLF_SYS_WIN
#include <regex>
#elif PLF_SYS_LINUX
#include <regex.h>
#endif
#include "sam2bed.h"
#include "BedLine.h"
#include "SortBed.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "RcoutRcerr.h"
using std::endl;
using std::stringstream;
using std::string;
using std::ifstream;
using std::ofstream;
using std::getline;
#define SAM_MAX_LINE_LENGTH 10000
//#define MAX_BUFFER_LINE  100000000
#define BED_LINE_BIT 60
SamToBed::SamToBed(char * ifilePath, char * ofilePath, int memSize,int down_sample,bool rmXS):
    MAX_BUFFER_LINE(memSize<128?memSize*12000000:150000000)// int or long should be considered
{
//    cout<<MAX_BUFFER_LINE<<endl;
//    cout.flush();
    this -> ifilePath = ifilePath;
    this -> ofilePath = ofilePath;
    this -> reads_count = 0;
    this->down_sample = down_sample;

    this->rmXS = rmXS;

    this->totalCounter = 0;
    this->saveCounter = 0;
    this->filtedCounter = 0;
    this->extLenCounter = 0;
    this->uniqueCounter = 0;
	this->xsCounter = 0;
}

int SamToBed::getReadsLen(string & CIGAR){
  int i = 0;
  int size = CIGAR.size();
  int seqlen = 0;
  for(int j = 0; j < size ; j++){
    if(CIGAR[j]<'0'||CIGAR[j]>'9'){
      switch (CIGAR[j]){
      case 'M':
      case 'N':
      case 'D':
      case '=':
      case 'X':
        seqlen += atoi(CIGAR.substr(i,j-i).c_str());
        i = j+1;
      case 'H':
      case 'I':
      case 'P':
      case 'S':
        break;
      default:
        break;
      }
    }
  }
  return seqlen;
}
int SamToBed::sam2bed(int pos_offset,int neg_offset,char ** chrList,int char_filter_size, bool sort,bool unique) {

    ifstream samifs(this -> ifilePath);
    ofstream * bedofs = NULL;
    SortBed* sortBed=NULL;


    int max_buffer_line = MAX_BUFFER_LINE;
    if(unique||sort){
        sortBed = new SortBed(this -> ofilePath,unique,max_buffer_line);
    }else{
        bedofs = new ofstream(this -> ofilePath);
    }
    string line;
    string tok;
    char strand;
    string chr;

    int  chr_start, chr_end, flag, mqs;

    string CIGAR;

    std::string pattern;
    if(char_filter_size>=1){

        pattern=chrList[0];
        if(char_filter_size>1){
            std::stringstream ss;
            ss << pattern ;
            for(int i = 1; i < char_filter_size; i++ ){
                ss << "|" << chrList[i] ;
            }
            ss >> pattern;
        }

    }else{
        pattern="";
    }

    cout<<"pattern:"<<pattern<<endl;
    cout.flush();

#ifdef PLF_SYS_WIN
    std::regex re(pattern);
    std::regex xsre(string(".*\tXS:i:.*"));
#elif PLF_SYS_LINUX
    const char * patt = pattern.c_str();
    regex_t reg;
    const size_t nmatch = 1;
    regmatch_t pm[1];
    regcomp(&reg,patt,REG_EXTENDED|REG_NOSUB);

    //const char * xspatt = ".*\tXS:i:.*";
    regex_t xsreg;
    regmatch_t xspm[1];
    regcomp(&xsreg,".*\tXS:i:.*",REG_EXTENDED|REG_NOSUB);
#endif
    
    string otherFlag="";
    while (getline(samifs,line))
    {

        if (line[0] == '@')
            continue;
        totalCounter ++;
        stringstream ss(line);
        ss >> tok;
        ss >> flag;

        ss >> chr;
		if (chr[0] == '*'){
            continue;
        }
#ifdef PLF_SYS_WIN
        if(std::regex_match(chr, re)){
//            cout<<"matchchr:"<<chr<<endl;
//            cout.flush();
            filtedCounter++;
            continue;
        }
#elif PLF_SYS_LINUX
        if(char_filter_size>=1&&regexec(&reg,chr.c_str(),nmatch,pm,REG_NOTBOL)!=REG_NOMATCH){
            filtedCounter++;
            continue;
        }
#endif
//        cout<<"unmatchchr:"<<chr<<endl;
//        cout.flush();


        ss >> chr_start;
        chr_start -= 1;
        //chr_end = chr_start + readlen;
        ss >> mqs;

        if ((flag & 0x10) == 0)
        {
            strand = '+';
            chr_start += pos_offset;
        }
        else
        {
            strand = '-';
            chr_start += neg_offset;
        }

        ss >> CIGAR;
        reads_count ++;
        chr_end = chr_start + getReadsLen(CIGAR);
        getline(ss,otherFlag);
#ifdef PLF_SYS_WIN
        if(rmXS && std::regex_match(otherFlag, xsre)){
            xsCounter ++;
            continue;
        }
#elif PLF_SYS_LINUX
        if(rmXS && regexec(&xsreg,otherFlag.c_str(),nmatch,pm,REG_NOTBOL)!=REG_NOMATCH){
            xsCounter ++;
            cout << xsCounter;
            continue;
        }
#endif
        if(bedofs){
            if(reads_count>down_sample){
                break;
            }
            (*bedofs) << chr << "\t" << chr_start << "\t" << chr_end << "\t" << tok << "\t" << mqs << "\t" << strand << endl;
            saveCounter ++;
        }else{
            stringstream ss;
            ss << tok <<"\t"<< mqs << "\t" << strand;
            string extend = ss.str();
            sortBed->insertBedLine(new BedLine(chr,chr_start,chr_end, extend));
        }
    }


    samifs.close();
    if(bedofs){
        bedofs->close();
        delete bedofs;
    }else{
        sortBed->mergeBed();
        saveCounter = sortBed->getSaveCounter();
        uniqueCounter = sortBed->getUniquedCounter();
        delete sortBed;
    }

#ifdef PLF_SYS_LINUX
    regfree(&reg);
    regfree(&xsreg);
#endif

    return reads_count;
}

int SamToBed::sam2bed_merge(int pos_offset,int neg_offset,char ** chrList,int char_filter_size,
                            bool sort,bool unique, int min_freg_len, int max_freg_len, bool save_ext_len) {

    ifstream samifs(this -> ifilePath);
    ofstream * bedofs = NULL;
    SortBed* sortBed = NULL;
    SortBed* sortExtBed =NULL;

    int max_buffer_line = MAX_BUFFER_LINE;
    if(save_ext_len){
        max_buffer_line = MAX_BUFFER_LINE / 2;
    }
    if(unique||sort){
        sortBed = new SortBed(this -> ofilePath,unique,max_buffer_line);
    }else{
        bedofs = new ofstream(this -> ofilePath);
    }
    if(save_ext_len){
        sortExtBed = new SortBed((std::string(this -> ofilePath)+std::string(".ext")).c_str(),unique,max_buffer_line);
    }

    string line = "";
    string line1 = "";
    string tok = "start";
    string tok1 = "start1";

    char strand;
  //char strand1;
    string chr;
  //char * chr1;

    int  chr_start, chr_end, flag, mqs;
    int  chr_start1, chr_end1, /*flag1,*/ mqs1;
    int  start, end, length,length1;
    string CIGAR = "";
    string CIGAR1 = "";
   // bool first = true;
    int freg_len;

    std::string pattern;
    if(char_filter_size>=1){

        pattern=chrList[0];
        if(char_filter_size>1){
            std::stringstream ss;
            ss << pattern ;
            for(int i = 1; i < char_filter_size; i++ ){
                ss << "|" << chrList[i] ;
            }
            ss >> pattern;
        }

    }else{
        pattern="";
    }

#ifdef PLF_SYS_WIN
    std::regex re(pattern);
    std::regex xsre(string(".*\tXS:i:.*"));
#elif PLF_SYS_LINUX
    const char * patt = pattern.c_str();
    regex_t reg;
    const size_t nmatch = 1;
    regmatch_t pm[1];
    regcomp(&reg,patt,REG_EXTENDED|REG_NOSUB);

    regex_t xsreg;
    regmatch_t xspm[1];
    regcomp(&xsreg,".*\tXS:i:.*",REG_EXTENDED|REG_NOSUB);
#endif
    while(getline(samifs,line)&&line[0]=='@'){
    }
    string otherFlag="";
    string otherFlag1="";
    while(getline(samifs,line1)){
        stringstream ss(line);
        stringstream ss1(line1);
        ss >> tok;
        ss1 >> tok1;
        if(tok != tok1){
            line = line1;
            continue;
        }
        totalCounter ++;
        ss1 >> flag;
        ss >> flag;
        ss >> chr;
        ss1 >> chr;
        if (chr[0] == '*'){
            continue;
        }
#ifdef PLF_SYS_WIN
        if(std::regex_match(chr, re)){
            filtedCounter++;
            continue;
        }
#elif PLF_SYS_LINUX
        if(char_filter_size>=1&&regexec(&reg,chr.c_str(),nmatch,pm,REG_NOTBOL)!=REG_NOMATCH){
            filtedCounter++;
            continue;
        }
#endif
        ss >> chr_start;
        chr_start -= 1;
        ss1 >> chr_start1;
        chr_start1 -= 1;

        ss >> mqs;
        ss1 >> mqs1;
        ss >> CIGAR;
        ss1 >> CIGAR1;
        getline(ss,otherFlag);
        getline(ss1,otherFlag1);
#ifdef PLF_SYS_WIN
        if(rmXS && (std::regex_match(otherFlag, xsre)
               || std::regex_match(otherFlag1, xsre))){
            xsCounter ++;

            continue;
        }
#elif PLF_SYS_LINUX
        if(rmXS && (regexec(&xsreg,otherFlag.c_str(),nmatch,pm,REG_NOTBOL)!=REG_NOMATCH
               || regexec(&xsreg,otherFlag1.c_str(),nmatch,pm,REG_NOTBOL)!=REG_NOMATCH)){
            xsCounter ++;
            continue;
        }
#endif
        length = getReadsLen(CIGAR);
        length1 = getReadsLen(CIGAR1);
        if ((flag & 0x10) == 0)
        {
            strand = '+';
            chr_start += pos_offset;
            chr_start1 += neg_offset;
            //chr_end = chr_start + length;
            chr_end1 = chr_start1 + length;
            start = chr_start;
            end = chr_end1;
        }
        else
        {
            strand = '-';
            chr_start1 += pos_offset;
            chr_start += neg_offset;
            chr_end = chr_start + length1;
            //chr_end1 = chr_start1 + length;
            start = chr_start1;
            end = chr_end;
        }

        freg_len=end-start;

        if(freg_len>=min_freg_len&&freg_len<=max_freg_len){
            reads_count ++;
            if(bedofs){
                if(reads_count>down_sample){
                    break;
                }
                (*bedofs) << chr << "\t" << start << "\t" << end << "\t" << tok << "\t" << mqs << "\t" << strand << endl;
                saveCounter ++;
            }else{
                stringstream ss;
                ss << tok <<"\t"<< mqs << "\t" << strand;
                string extend = ss.str();
                sortBed->insertBedLine(new BedLine(chr,start,end, extend));
            }
        }else if(save_ext_len){
            stringstream ss;
            ss << tok <<"\t"<< mqs << "\t" << strand;
            string extend = ss.str();
            sortExtBed->insertBedLine(new BedLine(chr,start,end,extend));
        }

    }

    samifs.close();

    if(bedofs){
        bedofs -> close();
        delete bedofs;
    }else{
        sortBed->mergeBed();
        saveCounter = sortBed->getSaveCounter();
        uniqueCounter = sortBed->getUniquedCounter();
        delete sortBed;
    }

    if(sortExtBed){
        sortExtBed->mergeBed();
        extLenCounter = sortExtBed->getSaveCounter();
        delete sortExtBed;
    }
    cout<<"finish"<<std::endl;


#ifdef PLF_SYS_LINUX
    regfree(&reg);
    regfree(&xsreg);
#endif

    return reads_count;
}

int SamToBed::getTotalCounter(){
    return totalCounter;
}
int SamToBed::getSaveCounter(){
    return saveCounter;
}
int SamToBed::getFiltedCounter(){
    return filtedCounter;
}
int SamToBed::getExtLenCOunter(){
    return extLenCounter;
}
int SamToBed::getUniqueCounter(){
    return uniqueCounter;
}
int SamToBed::getXsCounter(){
    return xsCounter;
}
