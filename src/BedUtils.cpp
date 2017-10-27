#include "BedUtils.h"
#ifdef PLF_SYS_WIN
#include <regex>
#elif PLF_SYS_LINUX
#include <regex.h>
#endif
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include "RcoutRcerr.h"
#include "BedLine.h"
#include "SortBed.h"
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::getline;
#ifdef R_EVN_FLAG

#else
#include<cstdlib>
#endif



BedUtils::BedUtils(const char * inputBedPath,const char * outputBedPath,const char* reportPath, int memSize,
                   bool mergePair,
                   int down_sample,
                   int pos_offset,int neg_offset,
                   char ** chrList ,int char_filter_size ,bool select,
                   bool isSortBed ,
                   bool unique ,
                   int min_freg_len, int max_freg_len,
                   bool report):
    inputBedPath(inputBedPath),
    outputBedPath(outputBedPath),
    reportPath(reportPath),
    downSample(down_sample),
    mergePair(mergePair),
    posOffset(pos_offset),
    negOffset(neg_offset),
    select(select),
    isSortBed(isSortBed),
    unique(unique),
    minFregLen(min_freg_len),
    maxFregLen(max_freg_len),
    isReport(report)
{

    bufferLineCounts = memSize<128?memSize*12000000:150000000;
    string pattern;
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
        isCheckPattern=true;
    }else{
        pattern="";
        isCheckPattern=false;
    }
    this->regexPattern = pattern;
    totalCounter=0;
    saveCounter=0;
    filtedCounter=0;
    extLenCounter=0;
    uniqueCounter=0;
    if(unique){
        this->isSortBed = true;
    }


}





BedUtils::~BedUtils(void)
{
}



void BedUtils::bedToBed(){

#ifdef PLF_SYS_WIN
    std::regex re(regexPattern);
#elif PLF_SYS_LINUX
    const char * patt = regexPattern.c_str();
    regex_t reg;
    const size_t nmatch = 1;
    regmatch_t pm[1];
    regcomp(&reg,patt,REG_EXTENDED|REG_NOSUB);
#endif
    SortBed * sortBed=NULL;
    ofstream *ofBed=NULL;
    ofstream *chrBed=NULL;
    ofstream * extBed=NULL;

    if(isSortBed){
        if(isReport){
            sortBed = new SortBed(outputBedPath,unique,bufferLineCounts,(string(reportPath)+string(".uniq")).c_str());
        }else{
            sortBed = new SortBed(outputBedPath,unique,bufferLineCounts,"");
        }
    }else{
        ofBed = new ofstream(outputBedPath)	;
    }
    if(isReport){
        if(isCheckPattern){
            chrBed = new ofstream((string(reportPath)+string(".chr")).c_str());
        }
        extBed = new ofstream((string(reportPath)+string(".extlen")).c_str());
    }






    vector<int> splist;
    if(downSample<INT_MAX){
        getLineCount();
    }
    bool isDownSample=false;
    if(downSample>0&&downSample<totalCounter){
        permut();
        isDownSample = true;
    }

    BedLine* (BedUtils::*getFregPt)(ifstream&);
    if(mergePair){
        getFregPt=&BedUtils::getFregMerge;
    }else{
        getFregPt=&BedUtils::getFreg;
    }
    BedLine * bedLine =NULL;

    int sampleCounter = -1;
    int len=0;
    bool discard = false;
    ifstream ifBed(inputBedPath);
    while(NULL!=(bedLine=(this->*getFregPt)(ifBed))){
        sampleCounter++;
        //cout<<sampleCounter<<endl;
        //discard down sampling
        if(isDownSample){
            if(idxqueue.empty()){
                break;
            }
            if(sampleCounter!=idxqueue.top()){
                continue;
            }else{
                idxqueue.pop();
            }
        }

        discard = false;
        //discard chr filting
        if(isCheckPattern){
#ifdef PLF_SYS_WIN
            if(std::regex_match(bedLine->chr, re)){
#elif PLF_SYS_LINUX
                if(regexec(&reg,bedLine->chr.c_str(),nmatch,pm,REG_NOTBOL)!=REG_NOMATCH){
#endif
                    if(!select){
                        discard = true;
                        outputBedLine(chrBed,bedLine);
                        filtedCounter++;
                    }
                }else{
                    if(select){
                        discard = true;
                        outputBedLine(chrBed,bedLine);
                        filtedCounter++;
                    }
                }
            }
            //check length
            len = bedLine->end-bedLine->start;
            if(len>maxFregLen||len<minFregLen){
                discard = true;
                outputBedLine(extBed,bedLine);
                extLenCounter++;
            }

            if(!discard){
                if(isSortBed){
                    sortBed->insertBedLine(bedLine);
                }else{
                    outputBedLine(ofBed,bedLine);
                    saveCounter++;
                }
            }

            delete bedLine;
        }


        ifBed.close();

        if(sortBed){
            sortBed->mergeBed();
            saveCounter = sortBed->getSaveCounter();
            if(unique){
                uniqueCounter = sortBed->getUniquedCounter();
            }
            delete sortBed;
        }

        if(ofBed){
            ofBed->close();
            delete ofBed;
        }
        if(chrBed){
            chrBed->close();
            delete chrBed;
        }
        if(extBed){
            extBed->close();
            delete extBed;
        }



    }



    void BedUtils::permut(){

        vector<double> value(totalCounter);
        vector<int> idx(totalCounter);
        for (int i = 0; i < totalCounter; ++i) {
            idx[i] = i;
#ifdef R_EVN_FLAG
            value[i]=(double)(R::runif(0,1));
#else
            value[i] = rand()/(double)INT_MAX;
#endif
        };
        sort(idx.begin(), idx.end(),
             [&value](int i1, int i2) {return value[i1] < value[i2];});
        if(downSample<totalCounter){
            idx.resize(downSample);
        }else{
            downSample=totalCounter;
        }
        for(int i=0;i<downSample;i++){
            idxqueue.push(idx[i]);
        }

    }

    int BedUtils::getLineCount(){
        ifstream ifBed(inputBedPath);
        string linestr;
        totalCounter = 0;

        while(getline(ifBed,linestr)){
            totalCounter++;
        }
        ifBed.close();
        if(this->mergePair){
            return totalCounter/2;
        }else{
            return totalCounter;
        }

    }


    BedLine * BedUtils::getFreg(ifstream& inBed){
        string tmpline;
        if(getline(inBed,tmpline)){
            BedLine *bedLine = new BedLine(tmpline,false) ;
            if(bedLine->strand=='+'){
                bedLine->start += posOffset;
                bedLine->end += posOffset;
            }else{
                bedLine->start += posOffset;
                bedLine->end += posOffset;
            }
            return bedLine;
        }else{
            return NULL;
        }
    }

    BedLine * BedUtils::getFregMerge(ifstream& inBed){
        string tmpline1;
        string tmpline2;
        if(getline(inBed,tmpline1)&&getline(inBed,tmpline2)){
            BedLine * bedline1 = new BedLine(tmpline1,false);
            BedLine * bedline2 = new BedLine(tmpline2,false);
            int start=0;
            int end=0;
            if(bedline1->strand=='+'){
                start = bedline1->start + posOffset;
            }else{
                end = bedline1->end + negOffset;
            }
            if(bedline2->strand=='+'){
                start = bedline2->start + posOffset;
            }else{
                end = bedline2->end + negOffset;
            }
            bedline1->start = start;
            bedline1->end = end;
            delete bedline2;
            return bedline1;
        }else{
            return NULL;
        }

    }



    void BedUtils::outputBedLine(ofstream * ofBed,BedLine*bedLine){
        if(ofBed){
            (*ofBed)<<bedLine->chr<<"\t"<<bedLine->start<<"\t"<<bedLine->end<<bedLine->extend<<endl;
        }
    }


    int BedUtils::getTotalCounter(){
        return totalCounter;
    }
    int BedUtils::getSaveCounter(){
        return saveCounter;
    }
    int BedUtils::getFiltedCounter(){
        return filtedCounter;
    }
    int BedUtils::getExtLenCOunter(){
        return extLenCounter;
    }
    int BedUtils::getUniqueCounter(){
        return uniqueCounter;
    }
