#include "LibComplexQC.h"
#include <fstream>
#include <sstream>
#include <map>
#include <cstdlib>
#include <cstring>
#include <iostream>


using std::ifstream;
using std::getline;
using std::stringstream;
using std::string;
using std::map;
using std::atoi;
using std::cout;
using std::endl;



LibComplexQC::LibComplexQC(string unsorted_bed,int max_reads):unsorted_bed(unsorted_bed),max_reads(max_reads){


}

LibComplexQC::LibComplexQC(string sorted_bed):sorted_bed(sorted_bed){
}

LibComplexQC::~LibComplexQC(void)
{
}

void LibComplexQC::calValSorted(){

    ifstream bedin(sorted_bed);
    cout<<sorted_bed<<endl;
    cout.flush();
    string line;


    string chr0="chrStart";
    string start0="-1";
    string end0="-1";
    string strand0="";

    string chr;
    string start;
    string end;
    string strand;

    int one=-1;
    int two=0;
    int total=-1;

    int counter=0;

    int totalcounter=0;
    while(getline(bedin,line)){
        cout<<line<<endl;
        cout.flush();
        counter ++;
        stringstream ss(line);
        ss>>chr>>start>>end;
        ss>>strand>>strand>>strand;
        if(start0!=start||end0!=end||chr0!=chr||strand0!=strand){
            switch(counter)
            {
            case 1:
                one++;
                total++;
                break;
            case 2:
                two++;
                total++;
                break;
            default:
                total++;
            break;
            };
            counter=0;
        }
        totalcounter++;
        chr0=chr;
        start0=start;
        end0=end;
        strand0=strand;
    }
    counter++;
    switch(counter)
    {
    case 1:
        one++;
        total++;
        break;
    case 2:
        two++;
        total++;
        break;
    default:
        total++;
    break;
    };
    NRF = one / (float)totalcounter;
    PBC1 = one / (float)total;
    if(two!=0){
        PBC2 = one / (float)two;
    }else{
        PBC2 = -1;
    }


    exact1 = one;
    exact2 = two;
    exactT = total;
    reads = totalcounter;

}


void  LibComplexQC::calValUnSorted(){

    ifstream bedin(unsorted_bed);
    map<string,int> count;
    string line;



    int totalcounter=0;

    string chr;
    string start;
    string end;
    string strand;
    while(totalcounter<max_reads && getline(bedin,line)){
        stringstream ss(line);
        ss>>chr>>start>>end;
        ss>>strand>>strand>>strand;

        ss.str("");
        ss.clear();
        ss<<chr<<"\t"<<start<<"\t"<<end<<"\t"<<strand;

        count[ss.str()]+=1;
        totalcounter ++;
    }

    int one=0;
    int two=0;
    int total=0;
    for(auto it=count.begin();it!=count.end();it++){
        switch (it->second)
        {
        case 1:
            one++;
            total++;
            break;
        case 2:
            two++;
            total++;
            break;
        default:
            total++;
        break;
        }
        //delete it->first;
        //count.erase(it++);
    }

    NRF = one / (float)totalcounter;
    PBC1 = one / (float)total;
    if(two!=0){
        PBC2 = one / (float)two;
    }else{
        PBC2 = -1;
    }
    exact1 = one;
    exact2 = two;
    exactT = total;
    reads = totalcounter;
}


float LibComplexQC::getNRF(){
    return NRF;
}


float LibComplexQC::getPBC1(){
    return PBC1;
}

float LibComplexQC::getPBC2(){
    return PBC2;
}


int LibComplexQC::getOne(){
    return exact1;

}
int LibComplexQC::getTwo(){
    return exact2;
}
int LibComplexQC::getTotal(){
    return exactT;
}
int LibComplexQC::getReads(){
    return reads;
}
