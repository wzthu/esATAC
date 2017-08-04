#include "BedLine.h"
#include<sstream>
#include<iostream>

BedLine::BedLine(std::string & linestr, int tag){
    std::stringstream ss(linestr);
    ss >> chr;
    ss >> start;
    ss >> end;
    std::getline(ss,extend);
    this->tag = tag;
}

BedLine::BedLine(std::string & chr, int start, int end, std::string & extend, int tag){
    this->chr = chr;
    this->start = start;
    this->end = end;
    this->extend = extend;
    this->tag = tag;
    this->extend = "\t"+this->extend;
}

BedLine::BedLine(const char* chr, int start, int end, const char* extend, int tag){
    this->chr = chr;
    this->start = start;
    this->end = end;
    this->extend = extend;
    this->tag = tag;
    this->extend = "\t"+this->extend;
}


BedLine::~BedLine(void)
{
}

BedLine::BedLine(void)
{
}

/*
bool BedLine::operator () (const BedLine *a,const BedLine *b) const{
    if(a->chr != b->chr){
        return a->chr > b->chr;
    }else if(a->start != b->start){
        return a->start > b->start;
    }else if(a->end != b->end){
        return a->end > b->end;
    }else{
        return false;// should be false when they are equal
    }
}
*/
bool BedLine::operator < (const BedLine &b) const{
  if(chr != b.chr){
    return chr > b.chr;
  }else if(start != b.start){
    return start > b.start;
  }else if(end != b.end){
    return end > b.end;
  }else{
    return false;// should be false when they are equal
  }
}


bool BedLine::operator == (const BedLine & bedLine) const{
    return bedLine.chr == this->chr &&
       bedLine.start == this->start &&
       bedLine.end == this->end;
}

bool BedLine::operator != (const BedLine & bedLine) const{
    return !(bedLine.chr == this->chr &&
       bedLine.start == this->start &&
       bedLine.end == this->end);
}
