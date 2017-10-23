#ifndef BEDLINE_H
#define BEDLINE_H
#include<string>
class BedLine
{
public:
    std::string chr;
    int start;
    int end;
    std::string extend;
    char strand;
    int tag;
    BedLine(std::string & linestr, int tag = 0);
    BedLine(std::string & chr, int start, int end, std::string & extend, int tag = 0);
    BedLine(const char* chr, int start, int end, const char* extend, int tag = 0);
    BedLine(std::string & linestr,bool strandAsTag);
    BedLine();
    ~BedLine(void);
    //bool operator () (const BedLine *a,const BedLine *b) const;
    bool operator < (const BedLine &b) const;
    bool operator == (const BedLine & bedLine) const;
    bool operator != (const BedLine & bedLine) const;
};
#endif
