#include <queue>
#include "BedLine.h"
#include <functional>
#include <vector>
class SortBed
{
private:
    std::string input_path;
    std::string output_path;
    std::string tmp_prefix;
    std::priority_queue<BedLine*,std::vector<BedLine*>,BedLine> bed_buf;
    int max_line;
    int tmp_count;
    void flush_bed_buf();
public:
    SortBed(const char * outputpath,int maxline);
    SortBed(const char * outputpath,int maxline,const char * tmpprefix);
    SortBed(const char * outputpath,const char * inputpath,int maxline);
    SortBed(const char * outputpath,const char * inputpath,int maxline,const char * tmpprefix);
    ~SortBed(void);
    void insertBedLine(BedLine * bedLine);
    void mergeBed();


};
