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
    std::string uniqued_path;
    std::priority_queue<BedLine,std::vector<BedLine>,std::less<BedLine>> bed_buf;
    int max_line;
    int tmp_count;
    bool unique;

    int uniquedCounter;
    int saveCounter;

    void flush_bed_buf();
public:
    SortBed(const char * output_path, bool unique, int max_line, const char* uniqued_path = "");/*
    SortBed(const char * output_path, bool unique, int max_line,const char * tmp_prefix, const char* uniqued_path = "");
    SortBed(const char * output_path, bool unique, const char * input_path,int max_line, const char* uniqued_path = "");
    SortBed(const char * output_path, bool unique, const char * input_path,int max_line,const char * tmp_prefix, const char* uniqued_path = "");*/
    ~SortBed(void);
    void insertBedLine(BedLine * bedLine);
    void mergeBed();
    int getUniquedCounter();
    int getSaveCounter();


};
