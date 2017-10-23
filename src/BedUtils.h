#include<climits>
#include<string>
#include<vector>
#include<queue>
#include"BedLine.h"
#include <fstream>
#include <functional>
using std::ofstream;
using std::ifstream;
using std::string;
using std::priority_queue;
using std::vector;
using std::greater;
class BedUtils
{
private:
    const char * inputBedPath;
    const char * outputBedPath;
    const char * reportPath;
    int bufferLineCounts;//
    int downSample;
    bool mergePair;
    int posOffset;
    int negOffset;
    bool select;
    std::string regexPattern;//

    bool isSortBed;
    bool unique;

    int minFregLen;
    int maxFregLen;

    bool isCheckPattern;
    bool isFiltLen;

    bool isReport;


    int totalCounter;
    int saveCounter;
    int filtedCounter;
    int extLenCounter;
    int uniqueCounter;

    void permut();
    int getLineCount();
    BedLine * getFreg(std::ifstream&);
    BedLine * getFregMerge(std::ifstream&);

    void recordFiltedChr(BedLine* bedLine);
    ;
    void recordLenExt(BedLine* bedLine);
    void outputBedLine(ofstream * ofBed,BedLine*bedline);
    std::priority_queue<int,vector<int>,greater<int>> idxqueue;

public:
    BedUtils(const char * inputBedPath,const char * outputBedPath,const char* reportPath, int memSize,
             bool mergePair = true,
             int down_sample=INT_MAX,
             int pos_offset=4,int neg_offset=-5,
             char ** chrList = NULL,int char_filter_size = 0,bool select=false,
             bool isSortBed = true,
             bool unique = true,
             int min_freg_len=0, int max_freg_len=INT_MAX,bool report=true);
    ~BedUtils(void);
    void bedToBed();
    void bedToBedMerge();
    int getTotalCounter();
    int getSaveCounter();
    int getFiltedCounter();
    int getExtLenCOunter();
    int getUniqueCounter();

};

