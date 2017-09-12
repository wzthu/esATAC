#include<string>
#include<fstream>
#include<climits>

class SamToBed
{
private:
  char * ifilePath;
  char * ofilePath;
  const long MAX_BUFFER_LINE;
  unsigned int reads_count;
  unsigned int down_sample;

  bool rmXS;

  int totalCounter;
  int saveCounter;
  int filtedCounter;
  int extLenCounter;
  int uniqueCounter;
  int xsCounter;

public:
  SamToBed(char * ifilePath, char * ofilePath, int memSize=8,int down_sample=UINT_MAX, bool rmXS=true);
  int sam2bed(int pos_offset=4,int neg_offset=-5,char ** chrList= NULL,int char_filter_size= 0,bool sort= true,bool unique= true);
  int sam2bed_merge(int pos_offset=4,int neg_offset=-5,char ** chrList = NULL,int char_filter_size = 0,
                    bool sort = true,bool unique = true,int min_freg_len=0, int max_freg_len=100, bool save_ext_len=false);
  int getTotalCounter();
  int getSaveCounter();
  int getFiltedCounter();
  int getExtLenCOunter();
  int getUniqueCounter();
  int getXsCounter();
private:
    int getReadsLen(std::string & CIGAR);
};
