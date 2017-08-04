#include<string>
#include<fstream>
class SamToBed
{
private:
  char * ifilePath;
  char * ofilePath;
  const long MAX_BUFFER_LINE;

public:
  SamToBed(char * ifilePath, char * ofilePath, int memSize=8);
  int sam2bed(int pos_offset,int neg_offset,char ** chrList,int char_filter_size,bool sort,bool unique);
  int sam2bed_merge(int pos_offset=4,int neg_offset=-5,char ** chrList = NULL,int char_filter_size = 0,
                    bool sort = true,bool unique = true,int min_freg_len=0, int max_freg_len=100, bool save_ext_len=false);
private:
    int getReadsLen(char * CIGAR);
};
