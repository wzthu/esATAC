#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class CutSiteCount
{
private:
  string readsfile;
  string motiffile;
  string matrixfile;
  int motif_len;
  int strand_len;

public:
  CutSiteCount(string readsfile, string motiffile, string matrixfile, int motif_len, int strand_len);
  int DoCutSiteCount();
};
