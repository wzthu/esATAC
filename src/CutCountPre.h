#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class CutCountPre
{
private:
  string readsIfile;  // your_reads_path/sample_reads.bed
  string readsOpath;  // your_divided_reads_path/

public:
  CutCountPre(string readsIfile, string readsOpath);
  int EXCutCount();
};
