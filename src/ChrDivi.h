#include <iostream>
#include <string>
#include <fstream>
using namespace std;
class ChrInfoDivi
{
private:
	string readsIfile;  // your_reads_path/sample_reads.bed
	string readsOpath;  // your_divided_reads_path/
public:
	ChrInfoDivi(string readsIfile, string readsOpath);
	int DoDivi();
};
