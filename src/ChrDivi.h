#include <iostream>
#include <string>
#include <Rcpp.h>
#include <fstream>

using namespace std;

class ChrInfoDivi
{
private:
	string readsIfile;  // your_reads_path/sample_reads.bed
	string readsOpath;  // your_divided_reads_path/
	string Outputname;  // output file name prefix
public:
	ChrInfoDivi(string readsIfile, string readsOpath, string Outputname);
    Rcpp::StringVector DoDivi();
};
