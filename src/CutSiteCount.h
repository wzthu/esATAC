#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <vector>
#include <stdlib.h>

using namespace std;

class CutSiteCount
{
private:
	string ForwReadsFile;
	string RevReadsFile;
	string MotifFile;
	string ForwMatrixFile;
	string RevMatrixFile;
	int motif_length;
	int strand_length;
	
public:
	CutSiteCount(string ForwReadsFile, string RevReadsFile, string MotifFile, string ForwMatrixFile, string RevMatrixFile, int motif_length, int strand_length);
	int ForwCutSiteCount();
	int RevCutSiteCount();
	
};