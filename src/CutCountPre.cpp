#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <Rcpp.h>
#include "CutCountPre.h"
#include "RcoutRcerr.h"
using namespace std;

CutCountPre::CutCountPre(string readsIfile, string readsOpath)
{
  this -> readsIfile = readsIfile;
  this -> readsOpath = readsOpath;
}


Rcpp::StringVector CutCountPre::EXCutCount()
{
  // input bed file
  string ipath = this -> readsIfile;
  // output path + prefix
  string opath = this -> readsOpath;

  ifstream readsifile(ipath.c_str(), ios::in);

  // vector to save cut site
  std::vector<int> cutsite;
  Rcpp::StringVector file_name;

  // parameters used in this program
  char line[100000] = {0};
  const char *sep = "\t ";

  // initialization
  if(!readsifile.getline(line, sizeof(line)))
  {
    cout << "WARNING: the input file is empty!" <<endl;
    return file_name;
  }
  string chr(strtok(line, sep));
  int start = atoi(strtok(NULL, sep)) + 1;
  int end = atoi(strtok(NULL, sep));
  string chr_flag;
  chr_flag = chr;
  string outputfile = opath + "_" + chr + ".cs";
  file_name.push_back(outputfile);
  cutsite.push_back(start);
  cutsite.push_back(end);

  while(readsifile.getline(line, sizeof(line)))
  {
    chr = strtok(line, sep);
    start = atoi(strtok(NULL, sep)) + 1;
    end = atoi(strtok(NULL, sep));
    if(chr == chr_flag)
    {
      cutsite.push_back(start);
      cutsite.push_back(end);
    }
    else
    {
      //********************write data to the output file********************
      ofstream readsofile(outputfile.c_str(), ios::out);
      sort(cutsite.begin(), cutsite.end());
      for(unsigned int i = 0; i < cutsite.size(); i++)
      {
        readsofile << cutsite[i] << "\n";
      }
      readsofile.close();
      cutsite.clear();
      //*********************************************************************
      chr_flag = chr;
      outputfile = opath + "_" + chr_flag + ".cs";
      file_name.push_back(outputfile);
      cutsite.push_back(start);
      cutsite.push_back(end);
    }
  }

  //cout << outputfile << endl;
  //********************write the last chromatin data to the output file********************
  ofstream readsofile(outputfile.c_str(), ios::out);
  sort(cutsite.begin(), cutsite.end());
  for(unsigned int i = 0; i < cutsite.size(); i++)
  {
    readsofile << cutsite[i] << "\n";
  }
  readsofile.close();
  //****************************************************************************************
  readsifile.close();

  return file_name;
}
