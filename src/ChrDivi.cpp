#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "ChrDivi.h"
#include "RcoutRcerr.h"
using namespace std;

ChrInfoDivi::ChrInfoDivi(string readsIfile, string readsOpath, string Outputname)
{
  this -> readsIfile = readsIfile;
  this -> readsOpath = readsOpath;
  this -> Outputname = Outputname;
}

int ChrInfoDivi::DoDivi()
{
  // input bed file
  string ipath = this -> readsIfile;
  // output path + prefix
  string opath = this -> readsOpath;
  string name = this -> Outputname;
  ifstream readsifile(ipath.c_str(), ios::in);

  // vector to save string
  vector<string> content;

  // parameters used in this program
  char line[100000] = {0};
  const char *sep = "\t ";

  // initialization
  if(!readsifile.getline(line, sizeof(line)))
  {
    cout<<"ERROR: the input file is empty!"<<endl;
    return 0;
  }
  string tmp_line = line;
  string chr(strtok(line, sep));
  string chr_flag;
  chr_flag = chr;
  string outputfile = opath + name + "_" + chr + ".bed";
  content.push_back(tmp_line);

  while(readsifile.getline(line, sizeof(line)))
  {
    tmp_line = line;
    chr = strtok(line, sep);
    if(chr == chr_flag)
    {
      content.push_back(tmp_line);
    }
    else
    {
      //********************write data to the output file********************
      ofstream readsofile(outputfile.c_str(), ios::out);
      int ctsize = content.size();
      for(int i = 0; i < ctsize; i++)
      {
        readsofile << content[i] << "\n";
      }
      readsofile.close();
      content.clear();
      //*********************************************************************
      chr_flag = chr;
      outputfile = opath + name + "_" + chr + ".bed";
      content.push_back(tmp_line);
    }
  }

  //********************write the last chromatin data to the output file********************
  ofstream readsofile(outputfile.c_str(), ios::out);
  int ctsize = content.size();
  for(int i = 0; i < ctsize; i++)
  {
    readsofile << content[i] << "\n";
  }
  readsofile.close();
  //****************************************************************************************

  readsifile.close();
  return 0;
}
