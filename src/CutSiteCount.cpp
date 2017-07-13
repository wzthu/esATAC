#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include "CutSiteCount.h"

using namespace std;

string IntToString(int);
vector<int> ReverseArray(vector<int> orig, unsigned short int b);

CutSiteCount::CutSiteCount(string readsfile, string motiffile, string matrixfile, int motif_len, int strand_len)
{
  this -> readsfile = readsfile;
  this -> motiffile = motiffile;
  this -> matrixfile = matrixfile;
  this -> motif_len = motif_len;
  this -> strand_len = strand_len;
}

int CutSiteCount::DoCutSiteCount()
{
  // These two files must be sorted
  string readsifile = this -> readsfile;
  string motififile = this -> motiffile;
  int motif_len = this -> motif_len;
  int strand_len = this -> strand_len;
  // Output file
  string matrixfile = this -> matrixfile;

  // matrix using to save cut site infomation
  int matrix_len = motif_len + strand_len * 2;
  vector<int> matrix (matrix_len, 0);


  ifstream readsfile(readsifile.c_str(), ios::in);
  ifstream motiffile(motififile.c_str(), ios::in);

  ofstream matrixOutput(matrixfile.c_str(), ios::out);

  char lineR[100000] = {0};
  char lineM[100000] = {0};
  const char *sep = "\t ";

  // read motif file, no value --> exit
  if(!motiffile.getline(lineM, sizeof(lineM)))
  {
    return 0;
  }
  string M_chr(strtok(lineM, sep));
  int M_start = atoi(strtok(NULL, sep));
  int M_end = atoi(strtok(NULL, sep));
  string M_strand(strtok(NULL, sep));
  int M_s = M_start - strand_len;
  int M_e = M_end + strand_len;

  // read reads file, no value --> exit
  if(!readsfile.getline(lineR, sizeof(lineR)))
  {
    return 0;

  }
  int R_start = atoi(strtok(lineR, sep));

  while(1)
  {
    if (R_start < M_s) // in this situation, R file read the next line
    {
      //cout << "step1" << endl;
      // read reads file, no value --> exit
      if(readsfile.getline(lineR, sizeof(lineR)))
      {
        R_start = atoi(strtok(lineR, sep));
      }
      else
      {
        // reads file end, write matrix
        if(M_strand == "p") // motif in 5'-->3'
        {
          for(int i = 0; i < matrix_len; i++) // write matrix to output file
          {
            matrixOutput << IntToString(matrix[i]) + "\t";
          }
          matrixOutput << "\n";
          fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
        }
        else  // motif in 3'-->5'
        {
          matrix = ReverseArray(matrix, matrix_len);
          for(int i = 0; i < matrix_len; i++) // write matrix to output file
          {
            matrixOutput << IntToString(matrix[i]) + "\t";
          }
          matrixOutput << "\n";
          fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
        }

        // for the rest motif, because of reads file end, write 0 matrix
        while(motiffile.getline(lineM, sizeof(lineM)))
        {
          if(M_strand == "p") // motif in 5'-->3'
          {
            for(int i = 0; i < matrix_len; i++) // write matrix to output file
            {
              matrixOutput << IntToString(matrix[i]) + "\t";
            }
            matrixOutput << "\n";
            fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
          }
          else  // motif in 3'-->5'
          {
            matrix = ReverseArray(matrix, matrix_len);
            for(int i = 0; i < matrix_len; i++) // write matrix to output file
            {
              matrixOutput << IntToString(matrix[i]) + "\t";
            }
            matrixOutput << "\n";
            fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
          }
        }
        // two file end
        return 0;
      }
    }
    else if((R_start >= M_s) && (R_start <= M_e))
    {
      //cout << "step2" << endl;
      matrix[R_start - M_s] += 1;// modify matrix value
      if(readsfile.getline(lineR, sizeof(lineR)))
      {
        R_start = atoi(strtok(lineR, sep));
      }
      else
      {
        // reads file end, write matrix
        if(M_strand == "p") // motif in 5'-->3'
        {
          for(int i = 0; i < matrix_len; i++) // write matrix to output file
          {
            matrixOutput << IntToString(matrix[i]) + "\t";
          }
          matrixOutput << "\n";
          fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
        }
        else  // motif in 3'-->5'
        {
          matrix = ReverseArray(matrix, matrix_len);
          for(int i = 0; i < matrix_len; i++) // write matrix to output file
          {
            matrixOutput << IntToString(matrix[i]) + "\t";
          }
          matrixOutput << "\n";
          fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
        }

        // for the rest motif, because of reads file end, write 0 matrix
        while(motiffile.getline(lineM, sizeof(lineM)))
        {
          if(M_strand == "p") // motif in 5'-->3'
          {
            for(int i = 0; i < matrix_len; i++) // write matrix to output file
            {
              matrixOutput << IntToString(matrix[i]) + "\t";
            }
            matrixOutput << "\n";
            fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
          }
          else  // motif in 3'-->5'
          {
            matrix = ReverseArray(matrix, matrix_len);
            for(int i = 0; i < matrix_len; i++) // write matrix to output file
            {
              matrixOutput << IntToString(matrix[i]) + "\t";
            }
            matrixOutput << "\n";
            fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
          }
        }
        // two file end
        return 0;
      }
    }
    else if (R_start > M_e) // a motif site end, write and read next motif
    {
      //cout << "step3" << endl;
      if(M_strand == "p") // motif in 5'-->3'
      {
        for(int i = 0; i < matrix_len; i++) // write matrix to output file
        {
          matrixOutput << IntToString(matrix[i]) + "\t";
        }
        matrixOutput << "\n";
        fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
      }
      else  // motif in 3'-->5'
      {
        matrix = ReverseArray(matrix, matrix_len);
        for(int i = 0; i < matrix_len; i++) // write matrix to output file
        {
          matrixOutput << IntToString(matrix[i]) + "\t";
        }
        matrixOutput << "\n";
        fill(matrix.begin(), matrix.end(), 0); // set matrix to 0
      }
      // read the next motif
      if(motiffile.getline(lineM, sizeof(lineM)))
      {
        M_chr = strtok(lineM, sep);
        M_start = atoi(strtok(NULL, sep));
        M_end = atoi(strtok(NULL, sep));
        M_strand = strtok(NULL, sep);
        M_s = M_start - strand_len;
        M_e = M_end + strand_len;
        if(readsfile.getline(lineR, sizeof(lineR)))
        {
          R_start = atoi(strtok(lineR, sep));
        }
      }
      else
      {
        return 0;
      }
    }
  }

  readsfile.close();
  motiffile.close();
  matrixOutput.close();
  return 0;

}


string IntToString(int i)
{
  string output;
  stringstream ss;
  ss << i;
  ss >> output;
  return output;
}

vector<int> ReverseArray(vector<int> orig, unsigned short int b)
{
  unsigned short int a = 0;
  int swap;
  for (a; a<--b; a++)
  {
    swap = orig[a];
    orig[a] = orig[b];
    orig[b] = swap;
  }
  return orig;
}






