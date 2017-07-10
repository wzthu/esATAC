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

CutSiteCount::CutSiteCount(string ForwReadsFile, string RevReadsFile, string MotifFile, string ForwMatrixFile, string RevMatrixFile, int motif_length, int strand_length)
{
	this -> ForwReadsFile = ForwReadsFile;
	this -> RevReadsFile = RevReadsFile;
	this -> MotifFile = MotifFile;
	this -> ForwMatrixFile = ForwMatrixFile;
	this -> RevMatrixFile = RevMatrixFile;
	this -> motif_length = motif_length;
	this -> strand_length = strand_length;
}


// For forward strand
int CutSiteCount::ForwCutSiteCount()
{
    // These two files must be sorted
    string readsifile = this -> ForwReadsFile;
    string motififile = this -> MotifFile;
    int motif_len = this -> motif_length;
    int strand_len = this -> strand_length;
    // Output file
    string matrixfile = this -> ForwMatrixFile;

    // matrix using to save cut site infomation
    int matrix_len = motif_len + strand_len * 2;
    vector<int> matrix (matrix_len, 0);


    ifstream readsfile(readsifile.c_str(), ios::in);
    ifstream motiffile(motififile.c_str(), ios::in);

    ofstream matrixOutput(matrixfile.c_str(), ios::out);

    int MAX_LINE_LENGTH = 100000;
    char lineR[MAX_LINE_LENGTH] = {0};
    char lineM[MAX_LINE_LENGTH] = {0};
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
    string R_chr(strtok(lineR, sep));
    // This code is for forward strand, so only R_start will be used
    int R_start = atoi(strtok(NULL, sep));
    int R_end = atoi(strtok(NULL, sep));

    while(1)
    {
        if (R_start < M_s) // in this situation, R file read the next line
        {
            //cout << "step1" << endl;
            // read reads file, no value --> exit
            if(readsfile.getline(lineR, sizeof(lineR)))
            {
                R_chr = strtok(lineR, sep);
                R_start = atoi(strtok(NULL, sep));
                R_end = atoi(strtok(NULL, sep));
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
                R_chr = strtok(lineR, sep);
                // This code is for forward strand, so only R_start will be used
                R_start = atoi(strtok(NULL, sep));
                R_end = atoi(strtok(NULL, sep));
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
                // file pointer back
                readsfile.seekg(-10000, ios::cur);
                if (readsfile == NULL)
                {
                    readsfile.clear();
                    readsfile.seekg(0, ios::beg);
                }
                else
                {
                    readsfile.getline(lineR, sizeof(lineR));
                }
                //after file pointer back operation, reads file must contain data
                if(readsfile.getline(lineR, sizeof(lineR)))
                {
                    R_chr = strtok(lineR, sep);
                    // This code is for forward strand, so only R_start will be used
                    R_start = atoi(strtok(NULL, sep));
                    R_end = atoi(strtok(NULL, sep));
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

// For reverse strand
int CutSiteCount::RevCutSiteCount()
{
    // These two files must be sorted
    string readsifile = this -> RevReadsFile;
    string motififile = this -> MotifFile;
    int motif_len = this -> motif_length;
    int strand_len = this -> strand_length;
    // Output file
    string matrixfile = this -> RevMatrixFile;

    // matrix using to save cut site infomation
    int matrix_len = motif_len + strand_len * 2;
    vector<int> matrix (matrix_len, 0);


    ifstream readsfile(readsifile.c_str(), ios::in);
    ifstream motiffile(motififile.c_str(), ios::in);

    ofstream matrixOutput(matrixfile.c_str(), ios::out);

    int MAX_LINE_LENGTH = 100000;
    char lineR[MAX_LINE_LENGTH] = {0};
    char lineM[MAX_LINE_LENGTH] = {0};
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
    string R_chr(strtok(lineR, sep));
    // This code is for reverse strand, so only R_end will be used
    int R_start = atoi(strtok(NULL, sep));
    int R_end = atoi(strtok(NULL, sep));

    while(1)
    {
        if (R_end < M_s) // in this situation, R file read the next line
        {
            //cout << "step1" << endl;
            // read reads file, no value --> exit
            if(readsfile.getline(lineR, sizeof(lineR)))
            {
                R_chr = strtok(lineR, sep);
                R_start = atoi(strtok(NULL, sep));
                R_end = atoi(strtok(NULL, sep));
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
        else if((R_end >= M_s) && (R_end <= M_e))
        {
            //cout << "step2" << endl;
            matrix[R_end - M_s] += 1;// modify matrix value
            if(readsfile.getline(lineR, sizeof(lineR)))
            {
                R_chr = strtok(lineR, sep);
                // This code is for reverse strand, so only R_end will be used
                R_start = atoi(strtok(NULL, sep));
                R_end = atoi(strtok(NULL, sep));
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
        else if (R_end > M_e) // a motif site end, write and read next motif
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
                // file pointer back
                readsfile.seekg(-10000, ios::cur);
                if (readsfile == NULL)
                {
                    readsfile.clear();
                    readsfile.seekg(0, ios::beg);
                }
                else
                {
                    readsfile.getline(lineR, sizeof(lineR));
                }
                //after file pointer back operation, reads file must contain data
                if(readsfile.getline(lineR, sizeof(lineR)))
                {
                    R_chr = strtok(lineR, sep);
                    // This code is for reverse strand, so only R_end will be used
                    R_start = atoi(strtok(NULL, sep));
                    R_end = atoi(strtok(NULL, sep));
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