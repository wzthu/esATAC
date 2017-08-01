#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "ChrDivi.h"
using namespace std;


ChrInfoDivi::ChrInfoDivi(string readsIfile, string readsOpath, string Outputname)
{
	this -> readsIfile = readsIfile;
	this -> readsOpath = readsOpath;
	this -> Outputname = Outputname;
}

int ChrInfoDivi::DoDivi()
{
	// Input reads bed file
	ifstream readsifile(this -> readsIfile.c_str(), ios::in);

	// Output reads bed file
	ofstream ch1O((this -> readsOpath + this -> Outputname + "_chr1.bed").c_str(), ios::out);
	ofstream ch2O((this -> readsOpath + this -> Outputname + "_chr2.bed").c_str(), ios::out);
	ofstream ch3O((this -> readsOpath + this -> Outputname + "_chr3.bed").c_str(), ios::out);
	ofstream ch4O((this -> readsOpath + this -> Outputname + "_chr4.bed").c_str(), ios::out);
	ofstream ch5O((this -> readsOpath + this -> Outputname + "_chr5.bed").c_str(), ios::out);
	ofstream ch6O((this -> readsOpath + this -> Outputname + "_chr6.bed").c_str(), ios::out);
	ofstream ch7O((this -> readsOpath + this -> Outputname + "_chr7.bed").c_str(), ios::out);
	ofstream ch8O((this -> readsOpath + this -> Outputname + "_chr8.bed").c_str(), ios::out);
	ofstream ch9O((this -> readsOpath + this -> Outputname + "_chr9.bed").c_str(), ios::out);
	ofstream ch10O((this -> readsOpath + this -> Outputname + "_chr10.bed").c_str(), ios::out);
	ofstream ch11O((this -> readsOpath + this -> Outputname + "_chr11.bed").c_str(), ios::out);
	ofstream ch12O((this -> readsOpath + this -> Outputname + "_chr12.bed").c_str(), ios::out);
	ofstream ch13O((this -> readsOpath + this -> Outputname + "_chr13.bed").c_str(), ios::out);
	ofstream ch14O((this -> readsOpath + this -> Outputname + "_chr14.bed").c_str(), ios::out);
	ofstream ch15O((this -> readsOpath + this -> Outputname + "_chr15.bed").c_str(), ios::out);
	ofstream ch16O((this -> readsOpath + this -> Outputname + "_chr16.bed").c_str(), ios::out);
	ofstream ch17O((this -> readsOpath + this -> Outputname + "_chr17.bed").c_str(), ios::out);
	ofstream ch18O((this -> readsOpath + this -> Outputname + "_chr18.bed").c_str(), ios::out);
	ofstream ch19O((this -> readsOpath + this -> Outputname + "_chr19.bed").c_str(), ios::out);
	ofstream ch20O((this -> readsOpath + this -> Outputname + "_chr20.bed").c_str(), ios::out);
	ofstream ch21O((this -> readsOpath + this -> Outputname + "_chr21.bed").c_str(), ios::out);
	ofstream ch22O((this -> readsOpath + this -> Outputname + "_chr22.bed").c_str(), ios::out);
	ofstream chXO((this -> readsOpath + this -> Outputname + "_chrX.bed").c_str(), ios::out);
	ofstream chYO((this -> readsOpath + this -> Outputname + "_chrY.bed").c_str(), ios::out);

	int MAX_LINE_LENGTH = 100000;
	char line[MAX_LINE_LENGTH] ;
	while(readsifile.getline(line, sizeof(line)))
	{
		string tmp_line = line;
		string tmp(strtok(line, "\t"));
		if(tmp == "chr1")
		{
			ch1O << tmp_line + "\n";
		}
		else if(tmp == "chr2")
		{
			ch2O << tmp_line + "\n";
		}
		else if(tmp == "chr3")
		{
			ch3O << tmp_line + "\n";
		}
		else if(tmp == "chr4")
		{
			ch4O << tmp_line + "\n";
		}
		else if(tmp == "chr5")
		{
			ch5O << tmp_line + "\n";
		}
		else if(tmp == "chr6")
		{
			ch6O << tmp_line + "\n";
		}
		else if(tmp == "chr7")
		{
			ch7O << tmp_line + "\n";
		}
		else if(tmp == "chr8")
		{
			ch8O << tmp_line + "\n";
		}
		else if(tmp == "chr9")
		{
			ch9O << tmp_line + "\n";
		}
		else if(tmp == "chr10")
		{
			ch10O << tmp_line + "\n";
		}
		else if(tmp == "chr11")
		{
			ch11O << tmp_line + "\n";
		}
		else if(tmp == "chr12")
		{
			ch12O << tmp_line + "\n";
		}
		else if(tmp == "chr13")
		{
			ch13O << tmp_line + "\n";
		}
		else if(tmp == "chr14")
		{
			ch14O << tmp_line + "\n";
		}
		else if(tmp == "chr15")
		{
			ch15O << tmp_line + "\n";
		}
		else if(tmp == "chr16")
		{
			ch16O << tmp_line + "\n";
		}
		else if(tmp == "chr17")
		{
			ch17O << tmp_line + "\n";
		}
		else if(tmp == "chr18")
		{
			ch18O << tmp_line + "\n";
		}
		else if(tmp == "chr19")
		{
			ch19O << tmp_line + "\n";
		}
		else if(tmp == "chr20")
		{
			ch20O << tmp_line + "\n";
		}
		else if(tmp == "chr21")
		{
			ch21O << tmp_line + "\n";
		}
		else if(tmp == "chr22")
		{
			ch22O << tmp_line + "\n";
		}
		else if(tmp == "chrX")
		{
			chXO << tmp_line + "\n";
		}
		else if(tmp == "chrY")
		{
			chYO << tmp_line + "\n";
		}

	}

	readsifile.close();
	ch1O.close();
	ch2O.close();
	ch3O.close();
	ch4O.close();
	ch5O.close();
	ch6O.close();
	ch7O.close();
	ch8O.close();
	ch9O.close();
	ch10O.close();
	ch11O.close();
	ch12O.close();
	ch13O.close();
	ch14O.close();
	ch15O.close();
	ch16O.close();
	ch17O.close();
	ch18O.close();
	ch19O.close();
	ch20O.close();
	ch21O.close();
	ch22O.close();
	chXO.close();
	chYO.close();

	return 0;
}
