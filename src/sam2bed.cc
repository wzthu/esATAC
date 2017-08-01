#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef PLF_SYS_WIN
#include <regex>
#elif PLF_SYS_LINUX
#include <regex.h>
#endif
#include "sam2bed.h"
#include "BedLine.h"
#include "SortBed.h"
#include <iostream>
#include <sstream>
#define SAM_MAX_LINE_LENGTH 10000
#define MAX_BUFFER_LINE  100000000
SamToBed::SamToBed(char * ifilePath, char * ofilePath){
  this -> ifilePath = ifilePath;
  this -> ofilePath = ofilePath;
}

int SamToBed::getReadsLen(char * CIGAR){
	int seqlen = 0;
		for(char *s = CIGAR; *s != '\0'; ){
			if(*s>'9'||*s<'0'){
				*s='\0';
				seqlen += atoi(CIGAR);
				CIGAR = ++s;
			}else{
			    s++;
			}
	    }
	return seqlen;
}
int SamToBed::sam2bed() {

  FILE *fp, *fp_out;

  fp = fopen(this -> ifilePath, "r");
  fp_out = fopen(this -> ofilePath, "w");

  char * line = NULL;
  char * tok;
  char strand;
  char * chr;

  int  chr_start, chr_end, flag, mqs;


  line = (char*)calloc(SAM_MAX_LINE_LENGTH, 1);
  char *CIGAR = NULL;

  while (fgets(line, SAM_MAX_LINE_LENGTH, fp))
  {
    if (line[0] == '@')
      continue;

    tok = strtok(line, "\t");
    flag = atoi(strtok(NULL, "\t"));

    chr = strtok(NULL, "\t");
    if (chr[0] != '*')
    {
      chr_start = atoi(strtok(NULL, "\t")) - 1;
      //chr_end = chr_start + readlen;
      mqs = atoi(strtok(NULL, "\t"));

      if ((flag & 0x10) == 0)
      {
        strand = '+';
      }
      else
      {
        strand = '-';
      }

	  CIGAR = strtok(NULL, "\t");
	  chr_end = chr_start + getReadsLen(CIGAR);
      fprintf(fp_out, "%s\t%d\t%d\t%s\t%d\t%c\n", chr, chr_start, chr_end, tok, mqs, strand);
    }
  }

  if (line)
    free(line);

  fclose(fp);
  fclose(fp_out);

  return 0;
}

int SamToBed::sam2bed_merge(int pos_offset,int neg_offset,char ** chrList,int char_filter_size,
                            bool sort,bool unique, int min_freg_len, int max_freg_len, bool save_ext_len) {

  FILE *fp, *fp_out=NULL;
    SortBed* sortBed;
    SortBed* sortExtBed;
  fp = fopen(this -> ifilePath, "r");
  int max_buffer_line = MAX_BUFFER_LINE;
  if(save_ext_len){
      max_buffer_line = MAX_BUFFER_LINE / 2;
  }
  if(unique||sort){
      sortBed = new SortBed(this -> ofilePath,unique,max_buffer_line);
  }else{
    fp_out = fopen(this -> ofilePath, "w");
  }
  if(save_ext_len){
      sortExtBed = new SortBed((std::string(this -> ofilePath)+std::string(".ext")).c_str(),unique,max_buffer_line);
  }
  char bedlineBuffer[SAM_MAX_LINE_LENGTH];
  char * line = (char*)calloc(SAM_MAX_LINE_LENGTH, 1);
  char * line1 = (char*)calloc(SAM_MAX_LINE_LENGTH, 1);
  char * tok = (char *)"start";
  char * tok1 = (char *)"start1";

  char strand;
  //char strand1;
  char * chr;
  //char * chr1;

  int  chr_start, chr_end, flag, mqs;
  int  chr_start1/*, chr_end1, flag1, mqs1*/;
  char *CIGAR = NULL;

  bool first = true;
  int freg_len;

  std::string pattern;
  if(char_filter_size>=1){

      pattern=chrList[0];
      if(char_filter_size>1){
          std::stringstream ss;
          ss << pattern ;
          for(int i = 1; i < char_filter_size; i++ ){
              ss << "|" << chrList[i] ;
          }
          ss >> pattern;
      }

  }else{
      pattern="";
  }

#ifdef PLF_SYS_WIN
  std::regex re(pattern);
#elif PLF_SYS_LINUX
  const char * patt = pattern.c_str();
  regex_t reg;
  const size_t nmatch = 1;
  regmatch_t pm[1];
  regcomp(&reg,patt,REG_EXTENDED|REG_NOSUB);
#endif
  //std::string pattern1("(chrM|chrUn.*|.*random.*)");
  //std::regex re(pattern);
  //std::string teststr="chr9_gl000199_random";
  //if(std::regex_match(teststr, re)){
//      std::cout <<"match"<<std::endl;
 // } else{
  //    std::cout << "not match"<<std::endl;
  //}
  //std::cout.flush();
      while (true)
      {
    		if(first){
    			if(!fgets(line, SAM_MAX_LINE_LENGTH, fp)){
    				break;
    			}
    			if (line[0] == '@'){
    				continue;
    			}
    			tok = strtok(line, "\t");
    			flag = atoi(strtok(NULL, "\t"));

    			chr = strtok(NULL, "\t");
#ifdef PLF_SYS_WIN
    			if(std::regex_match(std::string(chr), re)){
    			    continue;
    			}
#elif PLF_SYS_LINUX
    			if(char_filter_size>=1&&regexec(&reg,chr,nmatch,pm,REG_NOTBOL)!=REG_NOMATCH){
    			    continue;
    			}
#endif

    			//filted = false;

    			//for(int i = 0; i < char_filter_size; i++ ){
    			//    if(!strcmp(chr,chrList[i])){
    			//        filted = true;
    			//        break;
    			//    }
    			//}
    			//if(filted){
    			//    continue;
    			//}
    			if (chr[0] != '*')
    			{
    			  chr_start = atoi(strtok(NULL, "\t")) - 1;
    			  //chr_end = chr_start + readlen;
    			  mqs = atoi(strtok(NULL, "\t"));

    			  if ((flag & 0x10) == 0)
    			  {
    				strand = '+';
    				chr_start += pos_offset;
    			  }
    			  else
    			  {
    				strand = '-';
    				chr_start += neg_offset;
    			  }
    			  CIGAR = strtok(NULL, "\t");
    			  chr_end = chr_start + getReadsLen(CIGAR);

    			}
    			first = false;
    		}else{
    			if(!fgets(line1, SAM_MAX_LINE_LENGTH, fp)){
    				break;
    			}
    			tok1 = strtok(line1, "\t");

    			if(strcmp(tok1,tok)){
    				first = true;
    				continue;
    			}
    			flag = atoi(strtok(NULL, "\t"));

    			chr = strtok(NULL, "\t");
    			if (chr[0] != '*')
    			{
    			  chr_start1 = atoi(strtok(NULL, "\t")) - 1;

    			  mqs = atoi(strtok(NULL, "\t"));
    			  CIGAR = strtok(NULL, "\t");
    			  if ((flag & 0x10) == 0)
    			  {
    				chr_start1 += pos_offset;
    				chr_start = chr_start1;
    			  }
    			  else
    			  {
    				chr_start1 += neg_offset;
    				chr_end = chr_start1 + getReadsLen(CIGAR);
    			  }
    			  freg_len=chr_end-chr_start;

    			  if(freg_len>=min_freg_len&&freg_len<=max_freg_len){
        			  if(fp_out){
        			      fprintf(fp_out, "%s\t%d\t%d\t%s\t%d\t%c\n", chr, chr_start, chr_end, tok, mqs, strand);
        			  }else{
        			      sprintf(bedlineBuffer,"%s\t%d\t%c",tok, mqs, strand);
        			      sortBed->insertBedLine(new BedLine(chr,chr_start,chr_end,bedlineBuffer));
        			  }
    			  }else if(save_ext_len){
    			      sprintf(bedlineBuffer,"%s\t%d\t%c",tok, mqs, strand);
    			      sortExtBed->insertBedLine(new BedLine(chr,chr_start,chr_end,bedlineBuffer));
    			  }
    			  first = true;
    			}

    		}


  }


  if (line){
      free(line);
  }
 if (line1){
     free(line1);
 }


  fclose(fp);
 if(fp_out){
     fclose(fp_out);
 }else{
     sortBed->mergeBed();
     delete sortBed;
 }
 if(sortExtBed){
     sortExtBed->mergeBed();
     delete sortExtBed;
 }
 std::cout<<"finish"<<std::endl;

#ifdef PLF_SYS_LINUX
 regfree(&reg);
#endif

  return 0;
}
