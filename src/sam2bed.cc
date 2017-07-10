/***************************************************************

 The Subread and Rsubread software packages are free
software packages:

you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

Subread is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details.

Authors: Drs Yang Liao and Wei Shi

***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam2bed.h"

SamToBed::SamToBed(char * ifilePath, char * ofilePath){
  this -> ifilePath = ifilePath;
  this -> ofilePath = ofilePath;
}

int SamToBed::sam2bed() {

  FILE *fp, *fp_out;

  fp = fopen(this -> ifilePath, "r");
  fp_out = fopen(this -> ofilePath, "w");

  char * line = NULL;
  char * tok;
  char strand;
  char * chr;
  char * reads;
  int  chr_start, chr_end, flag, mqs;

  int MAX_LINE_LENGTH = 100000;



  line = (char*)calloc(MAX_LINE_LENGTH, 1);

  while (fgets(line, MAX_LINE_LENGTH, fp))
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
	  for(int i=0;i<5;i++){
		  reads = strtok(NULL, "\t");
	  }
	  chr_end = chr_start + strlen(reads);
      fprintf(fp_out, "%s\t%d\t%d\t%s\t%d\t%c\n", chr, chr_start, chr_end, tok, mqs, strand);
    }
  }

  if (line)
    free(line);

  fclose(fp);
  fclose(fp_out);

  return 0;
}
