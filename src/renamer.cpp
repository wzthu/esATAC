#include "renamer.h"
Renamer::Renamer(char * ifilePath,char * ofilePath){
        this -> ifilePath = ifilePath;
        this -> ofilePath = ofilePath;
    }
int Renamer::renameFasta(){
       std::ifstream fin(this -> ifilePath);
       std::ofstream fout(this -> ofilePath);
       char line[BUFFER_LENGTH];
       int i=0;
       while(fin.getline(line,BUFFER_LENGTH)){
           if(i%2==0){
               fout << (i/2+1) << std::endl;
           }else{
               fout << line << std::endl;
           }
           i++;
       }
    }
int Renamer::renameFastq(){
        std::ifstream fin(this -> ifilePath);
        std::ofstream fout(this -> ofilePath);
        char line[BUFFER_LENGTH];
        int i=0;
        int j=0;
        while(fin.getline(line,BUFFER_LENGTH)){
            if(i%2==0){
		if(i%4==0){
                	fout <<"@"<< (i/4+1) << std::endl;
		}else{
			fout <<"+"<< (i/4+1) << std::endl;
		}
            }else{
                fout << line << std::endl;
            }
            i++;
        }
    }

