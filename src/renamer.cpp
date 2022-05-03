#include "renamer.h"
Renamer::Renamer(char * ifilePath,char * ofilePath){
        this -> ifilePath = ifilePath;
        this -> ofilePath = ofilePath;
    }
void Renamer::renameFasta(){
    std::ifstream fin(this -> ifilePath);
    std::ofstream fout(this -> ofilePath);
    char line[BUFFER_LENGTH];
    long i=0L;
    while(fin.getline(line,BUFFER_LENGTH)){
        if(i%2L==0L){
            fout << (i/2L+1L) << std::endl;
        }else{
            fout << line << std::endl;
        }
        i++;
    }
}
void Renamer::renameFastq(){
    std::ifstream fin(this -> ifilePath);
    std::ofstream fout(this -> ofilePath);
    char line[BUFFER_LENGTH];
    long i=0L;
    while(fin.getline(line,BUFFER_LENGTH)){
        if(i%2L==0L){
            if(i%4L==0L){
                fout <<"@"<< (i/4L+1L) << std::endl;
            }else{
                fout <<"+"<< std::endl;
            }
        }else{
            fout << line << std::endl;
        }
        i++;
    }
}
void Renamer::renameInterleaveFastq(){
    std::ifstream fin(this -> ifilePath);
    std::ofstream fout(this -> ofilePath);
    char line[BUFFER_LENGTH];
    long i=0L;
    while(fin.getline(line,BUFFER_LENGTH)){
        if(i%2L==0L){
            if(i%4L==0L){
                fout <<"@"<< (i/8L+1L) << std::endl;
            }else{
                fout <<"+"<< std::endl;
            }
        }else{
            fout << line << std::endl;
        }
        i++;
    }
}

