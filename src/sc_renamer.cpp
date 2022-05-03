#include "sc_renamer.h"
SC_Renamer::SC_Renamer(char * ifilePath,char * ofilePath){
        this -> ifilePath = ifilePath;
        this -> ofilePath = ofilePath;
    }
void SC_Renamer::renameFasta(){
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
void SC_Renamer::renameFastq(){
    std::ifstream fin(this -> ifilePath);
    std::ofstream fout(this -> ofilePath);
    char line[BUFFER_LENGTH];
    char new_line[BUFFER_LENGTH];
    long i=0L;
    int t = 0;
    int t1 = 0;
    int t2 = 0;
    int k =0;
    while(fin.getline(line,BUFFER_LENGTH)){
        if(i%2L==0L){
            if(i%4L==0L){
                t = 0;
                while(line[t] != '\0' && line[t] != ' ') t++;
                t2 = t;
                while(line[t2] != '\0') t2++;
                t1 = t2;
                while(t1 > t && line[t1] != ':' && line[t1] != ' ') t1--;
                if (line[t] == '\0'){
                    for(int j=0;j<=t;j++){
                        new_line[j] = line[j];
                    }
                }else if(line[t] == ' '){
                    if(line[t1] == ':' || line[t1] == ' '){
                        k = 0;
                        for(int j = t1+1; j < t2; j++){
                            new_line[k] =  line[j];
                            k++;
                        }
                        new_line[k]=':';
                        k+=1;
                        for(int j = 1; j < t; j++){
                            new_line[k] = line[j];
                            k++;
                        }
                        new_line[k]='\0';
                    }
                }
                fout <<"@"<< new_line << std::endl;
            }else{
                fout <<"+"<< std::endl;
            }
        }else{
            fout << line << std::endl;
        }
        i++;
    }
}
void SC_Renamer::renameInterleaveFastq(){
    std::ifstream fin(this -> ifilePath);
    std::ofstream fout(this -> ofilePath);
    char line[BUFFER_LENGTH];
    char new_line[BUFFER_LENGTH];
    long i=0L;
    int t = 0;
    int t1 = 0;
    int t2 = 0;
    int k =0;
    while(fin.getline(line,BUFFER_LENGTH)){
        if(i%2L==0L){
            if(i%4L==0L){
                t = 0;
                while(line[t] != '\0' && line[t] != ' ') t++;
                t2 = t;
                while(line[t2] != '\0') t2++;
                t1 = t2;
                while(t1 > t && line[t1] != ':' && line[t1] != ' ') t1--;
                if (line[t] == '\0'){
                    for(int j=0;j<=t;j++){
                        new_line[j] = line[j];
                    }
                }else if(line[t] == ' '){
                    if(line[t1] == ':' || line[t1] == ' '){
                        k = 0;
                        for(int j = t1+1; j < t2; j++){
                            new_line[k] =  line[j];
                            k++;
                        }
                        new_line[k]=':';
                        k+=1;
                        for(int j = 1; j < t; j++){
                            new_line[k] = line[j];
                            k++;
                        }
                        new_line[k]='\0';
                    }
                }
                fout <<"@"<< new_line << std::endl;
            }else{
                fout <<"+"<< std::endl;
            }
        }else{
            fout << line << std::endl;
        }
        i++;
    }
}

