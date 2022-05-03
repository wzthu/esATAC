#include<string>
#include<fstream>
#include<iostream>
class SC_Renamer
{
private:
    char * ifilePath;
    char * ofilePath;
    const int BUFFER_LENGTH=10000;
public:
    SC_Renamer(char * ifilePath,char * ofilePath);
    void renameFasta();
    void renameFastq();
    void renameInterleaveFastq();
};
