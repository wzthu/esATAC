#include<string>
#include<fstream>
class Renamer
{
private:
    char * ifilePath;
    char * ofilePath;
    const int BUFFER_LENGTH=10000;
public:
    Renamer(char * ifilePath,char * ofilePath);
    void renameFasta();
    void renameFastq();
    void renameInterleaveFastq();
};
