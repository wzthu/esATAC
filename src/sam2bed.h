#include<string>
#include<fstream>
class SamToBed
{
private:
  char * ifilePath;
  char * ofilePath;
  int readlen;
public:
  SamToBed(char * ifilePath, char * ofilePath, int readlen);
  int sam2bed();
};
