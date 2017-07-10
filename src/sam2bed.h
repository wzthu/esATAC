#include<string>
#include<fstream>
class SamToBed
{
private:
  char * ifilePath;
  char * ofilePath;

public:
  SamToBed(char * ifilePath, char * ofilePath);
  int sam2bed();

};
