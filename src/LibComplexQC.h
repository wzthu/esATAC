#include<string>
using std::string;
class LibComplexQC
{
private:
    string unsorted_bed;
    string sorted_bed;
    int max_reads;
    float NRF;
    float PBC1;
    float PBC2;
    int exact1;
    int exact2;
    int exactT;
    int reads;
public:
    LibComplexQC(string unsorted_bed,int max_reads);
    LibComplexQC(string sorted_bed);
    ~LibComplexQC(void);
    void calValSorted();
    void calValUnSorted();
    float getNRF();
    float getPBC1();
    float getPBC2();
    int getOne();
    int getTwo();
    int getTotal();
    int getReads();
};
