#ifdef R_EVN_FLAG
#include <Rcpp.h>
using Rcpp::Rcout;
using Rcpp::Rcerr;
#define cout Rcout
#define cerr Rcerr
#else
#include<iostream>
using namespace std;
using std::cout;
using std::cerr;
#endif
