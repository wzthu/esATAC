#ifndef _RCOUT_RCERR_H
#define _RCOUT_RCERR_H
#ifdef R_ENV_FLAG
#include <Rcpp.h>
using Rcpp::Rcout;
using Rcpp::Rcerr;
#define cout Rcout
#define cerr Rcerr
#else
using namespace std;
using std::cout;
using std::cerr;
#endif
#endif
