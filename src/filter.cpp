////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////This file using the following license////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// MIT License
//
// Copyright (c) 2021 Tim Stuart
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy
//     of this software and associated documentation files (the "Software"), to deal
//     in the Software without restriction, including without limitation the rights
//     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//     copies of the Software, and to permit persons to whom the Software is
//     furnished to do so, subject to the following conditions:
//
//         The above copyright notice and this permission notice shall be included in all
//         copies or substantial portions of the Software.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//         IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//         FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//         AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//         LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//         OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//         SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////////////////////

#include <Rcpp.h>
#include <zlib.h>
#include <iostream>
#include <fstream>


// [[Rcpp::export]]
int filterCells(
    std::string fragments,
    std::string outfile,
    std::vector<std::string> keep_cells,
    int buffer_length,
    bool verbose = true
) {
  // opening gzipped compressed stream
  gzFile ifileHandler = gzopen(fragments.c_str(), "rb");
  std::ofstream ofileHandler;
  ofileHandler.open(outfile.c_str());

  // return 1 if it can't find the file
  if (ifileHandler == NULL) {
    Rcpp::Rcerr << "can't open file" << std::flush;
    return 1;
  }

  // C based buffered string parsing
  char* cb_char;
  size_t line_counter {1};
  char *buffer = new char[buffer_length];

  // Hash Map storing the barcodes to keep
  std::unordered_set<std::string> index_hash(keep_cells.begin(), keep_cells.end());

  size_t num_whitelist_cells {0};
  {

    if (verbose) {
      num_whitelist_cells = index_hash.size();
      Rcpp::Rcerr << "Keeping " << num_whitelist_cells
                  << " cell barcodes"
                  << std::endl << std::flush;
    }
  }

  // char * to string extraction
  std::string cb_seq, line_seq;
  cb_seq.reserve(32);
  line_seq.reserve(buffer_length);

  bool eof_check;
  while ((eof_check = gzgets(ifileHandler, buffer, buffer_length)) !=0) {
    line_seq.clear();
    line_seq.append(buffer);

    if (line_seq.at(0) != '#') {
      break;
    }
  }

  if (!eof_check) {
    Rcpp::Rcerr << "Error: fragment file contains header only\n" << std::flush;
    gzclose(ifileHandler);
    return 1;
  }

  // looping over the fragments file
  do {
    line_seq.append(buffer);

    cb_char = strtok ( buffer, "\t" );
    for (auto i=1; i<=3; i++) {
      cb_char = strtok (NULL, "\t");

      if(i == 3) {
        cb_seq.clear();
        cb_seq.append(cb_char);
        if (index_hash.count(cb_seq) > 0) {
          ofileHandler << line_seq.c_str();
        }
      }
    }

    line_counter += 1;
    bool is_ten_mil = line_counter % 10000000 == 0;
    if (verbose) {
      if (is_ten_mil) {
        Rcpp::Rcerr << "\r                                                  ";
      }

      if (line_counter % 1000000 == 0) {
        Rcpp::Rcerr << "\rDone Processing " << line_counter / 1000000
                    << " million lines";
      }
    }

    if (is_ten_mil) {
      Rcpp::checkUserInterrupt();
    }

    line_seq.clear();
  } while(gzgets(ifileHandler, buffer, buffer_length) !=0 );

  //Cleanup
  gzclose(ifileHandler);
  ofileHandler.close();

  return 0;
}
