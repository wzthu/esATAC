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

// [[Rcpp::export]]
bool validateCells(
    std::string fragments,
    std::vector<std::string> cells,
    std::size_t find_n,
    std::size_t max_lines = 0,
    bool verbose = true
) {
  // opening gzipped compressed stream
  gzFile fileHandler = gzopen(fragments.c_str(), "rb");

  // determine if we read the whole file or the first n lines
  bool read_part {false};
  if (max_lines > 0) {
    read_part = true;
  }

  // return false it can't find the file
  if (fileHandler == NULL) {
    Rcpp::Rcerr << "can't open file" << std::flush;
    gzclose(fileHandler);
    return (false);
  }

  // C based buffered string parsing
  char* cb_char;
  size_t line_counter {1};
  size_t total_seen {0};
  uint32_t buffer_length = 256;
  char *buffer = new char[buffer_length];

  // Hash Map storing the barcodes to look for
  std::unordered_set<std::string> index_hash(cells.begin(), cells.end());

  // work out how many cells we need to see before returning true
  size_t ncell = index_hash.size();
  {
    if (verbose) {
      Rcpp::Rcerr << "Checking for " << ncell
                  << " cell barcodes"
                  << std::endl << std::flush;
    }
  }

  // char * to string extraction
  std::string cb_seq, line_seq;
  cb_seq.reserve(32);
  line_seq.reserve(buffer_length);

  // skip header if present
  bool eof_check;
  while ((eof_check = gzgets(fileHandler, buffer, buffer_length)) !=0) {
    line_seq.clear();
    line_seq.append(buffer);

    if (line_seq.at(0) != '#') {
      break;
    }
  }

  if (!eof_check) {
    Rcpp::Rcerr << "Error: fragment file contains header only\n" << std::flush;
    gzclose(fileHandler);
    return (false);
  }

  // looping over the fragments file
  do {
    cb_char = strtok ( buffer, "\t" );

    for (auto i=1; i<=3; i++) {
      cb_char = strtok (NULL, "\t");
      if(i == 3) {
        cb_seq.clear();
        cb_seq.append(cb_char);
        auto it = index_hash.find(cb_seq);
        if (it != index_hash.end()) {
          // cell exists in the set, remove from hash map
          index_hash.erase(it);
          total_seen++;
        }
      }
    }

    if (total_seen >= find_n) {
      gzclose(fileHandler);
      return(true);
    }

    line_counter += 1;
    if (read_part) {
      if (line_counter > max_lines) {
        gzclose(fileHandler);
        return(false);
      }
    }

    if (line_counter % 2000000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  } while(gzgets(fileHandler, buffer, buffer_length) !=0 );

  //Cleanup
  gzclose(fileHandler);

  return (false);
}
