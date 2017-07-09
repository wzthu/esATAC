#ifndef MAIN_DEMULTIPLEX_H
#define MAIN_DEMULTIPLEX_H

#include<fstream>
#include "demultiplex.h"
#include "userconfig.h"
namespace ar
{
 int demultiplex_sequences(const userconfig& config);
void write_demultiplex_statistics(std::ofstream& output,
                                  const userconfig& config,
                                  const demultiplex_reads* step);
}

#endif


