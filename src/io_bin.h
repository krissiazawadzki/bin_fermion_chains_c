#ifndef __IO_BIN__

#include <string>

#include "bin_sectors.h"
#include "QSmatrix.h"

void print_all_sectors();

void print_sector_friendly(bin_sector sec);

void print_all_sectors_friendly();


void print_target_hilbert_space(int Q, int dS, bool print_bin_confs=false, std::string folder_results="", int wfmt = 10, int precisionfmt=15);


#define __IO_BIN__
#endif
