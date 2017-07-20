#ifndef __PRIM__
#include <vector>
#include "bin_sectors.h"

using namespace std;

void primitive_basis(vector<vector<bin_sector> > &qsMatrix);

void O_SOUTH(bin_sector &sectorchild, bin_sector sectorparent);
void O_NORTH(bin_sector &sectorchild, bin_sector sectorparent);
void O_WEST(bin_sector &sectorchild, bin_sector sectorparent);
void O_EAST(bin_sector &sectorchild, bin_sector sectorparent);

int spin_flip(int binary_index, int nsites, int *flipped_states);

vector<int> spin_flip_long(int binary_index, int nsites, int binary_index_pt2=0, int nsites_pt2=0);

#define __PRIM__
#endif
