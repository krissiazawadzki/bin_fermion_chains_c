#ifndef __BIN_SECTORS__
#include <vector>
#include <map>

using namespace std;


struct bin_sector{
  int Q; // charge
  int dS; // 2*spin
  int np; // number of primitive states: dimension of reduced Hilbert space
  int nb; // number of total binary configurations inside sector
  int nub; // effective number of binary configurations for each sector
  int *bin_confs; //  array with the binary cconfigurations inside sector
  double *t_coefs; // array with the coefficients of the translation matrix
  int *np_binconfs; // for each state np, it gives the number of binary states that generated it
  vector<int> npv;
  bin_sector() :npv(4,0) {};
}; 
	

	
// accessing fields of the sector data structure
int Q(bin_sector sec);

int dS(bin_sector sec);

int np(bin_sector sec);

int nb(bin_sector sec);

int nub(bin_sector sec);

vector<int> npv(bin_sector sec);

int sum_npv(vector<int> npv);

int npN(bin_sector sec);

int npE(bin_sector sec);

int npS(bin_sector sec);

int npW(bin_sector sec);

vector<int> initial_index_in_bin_confs(bin_sector sec);


// obtaining relations between bin_sectors

bin_sector *parent(bin_sector sect, int p);

int gender(bin_sector sect, int p);

int index_invar(bin_sector sec_bra, int ket_index, int bra_index);

int bra_index(bin_sector sec, int index);

int ket_index(bin_sector sec, int index);

// printing a bin_sector
void print_sector(bin_sector sec, int ds, int curr=0, int print_translation_mat = 0);

void print_sector_translation(bin_sector sec);

// deleting a sector
void clear_sector(bin_sector &sec, bool delete_vectors = false);



#define __BIN_SECTORS__
#endif
