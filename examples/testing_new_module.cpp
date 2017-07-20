/* basic standard routines */
#include <cstdlib>
#include <iostream>

/* IO */
#include <fstream>
#include <iomanip>
#include <sstream>


/* math */
#include <cmath>

/* string */
#include <cstring>


/* custom libraries */
#include "spin_configurations.h"
#include "bin_sectors.h"
#include "qs_hilbert.h"
#include "QSmatrix.h"
#include "first_iter.h"

#include "io_bin.h"

int main(int argc, char *argv[]){
	
	int nsites = 2;
	int Q = nsites;
	int dS = nsites%2;
	
	if(argc == 2){
		nsites = atoi(argv[1]);
	}
	Q = nsites;
	dS = nsites%2;
	
	
	std::cout << "total number of sites: " << nsites << std::endl;
	std::cout << "target (Q,S) : (" << Q << ", " << dS << ")" << std::endl;
	
	int Nmax = nsites;
	set_itermax(Nmax);
	std::cout << "MAX NUM OF ITERATIONS "  << iter_max() << std::endl;	


	// initializing 
	init_QSmatrix();
	generate_single_site_hilbert_space();
	
	std::cout << "current iter " << iter() << std::endl;
	print_all_sectors();
	
	update_from_single_site_to_dimer();
		
	return 0;
}
