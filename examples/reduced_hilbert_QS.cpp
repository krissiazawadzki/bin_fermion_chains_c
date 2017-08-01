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
#include "prim.h"
/* custom aux libraries */
#include "int_array_util.h"


/* printing binary configurations info */
#include "io_bin.h"





int main(int argc, char *argv[]){
	
	int nsites = 2;
	int Q = nsites;
	int dS = nsites%2;
	std::string folder_results = "/home/krissia/Desktop/cond_matter_module/bin_spin_chains/rotation_matrices_QS_cpp";
	
	if(argc == 2){
		nsites = atoi(argv[1]);
	}
	Q = nsites;

	if(argc == 3){
		nsites = atoi(argv[1]);
		Q = atoi(argv[2]);
	}	
	dS = nsites%2;

	if(argc == 4){
		nsites = atoi(argv[1]);
		Q = atoi(argv[2]);
		dS = atoi(argv[3]);
	}
	if(argc == 5){
		nsites = atoi(argv[1]);
		Q = atoi(argv[2]);
		dS = atoi(argv[3]);
		folder_results = argv[4];
	}
		
	std::cout << "\n\ttotal number of sites: " << nsites << std::endl;
	std::cout << "\n\ttarget (Q,S) : (" << Q << ", " << dS << ")" << std::endl;
	
	int Nmax = nsites;
	set_itermax(Nmax);
	

	// initializing 
	init_QSmatrix();
	generate_single_site_hilbert_space();
	update_from_single_site_to_dimer();
	
	
	
	int n = iter();
	
	while(n <= iter_max() ){
		cout << "\n\tBuilding Hilbert space for " << n << " sites " << endl;
		cout << " \t  \t " << endl;
		
		set_iter(n);
		primitive_basis(QS_matrix());
		
		
		cout << " * * * * * " << endl;
		
		//print_all_sectors_friendly();
		
		//print_all_sectors();
		if(n == iter_max()){
			std::cout << "target Hilbert space info (Q,S) = (" << Q << ", " << dS << ")" << std::endl;
			print_target_hilbert_space(Q,dS, false, folder_results);
		}
		
		
		update_QSmatrix();
		update_q_dS();	
		
		n++;
		}
		
	update_QSmatrix();
	
	clear_QSmatrix(QS_matrix());	
	
	std::cout << "END OF CALCULATION" << std::endl;
	
	return 0;
}
