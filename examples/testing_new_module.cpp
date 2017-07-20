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

int main(int argc, char *argv[]){
	
	int nsites = 2;
	int Q = nsites;
	int dS = nsites%2;
	
	if(argc == 2){
		nsites = atoi(argv[1]);
	}
	
	std::cout << "total number of sites: " << nsites << std::endl;
	std::cout << "target (Q,S) : (" << Q << ", " << dS << ")" << std::endl;
	
	
	return 0;
}
