#include "prim.h"

/* dependences */
#include "spin_configurations.h"
#include "qs_hilbert.h"
#include "first_iter.h"

/* auxiliary packages */
#include "matrix_ops/matrix_ops.h"

#include <algorithm>
#include <iostream>
#include <cstring>
#include <cmath>


void primitive_basis(vector<vector<bin_sector> > &qsMatrix)
/*	void primitive_basis(vector<vector<sector> > &qsMatrix)
 * 
 * 	Function constructs primitive basis for current iteration. qsMatrix 
 * 	contains all (q, ds) sectors generated in previous iteration. For each
 * 	eigenstate in such sector, generate four basis states, of genders
 * 	N, E, S, W.
 * 
 * 	Inputs:
 * 		vector<vector<sector> > &qsMatrix - QS_matrix containing sectors 
 * 			from previous iteration
 * 
 * 	Output:
 * 		none
 * */
{
  int qmax_previous = qmax_prev();
  for(int q = qmax_previous; q >= -qmax_previous; q--){		
    int qind = q_index(q);		// makes indices non-negative
    if(qind >= 0 && qind < 2*QMAX+1){
      int dsmax_p = dsmax_prev(q);
      for(int ds = dsmax_p%2; ds <= dsmax_p; ds +=2){
		int np = qsMatrix[qind][ds].np;
		//int nb = qsMatrix[qind][ds].nb;
		if(np == 0){ 
			continue;
			} // if there is no state in sector Q,S we don't have to compute sector
		// NORTH
		if(qind-1 >= 0 ){ // assures the access of valid index in Qsmatrix
		  (qsMatrix[qind-1][ds]).Q = q+1;
		  (qsMatrix[qind-1][ds]).dS = ds;
		  O_NORTH(qsMatrix[qind-1][ds], qsMatrix[qind][ds]);
		  (qsMatrix[qind-1][ds]).npv[NORTH] = np;
		  //(qsMatrix[qind-1][ds]).nr += nub;
		  (qsMatrix[qind-1][ds]).np = sum_npv((qsMatrix[qind-1][ds]).npv);
			}
		// EAST	
		if(ds+1 < DSMAX){ // assures the access of valid index in Qsmatrix
		  (qsMatrix[qind][ds+1]).Q = q;
		  (qsMatrix[qind][ds+1]).dS = ds+1;
		  O_EAST(qsMatrix[qind][ds+1], qsMatrix[qind][ds]);
		  (qsMatrix[qind][ds+1]).npv[EAST] = np;
		  //(qsMatrix[qind][ds+1]).nr += nr;
		  (qsMatrix[qind][ds+1]).np = sum_npv((qsMatrix[qind][ds+1]).npv);
			}
		// SOUTH
		if(qind+1 < 2*QMAX+1){ // assures the access of valid index in Qsmatrix
		  (qsMatrix[qind+1][ds]).Q = q-1;
		  (qsMatrix[qind+1][ds]).dS = ds;
		  O_SOUTH(qsMatrix[qind+1][ds], qsMatrix[qind][ds]);
		  (qsMatrix[qind+1][ds]).npv[SOUTH] = np;
		  //(qsMatrix[qind+1][ds]).nr += np;
		  (qsMatrix[qind+1][ds]).np = sum_npv((qsMatrix[qind+1][ds]).npv);
			}
		// WEST
		if(ds-1 >=0){ // assures the access of valid index in Qsmatrix
		  (qsMatrix[qind][ds-1]).Q = q;
		  (qsMatrix[qind][ds-1]).dS = ds-1;
		  O_WEST(qsMatrix[qind][ds-1], qsMatrix[qind][ds]);
		  (qsMatrix[qind][ds-1]).npv[WEST] = np;
		  //(qsMatrix[qind][ds-1]).nr += nr;
		  (qsMatrix[qind][ds-1]).np = sum_npv((qsMatrix[qind][ds-1]).npv);
			}
		  }
		}
	  }
	  
}



int spin_flip(int binary_index, int nsites, int *flipped_states)
/*	*spin_flip(int *binary_state)
 * 
 * 	Function constructs all the configurations required to build states 
 * 	with definite spins when Sz = S - 1 is required.
 * 
 * 	Inputs:
 * 		int binary_index - index of binary notation;
 * 		int nsites - number of sites in the configuration
 * 
 * 	Output:
 * 		vector<int> - indexes of configurations with flipped_spins up
 * */	
{	
	int q_binstate;
	int *fns_binstate = new int[2*nsites];
	int *binstate = index_to_binary(binary_index, nsites, fns_binstate, q_binstate, 0);
	int *copy_flipped = new int[2*nsites];
	
	
	int numflips = 0;
	for(int i = 0; i < nsites; i++){
		int bit_up = binstate[2*i];
		int bit_down = binstate[2*i+1];
		if(bit_up == 1 && bit_down == 0){
			memcpy(copy_flipped, binstate, sizeof(int)*2*nsites);
			copy_flipped[2*i] = 0;
			copy_flipped[2*i+1] = 1;
			int index_flipped = binary_to_index(copy_flipped, nsites);
			flipped_states[numflips] = index_flipped;
			numflips++;
			}
		}
	
	delete[] fns_binstate;
	delete[] binstate;
	delete[] copy_flipped;
	
	return numflips;
}


int *my_push_back_int(int *oldarray, int oldsize, int element){
	int *newarray = new int[oldsize+1];
	memcpy(newarray, oldarray, sizeof(int)*oldsize);
	newarray[oldsize] = element;
	delete[] oldarray;
	return newarray;
}

double *my_push_back_double(double *oldarray, int oldsize, double element){
	double *newarray = new double[oldsize+1];
	memcpy(newarray, oldarray, sizeof(double)*oldsize);
	newarray[oldsize] = element;
	delete[] oldarray;
	return newarray;
}


void O_NORTH(bin_sector &sectorchild, bin_sector sectorparent){
	
	int nsites_curr = iter();
	
	int np_parent = sectorparent.np; 
	int nb_parent = sectorparent.nb;
	
	int np_child = sectorchild.np;
	int nb_child = sectorchild.nb;
	
	int *bin_confs = new int[nb_child + nb_parent];
	double *t_coefs = new double[nb_child + nb_parent];
	int *np_binconfs = new int[np_child + np_parent];

  for(int i = 0; i < nb_child; i++){
    bin_confs[i] = sectorchild.bin_confs[i];
    t_coefs[i] = sectorchild.t_coefs[i];
		}
  for(int p = 0; p < np_child; p++){
    np_binconfs[p] = sectorchild.np_binconfs[p];
		}
  
  /* I decided to avoid memcpy
   
   if(nb_child != 0){
   memcpy(bin_confs, (sectorchild.bin_confs) , sizeof(int)*nb_child);
   memcpy(t_coefs, (sectorchild.t_coefs) , sizeof(double)*nb_child);
   }
   if(np_child != 0){
   memcpy(np_binconfs, (sectorchild.np_binconfs) , sizeof(int)*np_child);
   }
   
   memcpy(bin_confs+nb_child, (sectorparent.bin_confs) , sizeof(int)*nb_parent);
   memcpy(t_coefs+nb_child, (sectorparent.t_coefs) , sizeof(double)*nb_parent);
   memcpy
	
   */
	for(int b = 0; b < nb_parent; b++){
		double tpb = sectorparent.t_coefs[b];
		int binary_index_child = pow(2, 2*nsites_curr-1) + pow(2, 2*nsites_curr -2) + sectorparent.bin_confs[b];
		t_coefs[nb_child + b] = tpb;
		bin_confs[nb_child + b] = binary_index_child;
		}
	for(int p = 0; p < np_parent; p++){
		np_binconfs[np_child + p] = sectorparent.np_binconfs[p];
		}
	
	delete[] (sectorchild.bin_confs);
	delete[] (sectorchild.t_coefs);
	delete[] (sectorchild.np_binconfs);
	
	(sectorchild.bin_confs) = bin_confs;
	(sectorchild.t_coefs) = t_coefs;
	(sectorchild.np_binconfs) = np_binconfs;

	sectorchild.nb += nb_parent;
	
}



void O_SOUTH(bin_sector &sectorchild, bin_sector sectorparent){

	
	int np_parent = sectorparent.np; 
	int nb_parent = sectorparent.nb;
	
	int np_child = sectorchild.np;
	int nb_child = sectorchild.nb;
	
	int *bin_confs = new int[nb_child + nb_parent];
	double *t_coefs = new double[nb_child + nb_parent];
	int *np_binconfs = new int[np_child + np_parent];
	
	for(int i = 0; i < nb_child; i++){
		bin_confs[i] = sectorchild.bin_confs[i];
		t_coefs[i] = sectorchild.t_coefs[i];
		}
	for(int p = 0; p < np_child; p++){
		np_binconfs[p] = sectorchild.np_binconfs[p];
		}
	
	/* I decided to avoid memcpy
	
	if(nb_child != 0){
		memcpy(bin_confs, (sectorchild.bin_confs) , sizeof(int)*nb_child);
		memcpy(t_coefs, (sectorchild.t_coefs) , sizeof(double)*nb_child);	
	}
	if(np_child != 0){
		memcpy(np_binconfs, (sectorchild.np_binconfs) , sizeof(int)*np_child);
	}
	
	memcpy(bin_confs+nb_child, (sectorparent.bin_confs) , sizeof(int)*nb_parent);
	memcpy(t_coefs+nb_child, (sectorparent.t_coefs) , sizeof(double)*nb_parent);
	memcpy(np_binconfs+np_child, (sectorparent.np_binconfs) , sizeof(int)*np_parent);
	* */
	
	for(int b = 0; b < nb_parent; b++){
		double tpb = sectorparent.t_coefs[b];
		int binary_index_child = sectorparent.bin_confs[b];
		t_coefs[nb_child + b] = tpb;
		bin_confs[nb_child + b] = binary_index_child;
		}
	for(int p = 0; p < np_parent; p++){
		np_binconfs[np_child + p] = sectorparent.np_binconfs[p];
		}
	
	delete[] (sectorchild.bin_confs);
	delete[] (sectorchild.t_coefs);
	delete[] (sectorchild.np_binconfs);
	
	(sectorchild.bin_confs) = bin_confs;
	(sectorchild.t_coefs) = t_coefs;
	(sectorchild.np_binconfs) = np_binconfs;

	sectorchild.nb += nb_parent;
	
}


vector<int> spin_flip_vector(int binary_index, int nsites)
/*	*spin_flip(int *binary_state)
 * 
 * 	Function constructs all the configurations required to build states 
 * 	with definite spins when Sz = S - 1 is required.
 * 
 * 	Inputs:
 * 		int binary_index - index of binary notation;
 * 		int nsites - number of sites in the configuration
 * 
 * 	Output:
 * 		vector<int> - indexes of configurations with flipped_spins up
 * */	
{	
	int q_binstate;
	int *fns_binstate = new int[2*nsites];
	int *binstate = index_to_binary(binary_index, nsites, fns_binstate, q_binstate, 0);
	int *copy_flipped = new int[2*nsites];
	
	vector<int> flipped_states;
	for(int i = 0; i < nsites; i++){
		int bit_up = binstate[2*i];
		int bit_down = binstate[2*i+1];
		if(bit_up == 1 && bit_down == 0){
			memcpy(copy_flipped, binstate, sizeof(int)*2*nsites);
			copy_flipped[2*i] = 0;
			copy_flipped[2*i+1] = 1;
			int index_flipped = binary_to_index(copy_flipped, nsites);
			flipped_states.push_back(index_flipped);
			}
		}
	
	delete[] fns_binstate;
	delete[] binstate;
	delete[] copy_flipped;
	
	return flipped_states;
}



void O_WEST(bin_sector &sectorchild, bin_sector sectorparent){
	
	int nsites_curr = iter();
	int nsites_prev = nsites_curr-1;
	int dS = sectorparent.dS;
	int Q = sectorparent.Q;
	

	// first, we simply copy all parent states with a spin down
	int np_parent = sectorparent.np; 
	int nb_parent = sectorparent.nb;
	
	int np_child = sectorchild.np;
	int nb_child = sectorchild.nb;
	
	int maxnumflipsall = nb_parent * ( Q + nsites_curr + dS + 2 ) / 2;
	
	int *bin_confs = new int[nb_child + nb_parent + maxnumflipsall];
	double *t_coefs = new double[nb_child + nb_parent + maxnumflipsall];
	int *np_binconfs = new int[np_child + np_parent + maxnumflipsall];
	
	
	for(int i = 0; i < nb_child; i++){
		bin_confs[i] = sectorchild.bin_confs[i];
		t_coefs[i] = sectorchild.t_coefs[i];
		}
	for(int p = 0; p < np_child; p++){
		np_binconfs[p] = sectorchild.np_binconfs[p];
		}
	
	/* I decided to avoid memcpy
	
	if(nb_child != 0){
		memcpy(bin_confs, (sectorchild.bin_confs) , sizeof(int)*nb_child);
		memcpy(t_coefs, (sectorchild.t_coefs) , sizeof(double)*nb_child);	
	}
	if(np_child != 0){
		memcpy(np_binconfs, (sectorchild.np_binconfs) , sizeof(int)*np_child);
	}
	* */
	
	
	
	
	
	vector<int> beg_states = initial_index_in_bin_confs(sectorparent);
	int beg_states_sum = 0;
	int nbtotal = 0;
	for(int p = 0; p < np_parent; p++){
		
		int nb_confs_p = sectorparent.np_binconfs[p];
		
		int maxnumflips = nb_confs_p * ( Q + nsites_curr + dS + 2 ) / 2;
		
		int *bin_confs_aux = new int[nb_confs_p + maxnumflips];
		double *t_coefs_aux = new double[nb_confs_p + maxnumflips];
		
		
		// in first step we just copy the binary configurations of the parent sector with a spin up
		for(int b = 0; b < nb_confs_p; b++){
			int binary_index_child_down = pow(2, 2*nsites_curr -2) + sectorparent.bin_confs[beg_states[p]+b];
			double tpb_down =  sqrt(dS / (dS+1.0)) * sectorparent.t_coefs[beg_states[p]+b];
			t_coefs_aux[b] = tpb_down; 
			bin_confs_aux[b] = binary_index_child_down;
			}
			
			
		// in second step, we construct the flipped states for each configuration
		map <int, double> flipped_map;
		
		//int *flipped_states = new int[maxnumflips];
		
		
		for(int b = 0; b < nb_confs_p; b++){
			int binary_index = sectorparent.bin_confs[beg_states[p]+b];
			double tpb_up = -  sqrt(1.0 / (dS+1.0)) * sectorparent.t_coefs[beg_states[p]+b];
			
			
			//int nflips = spin_flip(binary_index, nsites_prev, flipped_states);
			vector<int> flipped_states = spin_flip_vector(binary_index, nsites_prev);
			int nflips = flipped_states.size();
			
			//int binary_index_child_down = pow(2, 2*nsites_curr -1) + flipped_states[0];
			//flipped_map[binary_index_child_down] = tpb_up;
			
			for(int f = 0; f < nflips; f++){
				// if key doesn't exist on map, add new binconf with coefficient tpb_up
				int binary_index_child_down = pow(2, 2*nsites_curr -1) + flipped_states[f];
				if(flipped_map.find(binary_index_child_down) == flipped_map.end() ){
					flipped_map[binary_index_child_down] = tpb_up;
					}
				// if key exists 
				else{
					double coef_old = flipped_map[binary_index_child_down];
					flipped_map[binary_index_child_down] = coef_old + tpb_up;
					}
				}
			flipped_states.clear();
			}
		// now, we must check the number of coefficients that are not zero
		//cout << "before map" << endl;
		//print_array_int(bin_confs_aux, nb_confs_p);
		//print_array(t_coefs_aux, nb_confs_p);
		//cout << " " << endl;
		int bf = 0;
		for(map<int,double>::iterator it = flipped_map.begin(); it != flipped_map.end(); ++it) {
			int binkey = it->first;
			double coefvalue = it->second;
			if(abs(coefvalue) > 1e-20){
				t_coefs_aux[nb_confs_p + bf] = coefvalue / sqrt(dS); 
				bin_confs_aux[nb_confs_p + bf] = binkey;
				bf++;
			}
		}
                /*
		cout << "output of map (coefs)" << endl;
		print_array_int(bin_confs_aux, nb_confs_p + bf );
		print_array(t_coefs_aux, nb_confs_p + bf );
		cout << " " << endl;
                 */
		
		flipped_map.clear();
		//delete[] flipped_states;
		
		// check this!!
		for(int bci = 0; bci < nb_confs_p + bf; bci++){
			bin_confs[nb_child + beg_states_sum + bci] = bin_confs_aux[bci];
			t_coefs[nb_child + beg_states_sum + bci] = t_coefs_aux[bci];
			}
		/*	
		memcpy(bin_confs + nb_child + beg_states_sum, bin_confs_aux , sizeof(int)*(nb_confs_p + bf) );
		memcpy(t_coefs + nb_child + beg_states_sum, t_coefs_aux , sizeof(double)*(nb_confs_p + bf) );	
		*/
		
		delete[] bin_confs_aux;
		delete[] t_coefs_aux;
		
		np_binconfs[np_child + p] = (nb_confs_p + bf);
		sectorchild.nb+= (nb_confs_p + bf);
		beg_states_sum+= (nb_confs_p + bf);
		nbtotal+= (nb_confs_p + bf);
		}
	

	//memcpy(np_binconfs+np_child, np_binconfs_aux , sizeof(int)*np_parent);

	
	delete[] (sectorchild.bin_confs);
	delete[] (sectorchild.t_coefs);
	delete[] (sectorchild.np_binconfs);
	
	(sectorchild.bin_confs) = new int[nb_child + nbtotal];
	(sectorchild.t_coefs) = new double[nb_child + nbtotal];
	(sectorchild.np_binconfs) = new int[np_child + np_parent];
	
	for(int p = 0; p < np_child + np_parent; p++){
		(sectorchild.np_binconfs)[p] = np_binconfs[p];
		}
	for(int b = 0; b < nb_child + nbtotal; b++){
		(sectorchild.bin_confs)[b] = bin_confs[b];
		(sectorchild.t_coefs)[b] = t_coefs[b];
		}
	
	
	delete[] bin_confs;
	delete[] t_coefs;
	delete[] np_binconfs;
	
	beg_states.clear();	
	
	
}





void O_EAST(bin_sector &sectorchild, bin_sector sectorparent){
	
	int nsites_curr = iter();
	
	int np_parent = sectorparent.np; 
	int nb_parent = sectorparent.nb;
	
	int np_child = sectorchild.np;
	int nb_child = sectorchild.nb;
	
	int *bin_confs = new int[nb_child + nb_parent];
	double *t_coefs = new double[nb_child + nb_parent];
	int *np_binconfs = new int[np_child + np_parent];
	
	
	for(int i = 0; i < nb_child; i++){
		bin_confs[i] = sectorchild.bin_confs[i];
		t_coefs[i] = sectorchild.t_coefs[i];
		}
	for(int p = 0; p < np_child; p++){
		np_binconfs[p] = sectorchild.np_binconfs[p];
		}
	
	/* I decided to avoid memcpy
	
	if(nb_child != 0){
		memcpy(bin_confs, (sectorchild.bin_confs) , sizeof(int)*nb_child);
		memcpy(t_coefs, (sectorchild.t_coefs) , sizeof(double)*nb_child);	
	}
	if(np_child != 0){
		memcpy(np_binconfs, (sectorchild.np_binconfs) , sizeof(int)*np_child);
	}
	* */
	
	
	for(int b = 0; b < nb_parent; b++){
		double tpb = sectorparent.t_coefs[b];
		int binary_index_child = pow(2, 2*nsites_curr-1) + sectorparent.bin_confs[b];
		t_coefs[nb_child + b] = tpb;
		bin_confs[nb_child + b] = binary_index_child;
		}
	for(int p = 0; p < np_parent; p++){
		np_binconfs[np_child + p] = sectorparent.np_binconfs[p];
		}
	
	delete[] (sectorchild.bin_confs);
	delete[] (sectorchild.t_coefs);
	delete[] (sectorchild.np_binconfs);
	
	
	(sectorchild.bin_confs) = bin_confs;
	(sectorchild.t_coefs) = t_coefs;
	(sectorchild.np_binconfs) = np_binconfs;

	sectorchild.nb += nb_parent;

}
