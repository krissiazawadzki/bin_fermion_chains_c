/* custom aux libraries */
#include "int_array_util.h"

/* auxiliary packages */
#include "matrix_ops/matrix_ops.h"

#include "translation_eigenstates_of_S.h"


double **translation_matrix(bin_sector &sec, int &Tmat_p, int &Tmat_b)
/*
 * 	description: 
 * 		Function generates the translation matrix from eigenstates of S
 * 		operator created iteratively
 * 
 * 	inputs:
 * 		- sec: sector data structure
 * 		- Tmat_p: dimension p of translation matrix, i.e. number of eigenstates of sector
 * 		- Tmat_b: dimension b of translation matrix, i.e. number of binary configurations needed to construct the eigenstates of S
 * 
 * 	outputs:
 * 		- T_pb: translation matrix with the Clebsch-Gordan coefficients
 * 
 * 		- notice that the field nub of sector is changed to have the number of unique binary configurations found
 * */

{
	
	int nP = sec.np;
	int nB = sec.nb;
	
	int *inv_u_bin_confs = new int[nB];
	for(int i = 0; i < nB; i++){
		inv_u_bin_confs[i] = sec.bin_confs[i];
		}
	
	int *u_bin_confs;
	int nBU = find_unique_intarray(inv_u_bin_confs, nB, u_bin_confs);
	
	// inversion map
	map_unique_indexes(sec.bin_confs, nB, u_bin_confs, nBU, inv_u_bin_confs);
	
	Tmat_p = nP;
	Tmat_b = nBU;
	
	double **T_pb = zero_matrix(nP, nBU);
	int startp = 0;
	
	for(int p = 0; p < nP; p++){
		int nb_p = sec.np_binconfs[p];
		
		for(int b = 0; b < nb_p; b++){
			int ib = inv_u_bin_confs[startp + b];
			T_pb[p][ib] += sec.t_coefs[startp + b];
			}
			startp+= nb_p;
		}
	
	sec.nub = nBU;	
	delete[] inv_u_bin_confs;
	
	return T_pb;
}
