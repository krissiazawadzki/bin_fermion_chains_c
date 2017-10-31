#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>

#include "io_bin.h"
#include "spin_configurations.h"
#include "qs_hilbert.h"

#include "translation_eigenstates_of_S.h"

/* custom aux libraries */
#include "int_array_util.h"

/* auxiliary packages */
#include "matrix_ops/matrix_ops.h"
#include "file_ops/file_ops.h"
#include "file_ops/mkdir.h"
#include "string_ops/string_ops.h"


using namespace std;

void print_all_sectors()
/*	void print_all_sectors()
 *
 * 	Function prints all sectors (Q,S) in current iteration. 
 *
 * 	Inputs:
 * 		none
 * 	Outputs:
 * 		none
 *
 * */
{ 
	if(iter() == 1){
		for(unsigned int qind = 0; qind < QS_matrix().size() ; qind++){
			for(unsigned int ds = 0; ds < QS_matrix()[0].size(); ds++){
				if(QS_matrix()[qind][ds].nub == 0){
					continue;
				}
					bin_sector sect = QS_matrix()[qind][ds];
					print_sector(sect, ds);
				cout << "******************************************************" << endl;
			}
		}
	}
	else{	
		int qmax_current = qmax_curr();
		for(int q = qmax_current; q >= -qmax_current; q--){
			int qind = q_index(q);		// insures indices are non-negative
			if(qind >= 0 && qind < 2*QMAX+1){
				int dsmax_c = dsmax_curr(q);
				for(int ds = dsmax_c%2; ds <= dsmax_c; ds +=2){
					if(QS_matrix()[qind][ds].np == 0){
						continue;
					}
					bin_sector sect = QS_matrix()[qind][ds];
					
					
					// constructing translation matrix
					
					int nB = 0, nP = 0;
					double **Tmat = translation_matrix(sect, nP, nB);
					print_sector(sect, ds);
					std::cout << "translation matrix" << std::endl;
					print_matrix(Tmat, nP, nB);
					std::cout << "" << std::endl;
					delete_matrix(Tmat);
					
				}
			cout << "******************************************************" << endl;
			}
		}			
	}
}




void print_target_hilbert_space(int nsites, int Q, int dS, bool print_bin_confs, std::string folder_results, int wfmt, int precisionfmt, bool print_matrix_in_terminal)
/*	void print_target_hilbert_space
 *
 * 	Function prints specific information of the reduced Hilbert space (Q,S) for a lattice with nsites 
 * 	Results are stored in folder_results
 *
 * 	Inputs:
 * 		none
 * 	Outputs:
 * 		none
 *
 * */
{ 
	
	// first, check if the Hilbert space requested exists
	bool requested_QS_ok = false;
	
	int qmax_current = qmax_curr();
	for(int q = qmax_current; q >= -qmax_current; q--){
		int qind = q_index(q);	
		if(qind >= 0 && qind < 2*QMAX+1){
			int dsmax_c = dsmax_curr(q);
			for(int ds = dsmax_c%2; ds <= dsmax_c; ds +=2){
				if(q == Q && ds == dS){
					requested_QS_ok = true;
				}		
			}
		}
	}			
	if(requested_QS_ok){
		int qind = q_index(Q);	
		bin_sector sec = QS_matrix()[qind][dS];
		
		
		int nb = sec.nb;
		int *inv_u_bin_confs = new int[nb];
		for(int i = 0; i < nb; i++){
			inv_u_bin_confs[i] = sec.bin_confs[i];
		}
		
		int *u_bin_confs;
		int nBU = find_unique_intarray(inv_u_bin_confs, nb, u_bin_confs);
		
		std::cout << "creating " << folder_results << " ... " << (makePath(folder_results) ? "OK" : "failed") << std::endl;
  
		std::string file_bin_confs = folder_results+"/binary_ids_in_Hilbert_nsites="+num2str(iter())+"_Nf="+num2str(Q)+"_dS="+num2str(dS)+".txt";
		std::string file_rotation = folder_results+"/rotation_matrix_nsites="+num2str(iter())+"_Nf="+num2str(nsites+Q)+"_dS="+num2str(dS)+".txt";
		
		print_array_to_file_cpp(file_bin_confs, u_bin_confs, nBU, "w", wfmt, precisionfmt);
		
		std::cout << "\n binconfs in sector :\t";
		for(int i = 0; i < nBU; i++){
			std::cout << u_bin_confs[i] << "\t";
			}
		std::cout << std::endl;

		if(print_bin_confs){
			std::cout << "binary confs inside reduced Hilbert space" << std::endl;
			for(int i = 0; i < nBU; i++){
				int q_bin;
				int nsites = iter();
				int *fns_bin = new int[2*nsites];
				int index_bin = u_bin_confs[i];
				int *binstate = index_to_binary(index_bin, nsites, fns_bin, q_bin, 0);
				cout << print_conf(binstate, nsites) << endl;
				delete[] binstate;
				delete[] fns_bin;	
				}
			std::cout << std::endl;
		}
		delete[] u_bin_confs;
		delete[] inv_u_bin_confs;

		std::cout << "\n nP: " << sec.np << "\t nP: " << nBU << std::endl;
	
		int nB = 0, nP = 0;
		double **Tmat = translation_matrix(sec, nP, nB);
		
		print_matrix_to_file_scientific_cpp(file_rotation, Tmat, nP, nB, "w", 5, 15);
		
		std::cout << std::endl;
		std::cout << "Rotation (Q,Sz)-> (Q,S)" << std::endl;

		if(print_matrix_in_terminal){
			print_matrix(Tmat, nP, nB);
		}
		std::cout << "" << std::endl;
		
		// saving results in 
		
		delete_matrix(Tmat);
		
	}
	
}




void print_sector_friendly(bin_sector sec)
/* void print_sector(sector sec)
 * 
 * 	Function print sector sec fields.
 * 
 * 	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		none
 * 
 * */
{
	cout << " " << endl;
	int dSz = sec.dS;
	cout << "Q = " << sec.Q << ", dS = " << sec.dS << ", dSz = " << dSz << endl;
	cout << "np = " << sec.np << ", nb = " << sec.nb << ", nub = " << sec.nub << endl;

	
	if((sec.npv).size() == 4)
		cout << "np_N = " << npN(sec) << ", np_E = " << npE(sec) << ", np_S = " << npS(sec) << ", np_W = " << npW(sec) << endl;
	
	cout << "configurations belonging to sector : ";

	for(int i = 0; i < sec.nb; i++){
		cout << sec.bin_confs[i] << " " ;
	}
	cout << endl;
	cout << "number of configurations per state : ";
	for(int i = 0; i < sec.np; i++){
		cout << sec.np_binconfs[i] << " " ;
	}
	cout << endl;
	cout << "translation coeficients of states : ";
	for(int i = 0; i < sec.nb; i++){
		cout << sec.t_coefs[i] << " " ;
	}
	cout << endl;
	
	int accbin = 0;
	for(int p = 0; p < sec.np; p++){
		cout << "state " << p << " = ";
		for(int b = 0; b < sec.np_binconfs[p]; b++){
			int q_bin;
			int nsites = iter();
			int *fns_bin = new int[2*nsites];
			int index_bin = sec.bin_confs[accbin];
			int *binstate = index_to_binary(index_bin, nsites, fns_bin, q_bin, 0);
			cout << sec.t_coefs[accbin];
			string state_str = "|";
			for(int site=0; site < nsites; site++){
				if(binstate[2*site] == 1 && binstate[2*site+1]==0)
					state_str+="u";
				if(binstate[2*site] == 1 && binstate[2*site+1] ==1)
					state_str+="c";
				if(binstate[2*site]==0 && binstate[2*site+1]==1)
					state_str+="d";
				if(binstate[2*site]==0 && binstate[2*site+1]==0)
					state_str+="0";
				}
				state_str+="> ";
			
			cout << state_str << endl;
			//cout << sec.t_coefs[accbin] << " " << print_conf(binstate, nsites) << endl;
			delete[] binstate;
			delete[] fns_bin;	
			accbin++; 
			}
			
		
		}
	
	

}

void print_all_sectors_friendly()
/*	void print_all_sectors()
 *
 * 	Function prints all sectors (Q,S) in current iteration. 
 *
 * 	Inputs:
 * 		none
 * 	Outputs:
 * 		none
 *
 * */
{ 
	if(iter() == 1){
		for(unsigned int qind = 0; qind < QS_matrix().size() ; qind++){
			for(unsigned int ds = 0; ds < QS_matrix()[0].size(); ds++){
				if(QS_matrix()[qind][ds].nub == 0){
					continue;
				}
					bin_sector sect = QS_matrix()[qind][ds];
					print_sector(sect, ds);
				cout << "******************************************************" << endl;
			}
		}
	}
	else{	
		int qmax_current = qmax_curr();
		for(int q = qmax_current; q >= -qmax_current; q--){
			int qind = q_index(q);		// insures indices are non-negative
			if(qind >= 0 && qind < 2*QMAX+1){
				int dsmax_c = dsmax_curr(q);
				for(int ds = dsmax_c%2; ds <= dsmax_c; ds +=2){
					if(QS_matrix()[qind][ds].nub == 0){
						continue;
					}
					bin_sector sect = QS_matrix()[qind][ds];
					print_sector_friendly(sect);
				}
			cout << "******************************************************" << endl;
			}
		}			
	}
}
