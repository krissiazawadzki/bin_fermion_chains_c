#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>

#include "io_bin.h"
#include "spin_configurations.h"
#include "qs_hilbert.h"

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
					print_sector(sect, ds);
				}
			cout << "******************************************************" << endl;
			}
		}			
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
