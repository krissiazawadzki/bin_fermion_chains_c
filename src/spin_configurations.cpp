#include "spin_configurations.h"
#include "matrix_ops.h"


#include <cstdlib>
#include <cfloat>
#include <cstdio>
#include <sstream>
#include <cmath>

using namespace std;


/*	
 * 	Function returns the number of possible configurations for a given 
 * 	spin chain of size nsites 
 * 
 * unsigned long 
 * */
int chain_possible_configurations(int nsites){
		return pow(4.0,nsites);
}

/* *index_to_binary(int index, int nsites, int *occupied_sites, int &charge, int relativebegin)
 *
 * 	Function returns the binary version of states
 *
 * 	Inputs:
 *      index - integer labelling the state
 * 	nsites - number of sites of the configuration
 *      occupied_sites - pointer to array to be filled with occupied sites
 *      charge - address of a int to store the charge of the state
 *      relative_begin - if the state belongs to a piece of the chain, relative begin is used
 *
 * 	Outputs:
 * 	state - binary version of state
 *
 * */
int *index_to_binary(int index, int nsites, int *occupied_sites, int &charge, int relativebegin){
  int *binstate = new int[2*nsites];
  int num_occ = 0;
  for (int i = 0; i < 2*nsites; i++){
    int digit = index / int(pow(2,i)) % 2;
    binstate[2*nsites-i-1] = int(digit);
    if(digit == 1){
      charge+=1;
      occupied_sites[num_occ] = 2*nsites-1-i+relativebegin;
      num_occ++;
    }
  }
  charge = num_occ;
  return binstate;
}

int binary_to_index(int *state, int nsites){
  int indexconf = 0;
  for(int i = 0; i < 2*nsites; i++){
    int power = pow(2,2*nsites-i-1);
    indexconf+= state[i] * power;
  }
  return indexconf;
}	


/* int Q_state(int *state, int nsites)
 *
 * 	Function returns the total charge for a given state
 *
 * 	Inputs:
 * 	state - binary version of state
 * 	nsites - number of sites of the configuration
 *
 * 	Outputs:
 * 	charge - value of Q
 *
 * */
int Q_state(int *state, int nsites){
		int Q = 0;
		for(int i = 0; i < 2*nsites; i++){
                  Q+= state[i];
                }
  return Q - nsites ;
}


/* int Sz_state(int *state, int nsites)
 * 
 * 	Function returns the total spin in z direction for a given state 
 * 	result is 2*Sz
 * 
 * 	Inputs:
 * 	state - binary version of state
 * 	nsites - number of sites of the configuration
 * 
 * 	Outputs:
 * 	spin z - value of Sz
 * 
 * */
int Sz_state(int *state, int nsites){
	int Sz = 0;
		for(int i = 0; i < nsites; i++){
			int ups = state[2*i];
			int downs = state[2*i+1];
			Sz+= ups - downs;
			}	
	return Sz;
}




string print_conf(int *state, int nsites){
	string conf = " | \t";
    for(int i = 0 ; i < 2*nsites; i++){
      conf+= static_cast<ostringstream*>( &(ostringstream() << state[i]) )->str()+ "\t";
    }
    conf+=">";
    return conf;    	
}




int* find_Q_charge(int* Qs, int nsites, int Q, int &numQs){
	numQs = 0;
	int total_confs = chain_possible_configurations(nsites);
	int *Qref = zero_matrix_int(1,total_confs)[0];
	for(int i = 0; i < total_confs; i++){
		if(Qs[i] == Q){
			numQs++;
			Qref[i] = 1;
			}
		}
	int *Qs_withQ = new_matrix_int(1,numQs)[0];
	int qi = 0;
	for(int i = 0; i < total_confs; i++){
		if(Qref[i] == 1){
			Qs_withQ[qi] = i;
			qi++;
			}
		}
	delete []Qref;
	return Qs_withQ;	
}

int* find_QSz_states(int nsites, int Q, int Sz, int &numQSzs, int *invQSz_idx){
	numQSzs = 0;
	int total_confs = chain_possible_configurations(nsites);
	int *QSzref = new int[total_confs];
	for(int i = 0; i < total_confs; i++){
		int Qi;
		int *occ_sites = new int[2*nsites];
		int *statei = index_to_binary(i, nsites, occ_sites, Qi, 0);
		int Szi = Sz_state(statei, nsites);
		invQSz_idx[i] = -1;
		if((Qi == Q) && (Szi == Sz)){
			QSzref[i] = 1;
			invQSz_idx[i] =numQSzs;
			numQSzs++; 
			}
		delete[] statei;
          delete[] occ_sites;
		}
	int *QSzs = new int[numQSzs];
	int qszi = 0;
	for(int i = 0; i < total_confs; i++){
		if(QSzref[i] == 1){
			QSzs[qszi] = i;
			qszi++;
			}
		}
	delete []QSzref;
	return QSzs;	
}
