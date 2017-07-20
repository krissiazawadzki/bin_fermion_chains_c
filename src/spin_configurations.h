#ifndef __SPIN_CONFIGURATIONS__

#include <cstring>
#include <iostream>

// returns the number of possible configurations for nsites
int chain_possible_configurations(int nsites);

// converts index of configuration to the binary version of this state
int *index_to_binary(int index, int nsites, int *occupied_sites, int &charge, int relativebegin=0);

// classify state according to the binary configuration it represents
int binary_to_index(int *state, int nsites);

// returns the charge of a state
int Q_state(int *state, int nsites);

// returns the spin 2Sz of the state
int Sz_state(int *state, int nsites);

std::string print_conf(int *state, int nsites);

//returns the indexes where the configurations have charge Q 
int* find_Q_charge(int* Qs, int nsites, int Q, int &numQs);

// returns the indexes where the configurations have charge Q and spin Sz
int* find_QSz_states(int nsites, int Q, int Sz, int &numQSzs, int *invQSz_idx);


#define __SPIN_CONFIGURATIONS__
#endif
