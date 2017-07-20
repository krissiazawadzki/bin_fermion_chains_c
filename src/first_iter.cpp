#include "first_iter.h"
#include "qs_hilbert.h"
#include "QSmatrix.h"

#include <cmath>

void generate_single_site_hilbert_space()
/*	void generate_single_site_hilbert_space()
 * 
 * 	Function initializes iteration cycle by N = 1, adding to QSmatrix
 *  the first three sectors [Q=-1,S=0], [Q=0,S=1] and [Q=1,S=0] filled 
 *  states:
 * 		|0> = empty
 * 
 * 		|1> = up
 * 
 * 		|2> = up down 
 * 
 * 	Inputs:
 * 		none
 * 	Outputs:	
 * 		none
 * */	
{

  set_iter(1);
 
  // creating sector Q = 1, S = 0, np = 1 
  QS_matrix()[q_index(1)][0].Q = 1;
  QS_matrix()[q_index(1)][0].dS = 0;
  QS_matrix()[q_index(1)][0].nb = 1;
  QS_matrix()[q_index(1)][0].np = 1;
  QS_matrix()[q_index(1)][0].nub = 1;
  QS_matrix()[q_index(1)][0].bin_confs = new int[1];
  QS_matrix()[q_index(1)][0].bin_confs[0] = 3;
  QS_matrix()[q_index(1)][0].t_coefs = new double[1];
  QS_matrix()[q_index(1)][0].t_coefs[0] = 1.0; 
  QS_matrix()[q_index(1)][0].np_binconfs = new int[1];
  QS_matrix()[q_index(1)][0].np_binconfs[0] = 1;
  
  // creating sector Q = 0, S = 1/2, np = 1
  QS_matrix()[q_index(0)][1].Q = 0;
  QS_matrix()[q_index(0)][1].dS = 1;
  QS_matrix()[q_index(0)][1].nb = 1;
  QS_matrix()[q_index(0)][1].np = 1;
  QS_matrix()[q_index(0)][1].nub = 1;
  QS_matrix()[q_index(0)][1].bin_confs = new int[1];
  QS_matrix()[q_index(0)][1].bin_confs[0] = 2;
  QS_matrix()[q_index(0)][1].t_coefs = new double[1];
  QS_matrix()[q_index(0)][1].t_coefs[0] = 1.0;
  QS_matrix()[q_index(0)][1].np_binconfs = new int[1];
  QS_matrix()[q_index(0)][1].np_binconfs[0] = 1;
 

  // creating sector Q = -1, S = 0, nr = 1, np = 0, E = 0
  QS_matrix()[q_index(-1)][0].Q = -1;
  QS_matrix()[q_index(-1)][0].dS = 0;
  QS_matrix()[q_index(-1)][0].nb = 1;
  QS_matrix()[q_index(-1)][0].np = 1;  
  QS_matrix()[q_index(-1)][0].nub = 1;  
  QS_matrix()[q_index(-1)][0].bin_confs = new int[1];
  QS_matrix()[q_index(-1)][0].bin_confs[0] = 0; 
  QS_matrix()[q_index(-1)][0].t_coefs = new double[1];
  QS_matrix()[q_index(-1)][0].t_coefs[0] = 1.0;
  QS_matrix()[q_index(-1)][0].np_binconfs = new int[1];
  QS_matrix()[q_index(-1)][0].np_binconfs[0] = 1;
  
  
  
  set_qmax_prev(1); // sets the max charge of iteration N = 1, single site states
  set_qmax_curr(qmax_prev()+1);

}


void update_from_single_site_to_dimer(){
	
  set_iter(2);
  
  set_dsmax_prev(1,0);
  set_dsmax_prev(0,1);
  set_dsmax_prev(-1,0);
  
  set_dsmax_curr(2,0);
  set_dsmax_curr(1,1);
  set_dsmax_curr(0,2);
  set_dsmax_curr(-1,1);
  set_dsmax_curr(-2,0);
 
	
	
}
