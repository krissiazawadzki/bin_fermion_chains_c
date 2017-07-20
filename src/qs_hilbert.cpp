#include "qs_hilbert.h"


#include <iostream>
#include <cstdlib>
#include <cfloat>
#include <cstring>

/* auxiliary packages */
#include "matrix_ops/matrix_ops.h"


/* First iteration */
int iteration = -1;

/* maximum value of charge reached at previous iteration */
int qmax_previous = 0;
/* maximum value of charge reached at current iteration*/
int qmax_current = 0;

/* vector of maximum spin associated with charge q at previous iteration*/
int dsmax_previous[2*QMAX+1];
/* vector of maximum spin associated with charge q at current iteration*/
int dsmax_current[2*QMAX+1];

void set_iter(int n)
/* 	void set_iter(int n)
 * 
 *	Function updtates current iteration cycle.
 * 
 * 	Inputs:
 * 		int n - index of the current iteration
 * 
 *  Outputs:
 * 		none
 * */
{
	iteration = n;
}

int iter(){
/*	int iter()
 * 
 * 	Function returns index of the current iteration.
 * 
 * 	Inputs:
 * 		none
 * 
 * 	Outputs:
 * 		int iteration - current iteration index
 * 
 * */
	return iteration;
}



void error(string message)
/*	void error(string message)
 * 
 * 	Function shows a error message when errors are detected.
 * 
 * 	Inputs:
 * 		string message - the message to be shown
 * 
 * 	Outputs:
 * 		none
 * 
 * */
{
  cerr << "***" << message.c_str();
  exit(1);
}


int q_index(int q)
/*	int q_index(int q)	
 * 	
 * 	Function converts charge value to matrix index QS_matrix.
 * 
 * 	Inputs:
 * 		int q - charge value
 * 
 * 	Outputs:
 * 		int q_index - converted index in QS_matrix
 * */
{
  return QMAX - q;
}

int index_q(int index)
/*	int index_q(int index)
 * 
 * 	Function converts from charge index in QS_matrix to charge value.
 * 
 * 	Inputs:
 * 		int index - index in QS_matrix
 * 
 * 	Outputs:
 * 		int index_q - charge value corresponding to index in QS_matrix
 * 
 * */
{
  return QMAX - index;
}


int ds_index(int ds)
/*	int q_index(int ds)	
 * 	
 * 	Function converts spin value to matrix index QS_matrix.
 * 
 * 	Inputs:
 * 		int ds - spin value
 * 
 * 	Outputs:
 * 		int ds_index - converted index in QS_matrix
 * */
{
  return DSMAX - ds;
}

int index_ds(int index)
/*	int index_q(int index)
 * 
 * 	Function converts from charge index in QS_matrix to charge value.
 * 
 * 	Inputs:
 * 		int index - index in QS_matrix
 * 
 * 	Outputs:
 * 		int index_q - charge value corresponding to index in QS_matrix
 * 
 * */
{
  return DSMAX - index;
}


int qsmax()
/*	int qmax()
 * 
 * 	Function returns maximum value of charge allowed in whole program.
 * 
 * 	Inputs:
 * 		none
 * 
 * 	Outputs:
 * 		int qmax - maximum charge value allowed
 * */
{
  return QMAX;
}

int dsmax()
/*	int dsmax()
 * 
 * 	Function returns maximum value of spin (2*S) allowed in whole program.
 * 
 * 	Inputs:
 * 		none
 * 
 * 	Outputs:
 * 		int dsmax - maximum spin value allowed (2*S)
 * 
 * */
{
  return DSMAX;
}




void set_qmax_prev(int q)
/*	void set_qmax_prev(int q)
 * 
 * 	Function deposits new maximum charge value reached at previous iteration.
 * 
 * 	Inputs:
 * 		int q - value of charge
 * 
 * 	Outputs:
 * 		none
 * 
 * */
{
  qmax_previous = q;
}

int qmax_prev()
/*	int qmax_prev()
 * 
 * 	Function returns maximum charge at previous iteration.
 * 
 * 	Inputs:
 * 		none
 * 
 * 	Outputs:
 * 		int qmax_prev - maximum charge at previous iteration
 * 
 * */
{
  return qmax_previous;
}

void set_qmax_curr(int q)
/*	void set_qmax_curr(int q)
 * 
 *	Function deposits new qmax_current.
 * 
 * 	Inputs:
 * 		int q - new charge value
 * 
 * 	Ouputs:
 * 		none
 * */
{
  qmax_current = q;
}

int qmax_curr()
/*	int qmax_curr()	
 * 
 * 	Function returns maximum charge at previous iteration.
 * 
 * 	Inputs:
 * 		none
 * 
 * 	Outputs:
 * 		int qmax_current - maximum charge value at current iteration
 * 
 * */
{
  return qmax_current;
}



void set_dsmax_prev(int q, int dsmax)
/*	void set_dsmax_prev(int q, int dsmax)
 * 
 * 	Function defines maximum spin for charge q at previous iteration.
 * 
 * 	Inputs:
 * 		int q - charge 
 * 		int dsmax - maximun spin for q at previous iteration
 * 
 * 	Outputs:
 * 		none
 * */
{
  dsmax_previous[QMAX-q] = dsmax;
}



int dsmax_prev(int q)
/*	int dsmax_prev(int q)	
 * 
 *	Function returns maximum spin at previous iteration for charge q.
 * 
 * 	Inputs:
 * 		int q - charge
 * 
 * 	Outputs:
 * 		int dsmax_prev - maximum spin at previous iteration for q
 * */
{
  return dsmax_previous[QMAX-q];
}

void set_dsmax_curr(int q, int dsmax)
/*	void set_dsmax_curr(int q, int dsmax)
 * 
 * 	Function defines maximum spin for charge q at current iteration.
 * 
 * 	Inputs:
 * 		int q - charge 
 * 		int dsmax - maximun spin for q at current iteration
 * 
 * 	Outputs:
 * 		none
 * */
{
  dsmax_current[QMAX-q] = dsmax;
}

int dsmax_curr(int q)
/*	int dsmax_curr(int q)	
 * 
 *	Function returns maximum spin at current iteration for charge q.
 * 
 * 	Inputs:
 * 		int q - charge
 * 
 * 	Outputs:
 * 		int dsmax_curr - maximum spin at current iteration for q
 * */
{
  return dsmax_current[QMAX-q];
}

void update_q_dS()
/*	void update_q_dS()
 * 
 * 	Function updates maximum charge and spin values from previous to 
 * 	current iteration. We also set ground state energy as DBL_MAX to 
 * 	compute next iteration ground state energy.
 * 
 * 	Inputs:
 * 		none
 * 
 * 	Outputs:
 * 		none 
 * */
{
	
	
	if(qmax_prev() < QMAX){
	  set_qmax_prev(qmax_curr()); //update qmax_prev only if it does not reach QMAX
	}
	
	if(qmax_curr() < QMAX)
	  set_qmax_curr(qmax_curr()+1); // update qmax only if does not reach QMAX  
	

	for(int q = qmax_curr(); q >= -qmax_curr(); q--){
	  set_dsmax_prev(q,dsmax_curr(q));
	}
	
    int qmax_previous = qmax_prev();	
	for(int q = qmax_previous; q >= -qmax_previous; q--){
		if(dsmax_prev(q)+1 <= DSMAX){ // check if this could be regarded to 2*DSMAX+1
			set_dsmax_curr(q,dsmax_prev(q)+1);
		}
		else
			set_dsmax_curr(q,dsmax_prev(q)-1);
	}	
	
	  
}





int indhamex(int lin, int col)
/*	int indham(int lin, int col)
 *
 * 	Function converts line and column indexes from a bidiomensional array
 *  representing the lower triangular matrix of Hamiltonian (H[0][0],...,
 *  H[lin][col],..., H[N][N]) in a index of linear array (A[0],..., A[i],
 *  ,...A[N*N]).
 *
 * 	Inputs:
 * 		int lin - line index in a matrix H[line][column]
 * 		int col - column index in a matrix H[line][column]
 *
 *  Ouputs:
 * 		int i - converted linear index
 *
 * */
{
  if(col > lin){
    int temp = col;
    col = lin;
    lin = temp;
  }
  int i = lin * (lin+1) / 2 + col;
  return i;
}

