#include "QSmatrix.h"
#include "qs_hilbert.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <algorithm>

vector<vector<bin_sector> > QSmatrix(2*QMAX+1,vector<bin_sector>(DSMAX+1));	

void init_QSmatrix()
/*	void init_QS_matrix()
 * 
 * 	Function initializes the matrix of sectors with all elements cleaned.
 * 
 * 	Inputs:
 * 		none
 * 
 * 	Outputs:
 * 		none
 * 
 *  */
{
  for(int i = 0; i < 2*QMAX+1; i++){
    for(int j = 0; j < DSMAX+1; j++){
      clear_sector(QSmatrix[i][j]);
    }
  }
}

void set_QSmatrix(const vector<vector<bin_sector> >& QScurrent)
/*	void set_QSmatrix(const vector<vector<sector> >& QScurrent)
 * 
 * 	Function updates QSmatrix content according to QScurrent matrix.
 * 
 * 	Inputs:
 * 		const vector<vector<sector> >& QScurrent - the new matrix 
 * 			we want in QSmatrix
 * 	Outputs:
 * 		none
 *  */
{
  QSmatrix = QScurrent;
}

void update_QSmatrix()
/*	void update_QSmatrix()
 * 
 * 	Function updates QSmatrix cleaning sectors of the previous iteration.
 * 
 * 	Inputs:
 * 		none
 * 	Outputs:
 * 		none
 * 
 * */
{
	int qmax_previous = qmax_prev();
	for(int q = qmax_previous; q >= -qmax_previous; q--){
		int qind = q_index(q);		// makes indices non-negative
		if(qind >= 0 && qind < 2*QMAX+1){
			int dsmax_p = dsmax_prev(q);
			for(int ds = dsmax_p%2; ds <= dsmax_p; ds +=2){
				clear_sector(QSmatrix[qind][ds], true);
			}
		}			
	}
}



void print_QSmatrix(const vector< vector<bin_sector> > &qs_matrix)
/*	void print_QSmatrix(const vector< vector<sector> > &qs_matrix)
 * 
 * 	Function prints the number of states in each sector of QSmatrix
 * 	if this number is non-zero.
 * 
 * 	Inputs:
 * 		const vector< vector<sector> > &qs_matrix - QSmatrix to be 
 * 			printed
 * 	Outputs:
 * 		none
 * 
 *  */
{	
	

  for(unsigned int i = 0; i < qs_matrix.size() ; i++){
    for(unsigned int j = 0; j < qs_matrix[0].size(); j++){
      if (QSmatrix[i][j].np != 0)
	cout << setw(12) << qs_matrix[i][j].np;
      else 
	cout << setw(12) << 0.0;
    }
    cout << "\n";      
  }
	
	  
}

void clear_QSmatrix(vector< vector<bin_sector> > &qs_matrix)
/*	void clear_QSmatrix(vector< vector<sector> > &qs_matrix)
 * 
 * 	Function erases all elements in QSmatrix.
 * 
 * 	Inputs:
 * 		vector< vector<sector> > &qs_matrix - the matrix we want to clear
 * 	Outputs:
 * 		none
 * 
 * */
{	
  for(unsigned int i = 0; i < qs_matrix.size() ; i++){
    qs_matrix[i].clear();
  }
  qs_matrix.clear();
}



int qs_nub(int q, int ds)
// returns no of binary confs in (q,ds)
{
  int qind = q_index(q);
  bin_sector qs = QSmatrix[qind][ds];
  return qs.nub;
}

int qs_np(int q, int ds)
// returns dimension of (q,ds)
{
  int qind = q_index(q);
  bin_sector qs = QSmatrix[qind][ds];
  return qs.np;
}

int qs_nb(int q, int ds)
// returns dimension of (q,ds)
{
  int qind = q_index(q);
  bin_sector qs = QSmatrix[qind][ds];
  return qs.nb;
}



int qs_nGender(int q, int ds, int g)
// returns number of states of gender g
{
  int qind = q_index(q);
  bin_sector qs = QSmatrix[qind][ds];

  return qs.npv[g];
}

int qs_ranking10gender(int q, int ds, int p)
// returns ranking of primitive state p among those of its gender in sector (q, s)
// and gender of p, encrypted as 10*ranking+gender
{
// primitive states are stored in S, E, W, N order

  int qind = q_index(q);
  bin_sector qs = QSmatrix[qind][ds];

  int ret = p - qs_np(q, ds);
  if(ret >= 0){
    cout << "p=" << p << " exceeds np = " << qs_np(q, ds) << " in qs_gender" 
	 << endl
	 << "    " << qs.npv[NORTH] << endl
	 << qs.npv[WEST] << "         " << qs.npv[EAST] << endl
	 << "    " << qs.npv[SOUTH] << endl;
    exit(1);
  }

  ret += qs.npv[NORTH];
  if(ret >= 0)
    return 10 * ret + NORTH;

  ret += qs.npv[WEST];
  if(ret >= 0)
    return 10 * ret + WEST;

  ret += qs.npv[EAST];
  if(ret >= 0)
    return 10 * ret + EAST;

  ret += qs.npv[SOUTH];
  if(ret >= 0)
    return 10 * ret + SOUTH;
  return -1;			// keeps compiler happy
}

int qs_gender(int q, int ds, int p)
// returns gender of primitive state p
{
  int ret = qs_ranking10gender(q, ds, p);
  
  return ret % 10;		// ret = 10*(offset position) + gender
}

int qs_ranking(int q, int ds, int p)
// returns ranking of primitive state among those of same gender
// eg, returns 0 if p is first WEST state or first NORTH state
{
  int ret = qs_ranking10gender(q, ds, p);

  return ret / 10;  		// ret = 10*(offset position) + gender
}
    
