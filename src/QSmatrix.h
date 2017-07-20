#ifndef __QSMATRIX__

#include "bin_sectors.h"

extern vector<vector<bin_sector> > QSmatrix;	

// function initializes all elements of QS matrix with empty sectors
void init_QSmatrix();

// function sets QSmatrix with matrix elements of the previous and current iteractions
void set_QSmatrix(const vector<vector<bin_sector> >& QScurrent);

// function returns QSmatrix 
inline vector<vector<bin_sector> >& QS_matrix(){
	return QSmatrix;
}

// function deletes elements from previous iteration
void update_QSmatrix();

// function prints QSmatrix
void print_QSmatrix(const vector< vector<bin_sector> > &qs_matrix);

// Function liberates matrix elements
void clear_QSmatrix(vector< vector<bin_sector> > &qs_matrix);


void qs_scan(int q, int ds);
int qs_nub(int q, int ds);
int qs_np(int q, int ds);
int qs_nb(int q, int ds);
int qs_nGender(int q, int ds, int g);
int qs_ranking10gender(int q, int ds, int p);
int qs_gender(int q, int ds, int p);
int qs_ranking(int q, int ds, int p);

#define __QSMATRIX__
#endif






