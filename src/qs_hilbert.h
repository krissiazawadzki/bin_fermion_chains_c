#ifndef __QS_HILBERT__

#include <string>

#define NORTH 0
#define EAST 1
#define SOUTH 2
#define WEST 3

using namespace std;


const int QMAX = 9;
const int DSMAX = 9;


void set_iter(int n);
int iter();
void set_itermax(int itermax);
int iter_max();

void error(string message);

int q_index(int q);
int index_q(int index);
int ds_index(int ds);
int index_ds(int index);

void set_qmax_prev(int q);
int qmax_prev();
void set_qmax_curr(int q);
int qmax_curr();

void set_dsmax_prev(int q, int dsmax);
int dsmax_prev(int q);
void set_dsmax_curr(int q, int dsmax);
int dsmax_curr(int q);


void update_q_dS();


int indhamex(int lin, int col);

int iter_max();
void set_itermax(int itermaxarg);

#define __QS_HILBERT__
#endif
