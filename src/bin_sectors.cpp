#include "bin_sectors.h"
#include "spin_configurations.h"
#include "qs_hilbert.h"

#include <iostream>
#include <numeric>
#include <cstdlib>

#include <cstdio>


// must re-check this

using namespace std;

/* sector matrix (contains all sectors, from current and previous iterations */
vector<vector<bin_sector> > matrix_sec(2*QMAX+1,vector<bin_sector>(DSMAX+1));	


int Q(bin_sector sec)
/*	int Q(sector sec)
 * 
 * 	Function returns charge of sector sec.
 * 
 * 	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		int Q - charge of sec
 * 
 * */
{
	return sec.Q;
}

int dS(bin_sector sec)
/*	int dS(sector sec)
 * 
 * 	Function returns spin (2*S) of sector sec.
 * 
 * 	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		int dS - spin of sec
 * 
 * */
{
	return sec.dS;
}
	


int np(bin_sector sec)
/*	int np(sector sec)
 * 
 * 	Function returns the number of primitive states in sector sec.
 * 
 *	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		int np - number of primitive states in sec
 * 	
 * */
{ 
	return sum_npv(sec.npv);
}



int nb(bin_sector sec)
/*	int nb(sector sec)
 * 
 * 	Function returns the number of binary configurations in sector sec.
 * 
 *	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		int nb - number of binary configurations in sec
 * 	
 * */
{
	return sec.nb;
}

int nub(bin_sector sec)
/*	int nr(sector sec)
 * 
 * 	Function returns the number of unique binary configurations in sector sec.
 * 
 *	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		int nr - number of binary states in sec
 * 	
 * */
{
	return sec.nub;
}


/*	binary_label_b(bin_sector sec, int b)
 * 
 * 	Function returns the binary index of b-th configuration of sector.
 * 
 *	Inputs:
 * 		sector sec - sector data structure
 * 		int b - b-th configuration inside sector
 * 	Outputs:
 * 		int indexb - index of the b-th configuration in binary notation
 * 	

{
	return sec.bin_confs[b];
}
*  * */


/*	int configuration_b(bin_sector sec, int b)
 * 
 * 	Function returns the binary index of b-th configuration of sector.
 * 
 *	Inputs:
 * 		sector sec - sector data structure
 * 		int b - b-th configuration inside sector
 * 	Outputs:
 * 		int indexbQS - index of the b-th configuration inside sector
 * 	

{
	int binary_index = sec.bin_confs[b];
	return sec.map_bin_labels[binary_index];
}
 * */

vector<int> npv(bin_sector sec)
/*	int npv(sector sec)
 * 
 * 	Function returns the vector containing the number of states came 
 * 	from each direction NORTH, EAST, SOUTH and WEST in npv[0], npv[1], 
 * 	npv[2] and npv[3], respectively. 
 * 
 *	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		vector<int> npv - vector counting number of states in each 
 * 			direction
 * 	
 * */
{
	return sec.npv;
}

int sum_npv(vector<int> npv)
/*	int np(sector sec)
 * 
 * 	Function sums over all primitive states in npv vector.
 * 
 *	Inputs:
 * 		vector<int> npv - vector containing primitive states came 
 * 			from each direction
 * 	Outputs:
 * 		int sum - total number of states counting over all directions.
 * 	
 * */
{
	int sum = accumulate( npv.begin(), npv.end(), 0 );
	return sum;
}

int npN(bin_sector sec)
/* npN(sector sec) 
 * 
 * 	Function returns the number of primitive states came from north.
 * 
 * 	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		int npN - states of sec came from north
 * 
 * */
{
	return sec.npv[NORTH];
}

int npE(bin_sector sec)
/* npE(sector sec) 
 * 
 * 	Function returns the number of primitive states came from east.
 * 
 * 	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		int npE - states of sec came from east
 * 
 * */
{
	return sec.npv[EAST];
}

int npS(bin_sector sec)
/* npS(sector sec) 
 * 
 * 	Function returns the number of primitive states came from south.
 * 
 * 	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		int npS - states of sec came from south
 * 
 * */
{
	return sec.npv[SOUTH];
}

int npW(bin_sector sec)
/* npN(sector sec) 
 * 
 * 	Function returns the number of primitive states came from west.
 * 
 * 	Inputs:
 * 		sector sec - sector data structure
 * 	Outputs:
 * 		int npW - states of sec came from west
 * 
 * */
{
	return sec.npv[WEST];
}







//bin_sector *parent(bin_sector sect, int p)
/*	bin_sector *parent(sector sect, int p)
 * 
 * 	Function returns a pointer to the parent of sect whose state is indexed 
 * 	by p.
 * 
 * 	Inputs:
 * 		sector sect - sector data structure
 * 		int p - index of state in sect
 * 	Outputs:
 * 		sector *parent - pointer to parent sector 
 * 
 * */

// TO BE ADDED IN QSMATRIX
/*
{
  int q = sect.Q;
  int ds = sect.dS;
  int qParent = q+1, dsParent = ds+1;
  
    if(p < npS(sect)){
      //gend = SOUTH;
      qParent = q + 1;
      dsParent = ds;
    }
    else if(p < npS(sect)+npE(sect)){
      //gend = EAST;
      qParent = q;
      dsParent = ds - 1;
    }
    else if(p < npS(sect)+npE(sect)+npW(sect)){
      //gend = WEST;
      qParent = q;
      dsParent = ds+1;
    }
    else if (p < npS(sect)+npE(sect)+npW(sect)+npN(sect)){
      //gend = NORTH;
      qParent = q -1 ;
      dsParent = ds;
    }
    return &QS_matrix()[q_index(qParent)][dsParent];
}
*/



int gender(bin_sector sect, int p)
/*	int gender(sector sect, int p)
 * 
 * 	Function returns the gender (NORTH, EAST, SOUTH, WEST) of state p
 * 	in sector sect.
 * 
 * 	Inputs:
 * 		sector sect - sector data structure
 * 		int p - index of state in sect
 * 	Outputs:
 * 		int gender - direction which p came from
 * 
 * */
{
	
	int gend = SOUTH;
	
	if(p < npS(sect))
		gend = SOUTH;
    else if(p < npS(sect)+npE(sect))
		gend = EAST;
    else if(p < npS(sect)+npE(sect)+npW(sect))
		gend = WEST;
    else if (p < npS(sect)+npE(sect)+npW(sect)+npN(sect))
		gend = NORTH;
		
    return gend;
}





void print_sector(bin_sector sec, int ds, int curr, int print_translation_mat)
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
	if(ds > DSMAX)
		dSz = - dSz;
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
		cout << "state " << p << " = " << endl;
		for(int b = 0; b < sec.np_binconfs[p]; b++){
			int q_bin;
			int nsites = iter() + 2;
			if(curr == 1)
				nsites = iter()+1;
			int *fns_bin = new int[2*nsites];
			int index_bin = sec.bin_confs[accbin];
			int *binstate = index_to_binary(index_bin, nsites, fns_bin, q_bin, 0);
			cout << sec.t_coefs[accbin] << " " << print_conf(binstate, nsites) << endl;
			delete[] binstate;
			delete[] fns_bin;	
			accbin++; 
			}
			
		
		}
	
	


}

void clear_sector(bin_sector &sec, bool delete_vectors)
/*	void clear_sector(sector &sec, bool delete_vectors)
 * 
 *	Function deletes all fields in data structure sector sec.
 * 
 * 	Inputs:
 * 		sector &sec - sector address
 * 		bool delete_vectors - boolean flag to sinalize when deleting 
 * 			vectors are required
 * 	Outputs:
 * 		none 
 * 
 * */
{
	sec.Q =  0;
	sec.dS = 0;

	
	if(delete_vectors) { // if delete_vectors = false no deleting is required
                delete[] sec.np_binconfs;
                delete[] sec.bin_confs;
                delete[] sec.t_coefs;
	}
        (sec.np_binconfs) = NULL;
        (sec.bin_confs) = NULL;
        (sec.t_coefs) = NULL;


	sec.np = 0;
	sec.nb = 0;
	sec.nub = 0;

	for(unsigned int i = 0; i < (sec.npv).size(); i++)
		sec.npv[i] = 0;
}



vector<int> initial_index_in_bin_confs(bin_sector sec)
/*

	Function returns a container with the first index to be taken in bin_confs vector when the binary configurations of the p-th state is required

	Inputs:
		sec - bin sector data structure
	Outputs:
		beg_index_of_states - vector containing the starting index
*/
{
	int np_ = sec.np;
	vector<int> beg_index_of_state(np_);
 	beg_index_of_state[0] = 0;
	for(int p = 1; p < np_; p++){
		beg_index_of_state[p] = beg_index_of_state[p-1] + sec.np_binconfs[p-1];
	}
	return beg_index_of_state;
}


/*
double **translation_matrix(bin_sector sect){
  int np_states = sect.np;
  // for this, we must determine the number of unique configurations inside the sector
  set<int> unique_binconfs(sect.bin_confs, sect.bin_confs + sect.nb);
  int nb_states = unique_binconfs.size();
  double **T_pb = new_matrix(np_states, nb_states);
  for(int p = 0; p < np_states; p++){
    int nb_pth = sect.np_binconfs[p]; // number of binary configurations defining
    for(int b = 0; b < nb_states; b++){
      
      }
    }

}
* */


void print_sector_translation(bin_sector sec)
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
	vector<int> beg_states = initial_index_in_bin_confs(sec);
	cout << "number of states in sector: " << sec.np << endl;
	for(int p = 0; p < sec.np; p++){	
		cout << sec.np_binconfs[p] << " ";
	}
	cout << endl;
	for(int b = 0; b < sec.nb; b++){
		cout << sec.bin_confs[b] << " " ;
	}
	cout << endl;
	for(int b = 0; b < sec.nb; b++){
		cout << sec.t_coefs[b] << " " ;
		}
	cout << endl;
		
	
	beg_states.clear();
	
}
