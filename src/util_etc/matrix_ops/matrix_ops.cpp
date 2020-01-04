#include "matrix_ops.h"
#include <iostream>
#include <cfloat>
#include <cstring>


#include <fstream>
#include <cstdlib>
#include <sstream>

#include <algorithm>
#include <iomanip>

using namespace std;


double **new_matrix(int n, int m) 
/*	double **new_matrix(int n, int m) 
 * 
 * 	Function allocates space for a bidimensional array (matrix). Each 
 * 	line has a pointer.
 * 
 * 	Inputs:
 * 		int n - line dimension
 * 		int m - column dimension
 * 
 * 	Outputs:
 * 		double ** matrix - allocated matrix	
 * 
 * */
{
    /* Aloca um vetor de ponteiros para as linhas da matriz */
    double **matrix = new double*[n];
    /* Aloca uma região contígua de memória para armazenar toda a matriz */
    matrix[0] = new double[n*m];
    /* Armazena um ponteiro para o início de cada linha */
    for(int i = 1; i < n; i++)
        matrix[i] = matrix[i-1] + m;
    return matrix;
}

double **zero_matrix(int n, int m) 
/*	void zero_matrix(int n, int m) 
 * 
 * 	Function initialize a matrix with all elements 0.
 * 
 * 	Inputs:
 * 		double **matrix - matrix to be initialized
 * 		int n - line dimension
 * 		int m - column dimension
 * 
 * 	Outputs:
 * 		none
 * */
{
    double** matrix = new_matrix(n,m);
    memset(matrix[0], 0, n*m*sizeof(double));   // assume IEEE 754
    return matrix;	
}

void delete_matrix(double **matrix) 
/*	void delete_matrix(double **matrix) 
 * 
 * 	Function deletes matrix and deallocates all lines.
 * 
 * 	Inputs:
 * 		double **matrix - matrix to be deallocated
 * 
 * 	Outputs:
 * 		none
 * 
 * */
{
    delete[] matrix[0];
    delete[] matrix;
}


int **new_matrix_int(int n, int m) 
/*	int **new_matrix(int n, int m) 
 * 
 * 	Function allocates space for a bidimensional array (matrix). Each 
 * 	line has a pointer.
 * 
 * 	Inputs:
 * 		int n - line dimension
 * 		int m - column dimension
 * 
 * 	Outputs:
 * 		int ** matrix - allocated matrix	
 * 
 * */
{
    /* Aloca um vetor de ponteiros para as linhas da matriz */
    int **matrix = new int*[n];
    /* Aloca uma região contígua de memória para armazenar toda a matriz */
    matrix[0] = new int[n*m];
    /* Armazena um ponteiro para o início de cada linha */
    for(int i = 1; i < n; i++)
        matrix[i] = matrix[i-1] + m;
    return matrix;
}


int** zero_matrix_int(int n, int m) 
/*	void zero_matrix_int(int **matrix, int n, int m) 
 * 
 * 	Function initialize a matrix with all elements 0.
 * 
 * 	Inputs:
 * 		int **matrix - matrix to be initialized
 * 		int n - line dimension
 * 		int m - column dimension
 * 
 * 	Outputs:
 * 		none
 * */
{
    int** matrix = new_matrix_int(n,m);
    memset(matrix[0], 0, n*m*sizeof(int));   // assume IEEE 754
    return matrix;
}

void delete_matrix_int(int **matrix) 
/*	void delete_matrix(int **matrix) 
 * 
 * 	Function deletes matrix and deallocates all lines.
 * 
 * 	Inputs:
 * 		int **matrix - matrix to be deallocated
 * 
 * 	Outputs:
 * 		none
 * 
 * */
{
    delete[] matrix[0];
    delete[] matrix;
}


void clean_matrix(double **matrix, int n, int m)
/*
	void clean_matrix(double **matrix, int n, int m)

	Function sets all elements of pre-existing matrix with 
	value 0.

	Inputs:
		double **matrix - pointer to matrix
		int n - matrix number of lines
		int m - matrix number of columns 

	Outputs:	
		none

*/
{
	memset(matrix[0], 0, n*m*sizeof(double)); 
}



std::vector<double> linspace_vector(double start, double stop, int numpoints, bool endpoint)
/*
 * 	Function generates an array containing values between start and stop linearly spaced.
 * 
 * 	Inputs:
 * 		double start - starting value of interval
 * 		double stop - final value
 * 		int numpoints - number of points to divide the interval [start, stop]
 * 		bool endpoint - if true, stop value is included [start, stop]
 * 						if false, stop value is supressed [start, stop)
 * 	
 * 	Outputs:
 * 		vector<double> vector containing the result;
 * 
 * */
{

	double step = (stop - start);
	if(endpoint == true){
		step /= numpoints -1;
		}
	else{
		step /= numpoints;
		}
		
	vector<double> lsvector(numpoints);
	lsvector[0] = start;
	for(int i = 1; i < numpoints; i++){
		lsvector[i] = lsvector[i-1] + step;
		}
	return lsvector;
}

double *linspace_array(double start, double stop, int numpoints, bool endpoint, bool repeat_midpoint)
/*
 * 	Function generates an array containing values between start and stop linearly spaced.
 * 
 * 	Inputs:
 * 		double start - starting value of interval
 * 		double stop - final value
 * 		int numpoints - number of points to divide the interval [start, stop]
 * 		bool endpoint - if true, stop value is included [start, stop]
 * 						if false, stop value is supressed [start, stop)
 * 		bool repeat_midpoint - if true, mid point is twice included in vector
 * 						default is false
 * 	
 * 	Outputs:
 * 		double *lsarray array containing the result;
 * 
 * 	NOTE: it must be deleted after usage
 * 
 * */
{
	if(repeat_midpoint){
		double midpoint = (stop - start) * 0.5;
		int halfpoints = numpoints / 2;
		double step = (midpoint - start) / (halfpoints-1);
		double *lsarray = new double [numpoints];
		lsarray[0] = start;
		for(int i = 1; i < halfpoints; i++){
			lsarray[i] = lsarray[i-1] + step;
			} 
		lsarray[halfpoints] = lsarray[halfpoints-1];	
		for(int i = halfpoints + 1; i < numpoints; i++){
			lsarray[i] = lsarray[i-1] + step;
			} 		
		return lsarray;
	}

	double step = (stop - start);
	if(endpoint == true){
		step /= numpoints -1;
		}
	else{
		step /= numpoints;
		}
		
	double *lsarray = new double [numpoints];
	lsarray[0] = start;
	for(int i = 1; i < numpoints; i++){
		lsarray[i] = lsarray[i-1] + step;
		}
	return lsarray;
}


/*
 * 
 * 
 * 		indexes conversion
 * 
 * 
 * */
 
int index_2D_to_1D(int line, int column, int n, int m){
	int index = line*m + column;
	return index;
}

void index_1D_to_2D(int index, int n, int m, int *line, int *column)
{
	*line = index / m ; // when call pass &line and &column
	*column = index % m;
}


/*	I/O functions to print a matrix
 *
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * */ 
template <typename T>
void print_matrix_cpp(T **matrix, int n, int m, int wfmt, int precisionfmt){
	
	int wfmtout = wfmt;
	
	if(precisionfmt > wfmt){
		wfmtout = precisionfmt + 1;
		}
	
	for(int i = 0; i< n; i++){
		for(int j = 0; j < m; j++){
			std::cout << std::left << std::setw(wfmtout) << std::setprecision(precisionfmt) << matrix[i][j] << "\t";
			}
		std::cout << std::endl;
		}	
}
template void print_matrix_cpp<double>(double**, int, int, int, int);
template void print_matrix_cpp<int>(int**, int, int, int, int);




template <typename T>
void print_linear_matrix_cpp(T *matrix, int n, int m, int wfmt, int precisionfmt){
	
	int wfmtout = wfmt;
	
	if(precisionfmt > wfmt){
		wfmtout = precisionfmt + 1;
		}
	
	for(int i = 0; i< n; i++){
		for(int j = 0; j < m; j++){
			std::cout << std::left << std::setw(wfmtout) << std::setprecision(precisionfmt) << matrix[i*m+j] << "\t";
			}
		std::cout << std::endl;
		}	
}
template void print_linear_matrix_cpp<double>(double*, int, int, int, int);
template void print_linear_matrix_cpp<int>(int*, int, int, int, int);




template <typename T>
void print_array_cpp(T *array, int n, int wfmt, int precisionfmt){
	
	int wfmtout = wfmt;
	
	if(precisionfmt > wfmt){
		wfmtout = precisionfmt + 1;
		}
	
	for(int i = 0; i< n; i++){
		std::cout << std::left << std::setw(wfmtout) << std::setprecision(precisionfmt) << array[i] << "\t";
	}	
	std::cout << std::endl;
}
template void print_array_cpp<double>(double*, int,  int, int);
template void print_array_cpp<int>(int*, int, int, int);




void print_matrix(double **matrix, int n, int m){
  for(int i=0;i<n;i++){
    for(int j=0; j<m; j++){
      printf("%.7f \t ", matrix[i][j]);
    }
    cout << "\n";     
  }
}

void print_matrix_int(int **matrix, int n, int m){
  for(int i=0;i<n;i++){
    for(int j=0; j<m; j++){
      printf("%d \t ", matrix[i][j]);
    }
    cout << "\n";     
  }
}


void print_array_int(int *arrayv, int len){
		for(int i = 0; i < len; i++)
			printf("%d \t ",arrayv[i]);
	cout << endl;
}


void print_array(double *arrayv, int len){
		for(int i = 0; i < len; i++)
  printf("%.4e \t ",arrayv[i]);
  cout << endl;
}


void print_line_array(double *arrayv, int len){
		for(int i = 0; i < len; i++)
  printf("%.4e \n ",arrayv[i]);
  cout << endl;
}





// IO FRUFRUS
string num2str(double number){
	ostringstream buff;
	buff<<number;
	return buff.str();
}

void print_linear_stored_matrix(double *matrix, int n, int m)
{
  for(int i=0;i<n;i++){
    for(int j=0; j<m; j++){
      printf("%.4e \t ", matrix[i*m+j]);
    }
    cout << "\n";
  }
}



// copy

void copy_piece_matrix(double **matdest, double **matsource, int destn, int destnm, int sourcen, int sourcem, int linec, int colc){
	for(int i = 0; i < sourcen; i++)
		for(int j = 0; j < sourcem; j++)
			matdest[i+linec][j+colc] = matsource[i][j];
}







// VECTOR ROUTINES
int find_intarray(std::vector<int> vint, int element){
  int founindex = -1;
  for(unsigned int i = 0; i < vint.size(); i++){
    if(vint[i] == element){
      founindex = i;
    }
  }

  return founindex;

}


int find_doublearray(std::vector<double> vdouble, double element){
  int founindex = -1;
  for(unsigned int i = 0; i < vdouble.size(); i++){
    if(vdouble[i] == element){
      founindex = i;
    }
  }
  
  return founindex;
  
}


