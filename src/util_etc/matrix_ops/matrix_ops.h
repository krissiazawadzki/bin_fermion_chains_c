#ifndef __MATRIX_OPS__

#include <string>
#include <vector>

// allocation and delete funcions for double matrices
double **new_matrix(int n, int m);
double **zero_matrix(int n, int m);
void delete_matrix(double **matrix);

// allocation and delete funcions for int matrices
int **new_matrix_int(int n, int m);
int **zero_matrix_int(int n, int m);
void delete_matrix_int(int **matrix);
void clean_matrix(double **matrix, int n, int m);

// linear spaced arrays - equivalent to python functions
std::vector<double> linspace_vector(double start, double stop, int numpoints, bool endpoint=false);
double *linspace_array(double start, double stop, int numpoints, bool endpoint=false, bool repeat_midpoint=false);


// indexes conversion
int index_2D_to_1D(int line, int column, int n, int m);
void index_1D_to_2D(int index, int n, int m, int *line, int *column);

// io funcs
template <typename T> void print_matrix_cpp(T**, int, int, int wfmt = 5, int precisionfmt=5);
template <typename T> void print_array_cpp(T*, int, int wfmt = 5, int precisionfmt=5);
template <typename T> void print_linear_matrix_cpp(T*, int, int, int wfmt = 5, int precisionfmt=5);

void print_matrix(double **matrix, int n, int m);
void print_matrix_int(int **matrix, int n, int m);

void print_array_int(int *arrayv, int len);
void print_array(double *arrayv, int len);
void print_line_array(double *arrayv, int len);


void copy_piece_matrix(double **matdest, double **matsource, int destn, int destnm, int sourcen, int sourcem, int linec, int colc);

int find_intarray(std::vector<int> vint, int element);

int find_doublearray(std::vector<double> vdouble, double element);


/* formating functions */
std::string num2str(double number);

#define __MATRIX_OPS__
#endif
