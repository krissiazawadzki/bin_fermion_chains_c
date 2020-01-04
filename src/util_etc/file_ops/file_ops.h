#ifndef __FILE_OPS__

#include <string>
#include <fstream>

double **loadtxt(const char * file_name, int &n, int &m);
void savetxt(char *file_name, int **matrix, int n, int m);

void print_linear_stored_matrix(double *matrix, int n, int m);
void print_matrix_to_file(const char *namefile, double **matrix, int n, int m, const char *mode="w");
void print_linear_stored_matrix_to_file(const char *namefile, double *matrix, int n, int m, const char *mode="w");


template <typename T> void print_linear_matrix_to_file_cpp(std::string, T*, int, int, std::string mode = "w", int wfmt = 5, int precisionfmt=5);

template <typename T> void print_matrix_to_file_cpp(std::string, T**, int, int, std::string mode = "w", int wfmt = 5, int precisionfmt=5);

template <typename T> void print_matrix_to_file_scientific_cpp(std::string, T**, int, int, std::string mode = "w", int wfmt = 5, int precisionfmt=5);


template <typename T> void print_array_to_file_cpp(std::string, T*, int, std::string mode = "w", int wfmt = 5, int precisionfmt=5);


#define __FILE_OPS__
#endif
