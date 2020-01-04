#include "file_ops.h"
#include "matrix_ops/matrix_ops.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

double **loadtxt(const char * file_name, int &n, int &m)
{
	std::ifstream arq;
    arq.open(file_name);
    arq >> n;
    arq >> m;

    double **matrix = zero_matrix(n,m);
    

    for( int i = 0; i < n; i++ ) {
    	for( int j = 0; j < m; j++ ) {
    		arq >> matrix[i][j]; 
    	}
    }
    return matrix;

}


void savetxt(const char *file_name, int **matrix, int n, int m){
    std::string out_file_name(file_name);
    
    size_t position = out_file_name.find(".");
    if(position == std::string::npos)
        position = out_file_name.length();
    out_file_name.replace(position,5,".txt");

    std::ofstream saida;
    saida.open(out_file_name.c_str());
    
    for(int i=0;i<n;i++){
        for(int j=0; j<m; j++){
            if (matrix[i][j] == n+1) matrix[i][j] = -1;
            saida << std::setw(4) << matrix[i][j] ;
        }
        saida << "\n";     
    }
}


void print_matrix_to_file(const char *namefile, double **matrix, int n, int m, const char *mode)
{
  
  FILE * (file);
  file = fopen(namefile,mode);
  
  
  for(int i = 0; i< n; i++){
    for(int j = 0; j < m; j++){
      fprintf(file, "%.4e \t", matrix[i][j]);
    }
    fprintf (file,"\n");
  }
  fclose (file);
}


void print_linear_stored_matrix_to_file(const char *namefile, double *matrix, int n, int m, const char *mode)
{
  
  FILE * (file);
  file = fopen(namefile,mode);
  
  
  for(int i = 0; i< n; i++){
    for(int j = 0; j < m; j++){
      fprintf(file, "%.4e \t", matrix[i*m + j]);
    }
    fprintf (file,"\n");
  }
  fclose (file);
}


template <typename T>
void print_linear_matrix_to_file_cpp(std::string filename, T *matrix, int n, int m, std::string mode, int wfmt, int precisionfmt){
	std::ofstream file;
	if(mode == "a"){
		file.open(filename.c_str(), std::ios_base::app);
	  }
	else{
		file.open(filename);
		}
		
		
	int wfmtout = wfmt;
	
	if(precisionfmt > wfmt){
		wfmtout = precisionfmt + 1;
		}
	
	for(int i = 0; i< n; i++){
		for(int j = 0; j < m; j++){
			file << std::left << std::setw(wfmtout) << std::setprecision(precisionfmt) << matrix[i*m+j] << "\t";
			}
		file << std::endl;
		}
	file.close();	
	
}
template void print_linear_matrix_to_file_cpp<double>(std::string, double*, int, int, std::string, int, int);
template void print_linear_matrix_to_file_cpp<int>(std::string, int*, int, int, std::string, int, int);



template <typename T>
void print_matrix_to_file_cpp(std::string filename, T **matrix, int n, int m, std::string mode, int wfmt, int precisionfmt)
{
	std::ofstream file;
	if(mode == "a"){
		file.open(filename.c_str(), std::ios_base::app);
	  }
	else{
		file.open(filename);
		}
		
		
	int wfmtout = wfmt;
	
	if(precisionfmt > wfmt){
		wfmtout = precisionfmt + 1;
		}
	
	for(int i = 0; i< n; i++){
		for(int j = 0; j < m; j++){
			file << std::left << std::setw(wfmtout) << std::setprecision(precisionfmt) << matrix[i][j] << "\t";
			}
		file << std::endl;
		}
	file.close();
}
template void print_matrix_to_file_cpp<double>(std::string, double**, int, int, std::string, int, int);
template void print_matrix_to_file_cpp<int>(std::string, int**, int, int, std::string, int, int);


template <typename T>
void print_matrix_to_file_scientific_cpp(std::string filename, T **matrix, int n, int m, std::string mode, int wfmt, int precisionfmt)
{
	std::ofstream file;
	if(mode == "a"){
		file.open(filename.c_str(), std::ios_base::app);
	  }
	else{
		file.open(filename);
		}
		
		
	int wfmtout = wfmt;
	
	if(precisionfmt > wfmt){
		wfmtout = precisionfmt + 1;
		}
	
	std::cout.precision(5);
	for(int i = 0; i< n; i++){
		for(int j = 0; j < m; j++){
			file << std::left << std::setw(wfmtout) << std::setprecision(precisionfmt) << std::scientific << matrix[i][j] << "\t";
			}
		file << std::endl;
		}
	file.close();
}
template void print_matrix_to_file_scientific_cpp<double>(std::string, double**, int, int, std::string, int, int);
template void print_matrix_to_file_scientific_cpp<int>(std::string, int**, int, int, std::string, int, int);



template <typename T>
void print_array_to_file_cpp(std::string filename, T *array, int n, std::string mode, int wfmt, int precisionfmt)
{
	std::ofstream file;
	if(mode == "a"){
		file.open(filename.c_str(), std::ios_base::app);
	  }
	else{
		file.open(filename);
		}

	int wfmtout = wfmt;
	
	if(precisionfmt > wfmt){
		wfmtout = precisionfmt + 1;
		}
	
	for(int i = 0; i< n; i++){
			file << std::left << std::setw(wfmtout) << std::setprecision(precisionfmt) << array[i] << std::endl;
		}
	file.close();
}

template void print_array_to_file_cpp<double>(std::string, double*, int, std::string, int, int);
template void print_array_to_file_cpp<int>(std::string, int*, int, std::string, int, int);
