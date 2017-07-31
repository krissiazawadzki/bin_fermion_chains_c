#include <cstdlib>
#include <iostream>

#include "int_array_util.h"


void bubble_sort_intarray(int *intArray, int dimArray)
/*
 * 
 * 	description:
 * 		function order an integer array using bubble sort and the input array is 
 * 		overwritten
 *
 * 	inputs:
 * 		- intArray: integer array
 * 		- dimArray: numer of elements
 * 
 * 	outputs:
 * 		- none
 * 
 * */
{
	for(int i = 0; i < dimArray; i++){
		for(int j = 0; j < dimArray-i-1; j++){
			if(intArray[j] > intArray[j+1]){
				int temp = intArray[j];
				intArray[j] = intArray[j+1];
				intArray[j+1] = temp;
				}
		}
	}
}


int find_unique_intarray(int *intArray, int dimArray, int *&unique_entries)
/*
 * 	description:
 * 		Function removes the duplicated entries of an integer array and store 
 * 		them in variable *&unique_entries, returning the number of unique entries
 * 		found
 * 
 * 	inputs:
 * 		- intArray: integer array
 * 		- dimArray: numer of elements
 * 		- unique_entries: array (not allocated previously) with the unique entries
 * 
 * 	outputs:
 * 		- dimUnique: number of unique entries found
 * */
{
	
	bubble_sort_intarray(intArray, dimArray);
	
	// counter starting from the first element
	int dimUnique = 1;
	int *aux_unique = new int[dimArray];
	aux_unique[0] = intArray[0];

	for(int i = 1; i < dimArray; i++){
		for(int j = dimUnique-1; j < dimUnique; j++){
			if(intArray[i] != aux_unique[j]){
					aux_unique[dimUnique] = intArray[i];
				dimUnique++;
				}
			}
		}
	std::cout << "number of unique things found " << dimUnique << std::endl; 
	// now, only copying unique entries	
	unique_entries = new int[dimUnique];
	for(int i = 0; i < dimUnique; i++){
		unique_entries[i] = aux_unique[i];
	}
		
	delete[] aux_unique;	
	return dimUnique;
}
