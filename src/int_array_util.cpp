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
	// now, only copying unique entries	
	unique_entries = new int[dimUnique];
	for(int i = 0; i < dimUnique; i++){
		unique_entries[i] = aux_unique[i];
	}
		
	delete[] aux_unique;	
	

	return dimUnique;
}


void map_unique_indexes(int *intArray, int dimArray, int *uniqueArray, int dimUnique, int *inv_Array_map)
/* 
 * 	description:
 * 
 * 	function creates a map that assign to each element of the original array an index with the id for the 
 * 	unique binary configuration
 * 
 * 	EX:
 * 
 *  original array :	9	7	1	8	3	7	9	0	7	1	
 * 	unique ordered array:	0	1	3	7	8	9	
 * 	ids unique ordered array:	0->0	1->1	3->2	7->3	8->4	9->5	
 * 	map: 5	3	1	4	2	3	5	0	3	1
 * 
 * 	inputs:
 * 		- intArray: integer array
 * 		- dimArray: numer of elements
 * 		- unique_entries: array found with routine find_unique_intarray
 * 		- dimUnique: number of elements in unique array
 * 		- inv_Array_map: integer array with dimArray elements allocated previously
 * 
 *	outputs:
 * 		- none
 * 
 * */
{
	
	for(int i = 0; i < dimArray; i++){
		for(int ui = 0; ui < dimUnique; ui++){
			if(intArray[i] == uniqueArray[ui]){
				inv_Array_map[i] = ui;
				}
		}	
	}
}
