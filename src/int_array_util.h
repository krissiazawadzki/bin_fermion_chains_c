#ifndef __INT_ARRAY_UTIL__

void bubble_sort_intarray(int *intArray, int dimArray);

int find_unique_intarray(int *intArray, int dimArray, int *&unique_entries);

void map_unique_indexes(int *intArray, int dimArray, int *uniqueArray, int dimUnique, int *inv_Array_map);

#define __INT_ARRAY_UTIL__
#endif
