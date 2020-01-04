#include "string_ops.h"
#include <iostream>
#include <sstream>
#include <iomanip>

bool is_substring_in_string(std::string substring, std::string full_string)
/*	
 * 	description:
 * 		Function check if substring is contained in full_string
 * 
 * 	inputs:
 * 		- substring: string we want to find into other
 * 		- full_string: string in which we want to find substring
 * 
 * 	outputs:
 * 		- bool: true or false
 * 
 * */
{
	return full_string.find(substring) != std::string::npos;
}


std::vector<std::string> split(std::string text, std::string delims)
/*
 * 	description:
 * 		Function splits a string with delimiters ans returns a vector
 * 		containing string elements separated by it as elements
 * 
 * 	inputs:
 * 		- text: string to be splitted
 * 		- delims: delimiters
 *	outputs:
 * 		- tokens: vector containing splitted string
 * 
 * */
{
    std::vector<std::string> tokens;
    std::size_t start = text.find_first_not_of(delims), end = 0;

    while((end = text.find_first_of(delims, start)) != std::string::npos)
    {
        tokens.push_back(text.substr(start, end - start));
        start = text.find_first_not_of(delims, end);
    }
    if(start != std::string::npos)
        tokens.push_back(text.substr(start));

    return tokens;
}



std::string double2string(double number, int fmt_precision)
/*	string double2string(float number)
 * 
 * 	Function converts double to a string and 
 * 	returns it in the new format.
 * 
 * 	Inputs:
 * 		number - number to be converted to string format
 * 	Outputs:
 * 		string buff - converted number
 * 
 * */
{
	std::stringstream buff;
	buff << std::fixed << std::setprecision(fmt_precision) << number;
	return buff.str();
}


std::string int2string(int number, int fmt_precision)
/*	string int2string(float number)
 * 
 * 	Function converts int to a string and 
 * 	returns it in the new format.
 * 
 * 	Inputs:
 * 		number - number to be converted to string format
 * 	Outputs:
 * 		string buff - converted number
 * 
 * */
{
	std::stringstream buff;
	buff << std::fixed << std::setprecision(fmt_precision) << number;
	return buff.str();
}


std::string double2string_scientific(double number, int fmt_precision)
/*	string double2string(float number)
 * 
 * 	Function converts double to a string and 
 * 	returns it in the new format in scientific notation.
 * 
 * 	Inputs:
 * 		number - number to be converted to string format
 * 	Outputs:
 * 		string buff - converted number
 * 
 * */
{
	std::stringstream buff;
	buff << std::scientific << std::setprecision(fmt_precision) << number;
	return buff.str();
}
