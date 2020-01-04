#ifndef __STRING_OPS__

#include <string>
#include <vector>


bool is_substring_in_string(std::string substring, std::string full_string);
std::vector<std::string> split(std::string text, std::string delims);


std::string double2string(double number, int fmt_precision=6);

std::string int2string(int number, int fmt_precision=6);

std::string double2string_scientific(double number, int fmt_precision=6);

#define __STRING_OPS__
#endif
