// string.cpp
// ****************************************************************************
using namespace std;
// ****************************************************************************
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>

#include "misc.h"
#include "string.h"

//----------------------------------------------------------------------------
int strcmp_case_insensitive(char *s1, char *s2)
{
	// continue comparing while both strings aren't finished
	while(*s1 != '\0' && *s2 != '\0')
	{
		if(tolower(*s1) > tolower(*s2)) return 1;
		if(tolower(*s1) < tolower(*s2)) return -1;
		s1++;
		s2++;
	}
	if(tolower(*s1) > tolower(*s2)) return 1;
	if(tolower(*s1) < tolower(*s2)) return -1;
	// otherwise the strings must be the same...
	return 0;
}

//----------------------------------------------------------------------------
char *copy_string(char *str1)
{
	int len;
	char *str2;
	
	len = strlen(str1);
	str2 = (char *) malloc(sizeof(char)*(len+1));
	if(str2 == NULL) quit_error((char*)"Malloc failed in copy_string()");
	strcpy(str2,str1);
	return str2;
}

//----------------------------------------------------------------------------
char *concatenate_strings(char *str1, char *str2)
{
	int len1, len2;
	char *str3;
	
	len1 = strlen(str1);
	len2 = strlen(str2);
	
	str3 = (char *) malloc(sizeof(char)*(len1+len2+1));
	if(str3 == NULL) quit_error((char*)"Malloc failed in copy_string()");
	strcpy(str3,str1);
	strcpy(str3+len1,str2);
	
	return str3;
}

//----------------------------------------------------------------------------
string convert_int_to_string(int number)
{
	ostringstream oss;
	oss << number;					// Works just like cout
	return oss.str();			// Return the underlying string
}

//----------------------------------------------------------------------------
string convert_int_to_zero_padded_string(int number)
{
	string mode_number_string, pad_zeros;
	mode_number_string = convert_int_to_string(number);
	if (number>9999) {
		quit_error("Too many modes to output, increase digits in output file names");
	} else if ( (number<9999) && (number>=1000) ) {
		pad_zeros = "";
	} else if ( (number<1000) && (number>=100) ) {
		pad_zeros = "0";
	} else if ( (number<100) && (number>=10) ) {
		pad_zeros = "00";
	} else if (number<10) {
		pad_zeros = "000";
	}
	return pad_zeros + mode_number_string;
}

//----------------------------------------------------------------------------
double char_to_double(char *input)
{
	std::stringstream ss(input);
	double result = 0;	
	ss >> result;
	return result;
}

//----------------------------------------------------------------------------
