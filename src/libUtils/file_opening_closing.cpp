// file_opening_closing.cpp

using namespace std;

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>

#include "misc.h"
#include "file_opening_closing.h"

//****************************************************************************
//open_file()
//****************************************************************************
 FILE *open_file(char *filename, char *mode)
{
	 return open_file(filename,mode,true);
}

//****************************************************************************
//open_file()
//****************************************************************************
FILE *open_file(char *filename, char *mode, bool verbose)
{
	FILE *fp;
	if (verbose == true) printf("Opening the file \"%s\" using the mode \"%s\"...\n", filename, mode);
	fp = fopen(filename,mode);
	if (fp==NULL) quit_error((char*)"Unable to open file \"%s\" using mode \"%s\".", filename, mode);
	
	return fp;
}

//****************************************************************************
//close_file()
//****************************************************************************
void close_file(FILE *fp, char *filename)
{
	close_file(fp,filename,true);
	return;
}

//****************************************************************************
//close_file()
//****************************************************************************
void close_file(FILE *fp, char *filename, bool verbose)
{
	if (verbose==true) printf("Closing the file \"%s\"...\n", filename);
	if (fp==NULL)
	{
		print_warning((char *)"Trying to close a file that has a NULL file pointer.");
		return;
	}
	else
		fclose(fp);
}


