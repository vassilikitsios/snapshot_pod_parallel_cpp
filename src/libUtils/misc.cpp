// misc.h

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>

using namespace std;

#include "misc.h"

//****************************************************************************
//print_warning()
//****************************************************************************
void print_warning(char *fmt, ...)
{
	va_list args;
	
	va_start(args,fmt); // must be run first
	fprintf(stdout, "WARNING: ");
	vfprintf(stdout, fmt, args); // special printf that takes in the variable argument list
	fprintf(stdout,"\n");
	va_end(args); // must be called at the end
	
	return;
}

//****************************************************************************
//quit_error()
//This takes in a variable length argument list in the same way that fprintf
//does.
//****************************************************************************
void quit_error(char *fmt, ...)
{
	FILE *fp[2] = {stderr,stdout};
	va_list args;
	int i;
	
	for(i=0;i<=1;i++)
	{
		va_start(args,fmt); // must be run first
		fprintf(fp[i], "TERMINATING PROGRAM: ");
		vfprintf(fp[i], fmt, args); // special printf that takes in the variable argument list
		fprintf(fp[i],"\n");
		va_end(args); // must be called at the end
	}
	
	exit(EXIT_FAILURE);
}

//****************************************************************************
//quit_error_mpi()
//This takes in a variable length argument list in the same way that fprintf
//does.
//****************************************************************************
void quit_error_mpi(int myrank, int nprocs, char *fmt, ...)
{
	FILE *fp[2] = {stderr,stdout};
	va_list args;
	int i;
	
	for(i=0;i<=1;i++)
	{
		va_start(args,fmt); // must be run first
		fprintf(fp[i], "TERMINATING PROGRAM: ");
		vfprintf(fp[i], fmt, args); // special printf that takes in the variable argument list
		fprintf(fp[i]," on processor %d of %d", myrank, nprocs);
		fprintf(fp[i],"\n");
		va_end(args); // must be called at the end
	}
	
	exit(EXIT_FAILURE);
}
