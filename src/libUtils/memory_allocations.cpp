// memory_allocaitons.cpp

using namespace std;

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>

#include "misc.h"
#include "memory_allocations.h"

//****************************************************************************
//malloc_1d_double_array()
//****************************************************************************
double *malloc_1d_array_double(long int Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	double *f = (double *) malloc(sizeof(double) * Nx);
	for (long int i=0; i<Nx; i++) f[i] = 0.0;
	
	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_double_array()");
	
	return f;
}

// ****************************************************************************
// realloc_1d_array_double()
// ****************************************************************************
double *realloc_1d_array_double(double *f_old, long int Nx)
{
	double *f;
	
	f = (double *) realloc(f_old,sizeof(double)*Nx);
	if(f==NULL)
	{
		fprintf(stderr,"In realloc_1d_array_double(%ld)...\n", Nx);
		quit_error((char*)"realloc failed");
	}
	
	return f;
}

/*********************************************************************************
malloc_1d_array_char()
*********************************************************************************/
char *malloc_1d_array_char(long int Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	char *f = (char *) malloc(sizeof(char) * Nx);
	
	if(!f) quit_error((char*)"Malloc failed in malloc_1D_char_array()");
	
	return f;
}

/*********************************************************************************
malloc_2d_array_char()
*********************************************************************************/
char **malloc_2d_array_char(long int Nx, long int Ny)
{
	if ( (Nx <= 0) || (Ny <= 0) ) quit_error((char*)"Can't allocate array of length <= 0");
	char **f;
	
	f = (char **) malloc(sizeof(char *) * Nx);
	if(!f) quit_error((char*)"Malloc failed in malloc_2D_array_char()");
	
	for(long int i=0 ; i<Nx ; i++) f[i] = malloc_1d_array_char(Ny);
	
	return f;
}

//****************************************************************************
//malloc_1d_int_array_fortran_indices()
//****************************************************************************
int *malloc_1d_array_int_fortran_indices(long int i0, long int i1)
{
	long int Nx = i1-i0+1;
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	int *f = (int *) malloc(sizeof(int) * Nx);
	
	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_int_array_fortran_indices()");
	
	return f - i0;
}

//****************************************************************************
//malloc_2d_int_array_fortran_indices()
//****************************************************************************
int **malloc_2d_array_int_fortran_indices(long int i0, long int i1, long int j0, long int j1)
{
	long int Nx = i1-i0+1;
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	int **f = (int **) malloc(sizeof(int*)*Nx);
	
	if( ! f ) quit_error((char*)"Malloc failed in malloc_2d_int_array_fortran_indices()");
	for(int i=0;i<Nx;i++)
		f[i]=malloc_1d_array_int_fortran_indices(j0, j1);
	
	return f - i0;
}

//****************************************************************************
//malloc_3d_int_array()
//****************************************************************************
int ***malloc_3d_array_int_fortran_indices(long int i0, long int i1, long int j0, long int j1, long int k0, long int k1)
{
	long int Nx = i1-i0+1;
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	int ***f = (int ***) malloc(sizeof(int**)*Nx);
	
	if(!f) quit_error((char*)"Malloc failed in malloc_3d_int_array_fortran_indices()");
	for(long int i=0 ; i<Nx ; i++)
	{
		f[i]=malloc_2d_array_int_fortran_indices(j0,j1,k0,k1);
	}
	
	return f - i0;
}

//****************************************************************************
//malloc_4d_int_array_fortran_indices()
//****************************************************************************
int ****malloc_4d_array_int_fortran_indices(long int I0, long int I1, long int i0, long int i1, long int j0, long int j1, long int k0, long int k1)
{
	long int Nx = I1-I0+1;
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	int ****f = (int ****) malloc(sizeof(int***)*Nx);
	
	if(!f) quit_error((char*)"Malloc failed in malloc_3d_int_array_fortran_vector_indices()");
	for(long int i=0 ; i<Nx ; i++)
	{
		f[i]=malloc_3d_array_int_fortran_indices(i0,i1,j0,j1,k0,k1);
	}
	
	return f - I0;
}

//****************************************************************************
//malloc_1d_double_array_fortran_indices()
//****************************************************************************
double *malloc_1d_array_double_fortran_indices(long int i0, long int i1)
{
	long int Nx = i1-i0+1;
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	double *f = (double *) malloc(sizeof(double) * Nx);
	
	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_double_array_fortran_indices()");
	
	return f - i0;
}

//****************************************************************************
//malloc_2d_double_array_fortran_indices()
//****************************************************************************
double **malloc_2d_array_double_fortran_indices(long int i0, long int i1, long int j0, long int j1)
{
	long int Nx = i1-i0+1;
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	double **f = (double **) malloc(sizeof(double*)*Nx);
	
	if( ! f ) quit_error((char*)"Malloc failed in malloc_2d_double_array_fortran_indices()");
	for(int i=0;i<Nx;i++)
		f[i]=malloc_1d_array_double_fortran_indices(j0, j1);
	
	return f - i0;
}

//****************************************************************************
//malloc_3d_double_array_fortran_indices()
//****************************************************************************
double ***malloc_3d_array_double_fortran_indices(long int i0, long int i1, long int j0, long int j1, long int k0, long int k1)
{
	long int Nx = i1-i0+1;
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	double ***f = (double ***) malloc(sizeof(double**)*Nx);
	
	if(!f) quit_error((char*)"Malloc failed in malloc_3d_double_array_fortran_indices()");
	for(long int i=0 ; i<Nx ; i++)
	{
		f[i]=malloc_2d_array_double_fortran_indices(j0,j1,k0,k1);
	}
	
	return f - i0;
}

//****************************************************************************
//malloc_3d_double_vector_array_fortran_indices()
//****************************************************************************
double ****malloc_4d_array_double_fortran_indices(long int I0, long int I1, long int i0, long int i1, long int j0, long int j1, long int k0, long int k1)
{
	long int Nx = I1-I0+1;
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	double ****f = (double ****) malloc(sizeof(double***)*Nx);
	
	if(!f) quit_error((char*)"Malloc failed in malloc_3d_double_array_fortran_indices()");
	for(long int i=0 ; i<Nx ; i++)
	{
		f[i]=malloc_3d_array_double_fortran_indices(i0,i1,j0,j1,k0,k1);
	}
	
	return f - I0;
}

//****************************************************************************
//malloc_1d_array_int()
//****************************************************************************
int *malloc_1d_array_int(long int Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	int *f = (int *) malloc(sizeof(int) * Nx);
	for (long int i=0; i<Nx; i++) f[i] = 0;
	
	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_int_array()");
	
	return f;
}

//****************************************************************************
//malloc_1d_array_long_int()
//****************************************************************************
long int *malloc_1d_array_long_int(long int Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	long int *f = (long int *) malloc(sizeof(long int) * Nx);
	for (long int i=0; i<Nx; i++) f[i] = 0;
	
	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_long_int_array()");
	
	return f;
}

//****************************************************************************
//malloc_1d_array_bool()
//****************************************************************************
bool *malloc_1d_array_bool(long int Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	bool *f = (bool *) malloc(sizeof(bool) * Nx);
	for (long int i=0; i<Nx; i++) f[i] = false;
	
	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_long_bool()");
	
	return f;
}

//****************************************************************************
//malloc_2d_double_array()
//****************************************************************************
double **malloc_2d_array_double(long int Nx, long int Ny)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	double **f = (double **) malloc(sizeof(double*) * Nx);
	
	if( ! f ) quit_error((char*)"Malloc failed in malloc_1d_double_array()");
	for(long int i=0 ; i<Nx ; i++)
		f[i] = malloc_1d_array_double(Ny);

	for (long int i=0; i<Nx; i++) 
		for (long int j=0; j<Ny; j++) 
			f[i][j] = 0.0;
	
	return f;
}

//****************************************************************************
//malloc_2d_int_array()
//****************************************************************************
long int **malloc_2d_array_long_int(long int Nx, long int Ny)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	long int **f = (long int **) malloc(sizeof(long int*)*Nx);
	
	if( ! f ) quit_error((char*)"Malloc failed in malloc_2d_int_long_array()");
	for(long int i=0;i<Nx;i++)
		f[i]=malloc_1d_array_long_int(Ny);
	
	for (long int i=0; i<Nx; i++) 
		for (long int j=0; j<Ny; j++) 
			f[i][j] = 0;
	
	return f;
}

//****************************************************************************
//malloc_2d_int_array()
//****************************************************************************
int **malloc_2d_array_int(long int Nx, long int Ny)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	int **f = (int **) malloc(sizeof(int*)*Nx);
	
	if( ! f ) quit_error((char*)"Malloc failed in malloc_2d_int_array()");
	for(int i=0;i<Nx;i++)
		f[i]=malloc_1d_array_int(Ny);
	
	for (long int i=0; i<Nx; i++) 
		for (long int j=0; j<Ny; j++) 
			f[i][j] = 0;
	
	return f;
}

//****************************************************************************
//malloc_3d_int_array()
//****************************************************************************
int ***malloc_3d_array_int(long int Nx, long int Ny, long int Nz)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	int ***f = (int ***) malloc(sizeof(int**) * Nx);
	
	if(!f) quit_error((char*)"Malloc failed in malloc_3d_int_array()");
	for(long int i=0 ; i<Nx ; i++) f[i]=malloc_2d_array_int(Ny, Nz);
		
	return f;
}

//****************************************************************************
//malloc_3d_double_array()
//****************************************************************************
double ***malloc_3d_array_double(long int Nx, long int Ny, long int Nz)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	double ***f = (double ***) malloc(sizeof(double**) * Nx);
	
	if(!f) quit_error((char*)"Malloc failed in malloc_3d_double_array()");
	for(long int i=0 ; i<Nx ; i++) f[i] = malloc_2d_array_double(Ny, Nz);
		
	return f;
}

//****************************************************************************
//malloc_4d_double_vector_array()
//****************************************************************************
double ****malloc_4d_array_double(int num_dimensions, long int Nx, long int Ny, long int Nz)
{
	if(num_dimensions <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	double ****f = (double ****) malloc(sizeof(double***) * num_dimensions);
	
	if(!f) quit_error((char*)"Malloc failed in malloc_3d_double_vector_array()");
	for(long int i=0 ; i < num_dimensions ; i++) f[i] = malloc_3d_array_double(Nx, Ny, Nz);
	
	return f;
}

//****************************************************************************
//malloc_4d_int_array()
//****************************************************************************
int ****malloc_4d_array_int(int num_dimensions, long int Nx, long int Ny, long int Nz)
{
	if(num_dimensions <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	int ****f = (int ****) malloc(sizeof(int***)*num_dimensions);
	
	if(!f) quit_error((char*)"Malloc failed in malloc_3d_int_vector_array()");
	for(long int i=0;i<num_dimensions;i++) f[i]=malloc_3d_array_int(Nx, Ny, Nz);
	
	return f;
}

//****************************************************************************
//free_1d_double_array()
//****************************************************************************
void free_1d_array_double(double *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_double_array()");
	free(f);

	return;
}

//****************************************************************************
//free_2d_double_array()
//****************************************************************************
void free_2d_array_double(double **f, long int Nx)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_2d_double_array()");
	for(long int i=0;i<Nx;i++)
		free_1d_array_double(f[i]);

	free(f);
	
	return;
}

//****************************************************************************
//free_3d_double_array()
//****************************************************************************
void free_3d_array_double(double ***f, long int Nx, long int Ny)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_3d_double_array()");
	for(long int i=0;i<Nx;i++)
		free_2d_array_double(f[i], Ny);

	free(f);
	
	return;
}

//****************************************************************************
//free_4d_double_array()
//****************************************************************************
void free_4d_array_double(double ****f, int num_dimensions, long int Nx, long int Ny)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_3d_double_array()");
	for(long int i=0;i<num_dimensions;i++)
		free_3d_array_double(f[i], Nx, Ny);
	
	free(f);
	
	return;
}


//****************************************************************************
//free_1d_int_array()
//****************************************************************************
void free_1d_array_int(int *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_int_array()");
	free(f);
	
	return;
}

//****************************************************************************
//free_2d_array_long_int()
//****************************************************************************
void free_2d_array_long_int(long int **f, long int Nx)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_2d_int_array()");
	for(long int i=0;i<Nx;i++)
		free_1d_array_long_int(f[i]);
	
	free(f);
	
	return;
}

//****************************************************************************
//free_2d_int_array()
//****************************************************************************
void free_2d_array_int(int **f, long int Nx)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_2d_int_array()");
	for(long int i=0;i<Nx;i++)
		free_1d_array_int(f[i]);
	
	free(f);
	
	return;
}

//****************************************************************************
//free_3d_int_array()
//****************************************************************************
void free_3d_array_int(int ***f, long int Nx, long int Ny)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_3d_int_array()");
	for(long int i=0;i<Nx;i++)
		free_2d_array_int(f[i], Ny);
	
	free(f);
	
	return;
}

//****************************************************************************
//free_4d_int_array()
//****************************************************************************
void free_4d_array_int(int ****f, int num_dimensions, long int Nx, long int Ny)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_3d_int_array()");
	for(long int i=0;i<num_dimensions;i++)
		free_3d_array_int(f[i], Nx, Ny);
	
	free(f);
	
	return;
}

//****************************************************************************
//free_1d_double_array_fortran_indices()
//****************************************************************************
void free_1d_array_double_fortran_indices(double *f,long int i0, long int i1)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_double_array()");
	free(f+i0);
	
	return;
}

//****************************************************************************
//free_2d_double_array_fortran_indices()
//****************************************************************************
void free_2d_array_double_fortran_indices(double **f,long int i0, long int i1, long int j0, long int j1)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_2d_double_array()");
	for(long int i=i0;i<=i1;i++)
		free_1d_array_double_fortran_indices(f[i],j0,j1);
	
	free(f+i0);
	
	return;
}

//****************************************************************************
//free_3d_double_array_fortran_indices()
//****************************************************************************
void free_3d_array_double_fortran_indices(double ***f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_3d_double_array()");
	for(long int i=i0;i<=i1;i++)
		free_2d_array_double_fortran_indices(f[i],j0,j1,k0,k1);
	
	free(f+i0);
	
	return;
}

//****************************************************************************
//free_4d_double_array_fortran_indices()
//****************************************************************************
void free_4d_array_double_fortran_indices(double ****f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1, long int l0, long int l1)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_4d_double_array()");
	for(long int i=i0;i<=i1;i++)
		free_3d_array_double_fortran_indices(f[i],j0,j1,k0,k1,l0,l1);
	
	free(f+i0);
	
	return;
}

//****************************************************************************
//free_1d_int_array_fortran_indices()
//****************************************************************************
void free_1d_array_int_fortran_indices(int *f,long int i0, long int i1)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_int_array()");
	free(f+i0);
	
	return;
}

//****************************************************************************
//free_2d_int_array_fortran_indices()
//****************************************************************************
void free_2d_array_int_fortran_indices(int **f,long int i0, long int i1, long int j0, long int j1)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_2d_int_array()");
	for(long int i=i0;i<=i1;i++)
		free_1d_array_int_fortran_indices(f[i],j0,j1);
	
	free(f+i0);
	
	return;
}

//****************************************************************************
//free_3d_int_array_fortran_indices()
//****************************************************************************
void free_3d_array_int_fortran_indices(int ***f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_3d_int_array()");
	for(long int i=i0;i<=i1;i++)
		free_2d_array_int_fortran_indices(f[i],j0,j1,k0,k1);
	
	free(f+i0);
	
	return;
}

//****************************************************************************
//free_4d_int_array_fortran_indices()
//****************************************************************************
void free_4d_array_int_fortran_indices(int ****f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1, long int l0, long int l1)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_4d_int_array()");
	for(long int i=i0;i<=i1;i++)
		free_3d_array_int_fortran_indices(f[i],j0,j1,k0,k1,l0,l1);
	
	free(f+i0);
	
	return;
}

//****************************************************************************
// free_1d_array_long_int()
//****************************************************************************
void free_1d_array_long_int(long int *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_array_integer()");
	free(f);
	
	return;
}

//****************************************************************************
// free_1d_array_bool()
//****************************************************************************
extern void free_1d_array_bool(bool *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_array_bool()");
	free(f);
	
	return;
}

//****************************************************************************
//free_1d_array_char()
//****************************************************************************
void free_1d_array_char(char *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_array_char()");
	free(f);
	
	return;
}

//****************************************************************************
//free_2d_array_char()
//****************************************************************************
void free_2d_array_char(char **f, long int Nx)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_2d_array_char()");
	for(long int i=0;i<Nx;i++)
		free_1d_array_char(f[i]);
	
	free(f);
	
	return;
}



