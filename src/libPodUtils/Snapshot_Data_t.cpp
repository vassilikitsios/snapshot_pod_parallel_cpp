
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>

#include "Utils.h"
#include "Snapshot_Data_t.h"


//-------------------------------------------------------------------------
// Non class functions
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
Snapshot_Data_t **malloc_2d_array_Snapshot_Data_t(long int Nx, long int Ny)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	Snapshot_Data_t **f = (Snapshot_Data_t **) malloc(sizeof(Snapshot_Data_t*) * Nx);
	
	if( ! f ) quit_error((char*)"Malloc failed in malloc_2d_array_Snapshot_Data_t()");
	for(int i=0;i<Nx;i++)
		f[i]=malloc_1d_array_Snapshot_Data_t(Ny);
	
	return f;
}

//-------------------------------------------------------------------------
Snapshot_Data_t *malloc_1d_array_Snapshot_Data_t(long int Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	Snapshot_Data_t *f = (Snapshot_Data_t *) malloc(sizeof(Snapshot_Data_t) * Nx);
	
	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_array_Snapshot_Data_t()");
	
	return f;
}

//-------------------------------------------------------------------------
void free_2d_array_Snapshot_Data_t(Snapshot_Data_t **f, long int Nx)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_2d_array_Snapshot_Data_t()");
	for(long int i=0;i<Nx;i++)
		free_1d_array_Snapshot_Data_t(f[i]);
	
	free(f);
	
	return;
}

//-------------------------------------------------------------------------
void free_1d_array_Snapshot_Data_t(Snapshot_Data_t *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_array_Snapshot_Data_t()");
	free(f);
	
	return;
}

//-------------------------------------------------------------------------
