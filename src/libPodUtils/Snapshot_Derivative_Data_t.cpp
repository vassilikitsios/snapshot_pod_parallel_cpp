
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>

#include "Utils.h"
#include "Snapshot_Derivative_Data_t.h"


//-------------------------------------------------------------------------
void Snapshot_Derivative_Data_t::apply_constant_test_field() 
{
	// u velocity set then v selected to satisfy continuity
	u = 1.0;			v = 1.0;				w = 0.0;
		
	du_dx = 0.0;		du_dy = 0.0;			du_dz = 0.0;
	dv_dx = 0.0;		dv_dy = -du_dx;			dv_dz = 0.0;
	dw_dx = 0.0;		dw_dy = 0.0;			dw_dz = 0.0;
	
	d2u_dx2 = 0.0;		d2v_dy2 = 0.0;			d2w_dz2 = 0.0;
	
	p = u;
}

//-------------------------------------------------------------------------
void Snapshot_Derivative_Data_t::apply_linear_test_field(double x, double y) 
{
	// u velocity set then v selected to satisfy continuity
	u = x + y;			v = -x - y;				w = 0.0;
	
	du_dx = 1.0;		du_dy = 1.0;			du_dz = 0.0;
	dv_dx = -1.0;		dv_dy = -du_dx;			dv_dz = 0.0;
	dw_dx = 0.0;		dw_dy = 0.0;			dw_dz = 0.0;
	
	d2u_dx2 = 0.0;		d2v_dy2 = 0.0;			d2w_dz2 = 0.0;
	
	p = u;
}

//-------------------------------------------------------------------------
void Snapshot_Derivative_Data_t::apply_quadratic_test_field(double x, double y) 
{
	// u velocity set then v selected to satisfy continuity
	u = x*x+x*y;		v = -2.0*x*y-0.5*y*y;	w = 0.0;
	
	du_dx = 2.0*x+y;	du_dy = x;				du_dz = 0.0;
	dv_dx = -2.0*y;		dv_dy = -du_dx;			dv_dz = 0.0;
	dw_dx = 0.0;		dw_dy = 0.0;			dw_dz = 0.0;
	
	d2u_dx2 = 2.0;		d2v_dy2 = -1.0;			d2w_dz2 = 0.0;
	
	p = u;
}

//-------------------------------------------------------------------------
// Non class functions
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
Snapshot_Derivative_Data_t **malloc_2d_array_Snapshot_Derivative_Data_t(long int Nx, long int Ny)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	Snapshot_Derivative_Data_t **f = (Snapshot_Derivative_Data_t **) malloc(sizeof(Snapshot_Derivative_Data_t*)*Nx);
	
	if( ! f ) quit_error((char*)"Malloc failed in malloc_2d_array_Snapshot_Derivative_Data_t()");
	for(int i=0;i<Nx;i++)
		f[i]=malloc_1d_array_Snapshot_Derivative_Data_t(Ny);
	
	return f;
}

//-------------------------------------------------------------------------
Snapshot_Derivative_Data_t *malloc_1d_array_Snapshot_Derivative_Data_t(long int Nx)
{
	if(Nx <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	Snapshot_Derivative_Data_t *f = (Snapshot_Derivative_Data_t *) malloc(sizeof(Snapshot_Derivative_Data_t) * Nx);
	
	if( !f ) quit_error((char*)"Malloc failed in malloc_1d_array_Snapshot_Derivative_Data_t()");
	
	return f;
}

//-------------------------------------------------------------------------
void free_2d_array_Snapshot_Derivative_Data_t(Snapshot_Derivative_Data_t **f, long int Nx)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_2d_array_Snapshot_Derivative_Data_t()");
	for(long int i=0;i<Nx;i++)
		free_1d_array_Snapshot_Derivative_Data_t(f[i]);
	
	free(f);
	
	return;
}

//-------------------------------------------------------------------------
void free_1d_array_Snapshot_Derivative_Data_t(Snapshot_Derivative_Data_t *f)
{
	if(!f)
		quit_error((char*)"Memory already deallocated in free_1d_array_Snapshot_Derivative_Data_t()");
	free(f);
	
	return;
}

//-------------------------------------------------------------------------
