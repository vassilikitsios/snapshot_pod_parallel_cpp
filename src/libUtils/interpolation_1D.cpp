#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>

using namespace std;

#include "misc.h"
#include "memory_allocations.h"
#include "matrix.h"
#include "interpolation.h"
#include "interpolation_1D.h"

//-------------------------------------------------------------------------
Interpolation_1D_t::Interpolation_1D_t(long int num_data_points, int order, int max_redundancy, double epsilon, double singularity_tolerance)
{
	this->num_data_points = num_data_points;
	this->order = order;
	this->max_redundancy = max_redundancy;
	this->epsilon = epsilon;
	this->singularity_tolerance = singularity_tolerance;
	p = order+1;
	m = max_redundancy + p;
	
	distance_to_points = malloc_Dist_t(num_data_points);
	D = malloc_1d_array_double(p*2 * m*3);
}

//-------------------------------------------------------------------------
Interpolation_1D_t::~Interpolation_1D_t()
{
	free_Dist_t(distance_to_points);
}

//-------------------------------------------------------------------------
void Interpolation_1D_t::form_interpolation_function(double *x_data, double x_interp_point)
{
	order_closest_points(x_data, x_interp_point);
	if (order>0) generate_stencil(x_interp_point);
}

//-------------------------------------------------------------------------
void Interpolation_1D_t::order_closest_points(double *x_data, double x_interp_point)
{
	for (long int i=0; i<num_data_points; i++) {
		distance_to_points[i].node_num = i;		
		distance_to_points[i].delta_x = x_data[i] - x_interp_point;
		distance_to_points[i].delta_mag = sqrt( SQ(distance_to_points[i].delta_x) );
	}
	qsort(distance_to_points, num_data_points, sizeof(Dist_t), (int (*)(const void*, const void*)) qsort_distance_compare_function);
}

//-------------------------------------------------------------------------
void Interpolation_1D_t::generate_stencil(double x_interp_point)
{	
	// Allocate memory
	double *B, *w, *wB, *wBT, *wBT_wB, *wBT_wB_inv, *wBT_w;	
	B = malloc_1d_array_double(m*p);			w = malloc_1d_array_double(m*m);
	wB = malloc_1d_array_double(m*p);			wBT = malloc_1d_array_double(p*m);	
	wBT_wB = malloc_1d_array_double(p*p);		wBT_wB_inv = malloc_1d_array_double(p*p);	
	wBT_w = malloc_1d_array_double(p*m);
	
	// Generate weights
	generate_B_matrix(B);
	generate_w_matrix(w);
	
	// Find W.B
	matrix_multiply(w, m, m, B, m, p, wB);
	
	// Find transpose(W.B)
	matrix_copy(wB, m, p, wBT);
	matrix_transpose(wBT, m, p);
	
	// Find transpose(W.B).(W.B)
	matrix_multiply(wBT, p, m, wB, m, p, wBT_wB);
	
	// Find inv(transpose(W.B).(W.B))
	matrix_invert(wBT_wB, wBT_wB_inv, p);
	
	// Find transpose(W.B).W
	matrix_multiply(wBT, p, m, w, m, m, wBT_w);
	
	// Find D = inv(transpose(W.B).(W.B)) . transpose(W.B).W
	matrix_multiply(wBT_wB_inv, p, p, wBT_w, p, m, D);
	
	//Deallocate memeory
	free_1d_array_double(B);
	free_1d_array_double(w);
	free_1d_array_double(wB);
	free_1d_array_double(wBT);
	free_1d_array_double(wBT_wB);
	free_1d_array_double(wBT_wB_inv);
	free_1d_array_double(wBT_w);
}

//-------------------------------------------------------------------------
void Interpolation_1D_t::generate_w_matrix(double *w)
{
	// Gaussian weighting function
	for (int i=0; i<m; i++) w[i*m+i] = exp(-SQ(distance_to_points[i].delta_mag/epsilon));
}

//-------------------------------------------------------------------------
void Interpolation_1D_t::generate_B_matrix(double *B)
{
	long int pos=0;
	for (int i=0; i<m; i++) {
		B[pos] = 1.0;
		pos++;
		for (int j=1; j<=order; j++) {
			B[pos] = pow(distance_to_points[i].delta_x, 1.0*j);
			pos++;
		}
	}
}

//-------------------------------------------------------------------------
void Interpolation_1D_t::generate_C_matrix(double *interpolation_data, double *C)
{
	long int node_num=0;
	for (int i=0; i<m; i++) {
		node_num = distance_to_points[i].node_num;
		C[i] = interpolation_data[node_num];
	}
}

//-------------------------------------------------------------------------
void Interpolation_1D_t::evaluate_interpolation_function(double *interpolation_data)
{		
	double *C, *P;
	C = malloc_1d_array_double(m);
	generate_C_matrix(interpolation_data, C);
	
	if (order>0) {		
		// Find P = D . C
		P = malloc_1d_array_double(p);
		matrix_multiply(D, p, m, C, m, 1, P);

		//Evaluate function and its derivatives
		f = P[0];	
		df_dx = P[1];
		d2f_dx2 = 2.0*P[2];
		
		free_1d_array_double(P);
	} else {
		f = C[0];
		df_dx = 0.0;
		d2f_dx2 = 0.0;
	}

	free_1d_array_double(C);
}

//-------------------------------------------------------------------------

