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

//============================================================================
// Non class functions
//============================================================================

//-------------------------------------------------------------------------
Dist_t* malloc_Dist_t(long int N)
{
	if(N <= 0) quit_error((char*)"Can't allocate array of length <= 0");
	Dist_t *dist = (Dist_t *) malloc(sizeof(Dist_t) * N);
	
	if( !dist ) quit_error((char*)"Malloc failed in malloc_Dist_t()");
	
	return dist;
}

//-------------------------------------------------------------------------
void free_Dist_t(Dist_t *dist)
{
	if(!dist)
		quit_error((char*)"Memory already deallocated in free_Dist_t()");
	free(dist);
	
	return;
}

//-------------------------------------------------------------------------
int qsort_distance_compare_function(Dist_t *dist1, Dist_t *dist2)
{
	if(dist1->delta_mag > dist2->delta_mag) return 1;
	if(dist1->delta_mag < dist2->delta_mag) return -1;
    return 0;
}


//============================================================================
// Interpolation_t Functions
//============================================================================

//-------------------------------------------------------------------------
Interpolation_t::Interpolation_t(long int num_data_points, int order, int max_redundancy, double epsilon, double singularity_tolerance, bool two_PI_periodic)
{
	this->num_data_points = num_data_points;
	this->order = order;
	this->max_redundancy = max_redundancy;
	this->epsilon = epsilon;
	this->singularity_tolerance = singularity_tolerance;
	this->two_PI_periodic = two_PI_periodic;
	p = (order+1)*(order+2)/2;
	if (max_redundancy==-1) {
		m = num_data_points;
		max_redundancy = m - p;
		cout << "calculated redundancy = " << max_redundancy << endl;
	} else {
		m = max_redundancy + p;
	}
	
	if (max_redundancy<0) quit_error("max_redundancy<0");
	
	distance_to_points = malloc_Dist_t(num_data_points);
	// D = malloc_1d_array_double(p * (m + p));		// This will suffice for normal interpolation, but for the continuity conservation a larger memory allocation is required so the larger one is adopted at all times
	D = malloc_1d_array_double(p*2 * m*3);
}

//-------------------------------------------------------------------------
Interpolation_t::~Interpolation_t()
{
	free_Dist_t(distance_to_points);
	free_1d_array_double(D);
}

//-------------------------------------------------------------------------
void Interpolation_t::form_interpolation_function_on_extruded_2D_data(double **x_data, double x_interp_point, double y_interp_point)
{
	order_closest_nodes_on_extruded_2D_data(x_data, x_interp_point, y_interp_point);
	if (order>0) generate_stencil(x_interp_point, y_interp_point);
}

//-------------------------------------------------------------------------
void Interpolation_t::order_closest_nodes_on_extruded_2D_data(double **x_data, double x_interp_point, double y_interp_point)
{
	double z_max = SMALL_NUMBER, z_min = BIG_NUMBER, z_tol;	
	for (int i=0; i<num_data_points; i++) {
		if(x_data[2][i] < z_min) z_min = x_data[2][i];
		if(x_data[2][i] > z_max) z_max = x_data[2][i];
	}
	z_tol = (z_max - z_min)/1.0e3;
	if (abs(z_tol) < ZERO) z_tol = ZERO;
	
	for (int i=0; i<num_data_points; i++) {
		distance_to_points[i].delta_x = BIG_NUMBER * 1.0;
		distance_to_points[i].delta_y = BIG_NUMBER * 1.0;
		distance_to_points[i].delta_mag = BIG_NUMBER * 1.0;
		distance_to_points[i].node_num = i;

		if (ABS(x_data[2][i]-z_max)<z_tol) {
			distance_to_points[i].delta_x = x_data[0][i] - x_interp_point;
			if(two_PI_periodic) { 
				if (abs(distance_to_points[i].delta_x) > abs(distance_to_points[i].delta_x - 2*PI)) {
					distance_to_points[i].delta_x = distance_to_points[i].delta_x - 2*PI;
				} else if (abs(distance_to_points[i].delta_x) > abs(distance_to_points[i].delta_x + 2*PI)) {
					distance_to_points[i].delta_x = distance_to_points[i].delta_x + 2*PI;
				}
			}	
			distance_to_points[i].delta_x = distance_to_points[i].delta_x;
			distance_to_points[i].delta_y = (x_data[1][i] - y_interp_point);
			distance_to_points[i].delta_mag = sqrt( SQ(distance_to_points[i].delta_x) + SQ(distance_to_points[i].delta_y) );
		}
		
	}
	qsort(distance_to_points, num_data_points, sizeof(Dist_t), (int (*)(const void*, const void*)) qsort_distance_compare_function);	
}

//-------------------------------------------------------------------------
void Interpolation_t::form_interpolation_function_on_2D_data(double **x_data, double x_interp_point, double y_interp_point)
{
	order_closest_nodes_on_2D_data(x_data, x_interp_point, y_interp_point);
	if (order>0) generate_stencil(x_interp_point, y_interp_point);
}

//-------------------------------------------------------------------------
void Interpolation_t::order_closest_nodes_on_2D_data(double **x_data, double x_interp_point, double y_interp_point)
{
	for (long int i=0; i<num_data_points; i++) {
		distance_to_points[i].delta_x = BIG_NUMBER * 1.0;
		distance_to_points[i].delta_y = BIG_NUMBER * 1.0;
		distance_to_points[i].delta_mag = BIG_NUMBER * 1.0;
		distance_to_points[i].node_num = i;
		
		distance_to_points[i].delta_x = x_data[0][i] - x_interp_point;
		if(two_PI_periodic) { 
			if (abs(distance_to_points[i].delta_x) > abs(distance_to_points[i].delta_x - 2*PI)) {
				distance_to_points[i].delta_x = distance_to_points[i].delta_x - 2*PI;
			} else if (abs(distance_to_points[i].delta_x) > abs(distance_to_points[i].delta_x + 2*PI)) {
				distance_to_points[i].delta_x = distance_to_points[i].delta_x + 2*PI;
			}
		}
		distance_to_points[i].delta_x = distance_to_points[i].delta_x;
		distance_to_points[i].delta_y = (x_data[1][i] - y_interp_point);
		distance_to_points[i].delta_mag = sqrt( SQ(distance_to_points[i].delta_x) + SQ(distance_to_points[i].delta_y) );
		
	}
	qsort(distance_to_points, num_data_points, sizeof(Dist_t), (int (*)(const void*, const void*)) qsort_distance_compare_function);
}

//-------------------------------------------------------------------------
void Interpolation_t::generate_stencil_and_test_stability(double x_interp_point, double y_interp_point)
{	
	// Determine best stencil
	double *B;
	long int num_singular_values_greater_than_tol = 1;

	if (max_redundancy==0) {
		m = p;
	} else {
		m = p-1;
		while ( (num_singular_values_greater_than_tol > 0) && (m<max_redundancy + p) && (distance_to_points[m].delta_mag < BIG_NUMBER * 1.0) ) {
			m++;
			B = malloc_1d_array_double(m*p);
			generate_B_matrix(B);
			num_singular_values_greater_than_tol=find_null_space_decomposition(B, m, p, singularity_tolerance);
			free_1d_array_double(B);
		}
		if (m==max_redundancy + p) {
			if (distance_to_points[m].delta_mag < BIG_NUMBER * 1.0) {
				cout << "WARNING: Singular stencil with maximum number of redundant points. num_singular_values_greater_than_tol = " << num_singular_values_greater_than_tol << endl;
			} else {
				cout << "WARNING: Singular stencil and all points in domain are used. num_singular_values_greater_than_tol = " << num_singular_values_greater_than_tol << endl;
			}
		}
		cout << "Stable stencil with, p, m, redundancy: " << p << " " << m << " " << m-p << endl;
	}	
	
	// Allocate memory
	double *w, *wB, *wBT, *wBT_wB, *wBT_wB_inv, *wBT_w;	
	B = malloc_1d_array_double(m*p);
	w = malloc_1d_array_double(m*m);
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
void Interpolation_t::generate_stencil(double x_interp_point, double y_interp_point)
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
void Interpolation_t::generate_w_matrix(double *w)
{
	// Gaussian weighting function
	for (int i=0; i<m; i++) w[i*m+i] = exp(-SQ(distance_to_points[i].delta_mag/epsilon));
}

//-------------------------------------------------------------------------
void Interpolation_t::generate_b_matrix(double *b)
{
	long int pos=0;
	for (int i=0; i<m; i++) {
		for (int j=0; j<=order; j++) {
			for (int k=0; k<=j; k++) {
				if (k==j) {
					if (k==0) {
						b[pos] = 1.0;
					} else {
						b[pos] = pow(distance_to_points[i].delta_y, 1.0*k);
					}
				} else {
					if (k==0) {
						b[pos] = pow(distance_to_points[i].delta_x, 1.0*(j-k));
					} else {
						b[pos] = pow(distance_to_points[i].delta_x, 1.0*(j-k)) * pow(distance_to_points[i].delta_y, 1.0*k);
					}
				}
				pos++;
			}
		}
	}
}

//-------------------------------------------------------------------------
void Interpolation_t::generate_b_x_matrix(double *b_x)
{
	long int pos=0;
	for (int i=0; i<m; i++) {
		for (int j=0; j<=order; j++) {
			for (int k=0; k<=j; k++) {
				if (k==j) {
					b_x[pos] = 0.0;
				} else {
					if (k==0) {
						b_x[pos] = 1.0 * j * pow(distance_to_points[i].delta_x, 1.0*(j-1));
					} else {
						b_x[pos] = 1.0 * (j-k) * pow(distance_to_points[i].delta_x, 1.0*(j-k-1)) * pow(distance_to_points[i].delta_y,k);
					}
				}
				pos++;
			}
		}
	}
}

//-------------------------------------------------------------------------
void Interpolation_t::generate_b_y_matrix(double *b_y)
{
	long int pos=0;
	for (int i=0; i<m; i++) {
		for (int j=0; j<=order; j++) {
			for (int k=0; k<=j; k++) {
				if (k==0) {
					b_y[pos] = 0.0;
				} else if (k==1) {
					b_y[pos] = pow(distance_to_points[i].delta_x, 1.0*j);
				} else if (k==j-1) {
					b_y[pos] = 1.0 * k * pow(distance_to_points[i].delta_y, k-1);
				} else {
					b_y[pos] = 1.0 * k * pow(distance_to_points[i].delta_x, 1.0*(j-k-1)) * pow(distance_to_points[i].delta_y, k-1);
				}
				pos++;
			}
		}
	}
}

//-------------------------------------------------------------------------
void Interpolation_t::generate_B_matrix(double *B)
{
	generate_b_matrix(B);
}

//-------------------------------------------------------------------------
void Interpolation_t::evaluate_interpolation_function(double *interpolation_data)
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
		df_dy = P[2];
		d2f_dx2 = 2.0*P[3];
		d2f_dxdy = P[4];
		d2f_dy2 = 2.0*P[5];
		
		free_1d_array_double(P);
	} else {
		f = C[0];
		df_dx = 0.0;
		df_dy = 0.0;
		d2f_dx2 = 0.0;
		d2f_dxdy = 0.0;
		d2f_dy2 = 0.0;
	}

	free_1d_array_double(C);
}

//-------------------------------------------------------------------------
void Interpolation_t::generate_C_matrix(double *interpolation_data, double *C)
{
	long int node_num=0;
	for (int i=0; i<m; i++) {
		node_num = distance_to_points[i].node_num;
		C[i] = interpolation_data[node_num];
	}
}

//-------------------------------------------------------------------------

