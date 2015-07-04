#include <cmath>
#include <iostream>

#include "memory_allocations.h"
#include "misc.h"
#include "cubic_spline.h"

//----------------------------------------------------------------------------
Cubic_Spline_t::Cubic_Spline_t(long int num_points)
{
	this->num_points = num_points;
	a = malloc_1d_array_double(num_points+1);
	c = malloc_1d_array_double(num_points+1);
	b = malloc_1d_array_double(num_points);
	d = malloc_1d_array_double(num_points);
	points = malloc_1d_array_double(num_points);
}

//----------------------------------------------------------------------------
Cubic_Spline_t::~Cubic_Spline_t()
{
	free_1d_array_double(a);
	free_1d_array_double(c);
	free_1d_array_double(b);
	free_1d_array_double(d);
	free_1d_array_double(points);
}

//----------------------------------------------------------------------------
void Cubic_Spline_t::calculate_interpolation_function(double *points_in, double *data)
{
	// Need to solve the system Ac=Q for 'c' to then determine 'b' and 'd' coefficients	
	
	// Allocate memory
	double *h, *diag1, *diag2, *diag3, *Q;
	h = malloc_1d_array_double(num_points);
	diag1 = malloc_1d_array_double(num_points);
	diag2 = malloc_1d_array_double(num_points);	
	diag3 = malloc_1d_array_double(num_points);
	Q = malloc_1d_array_double(num_points);
	
	// Copt across the points
	for (long int i=0; i<num_points; i++) points[i] = points_in[i];
	
	// Calculate spacinings between spline points
	for (long int i=0; i<num_points-1; i++)	{
		h[i] = points[i+1] - points[i];
//		fprintf(stdout, "h = %f\n",h[i]);
	}
	
	// Assign 'a' coefficients
	for (long int i=0; i<num_points; i++) {
		a[i] = data[i];
//		fprintf(stdout, "a = %f\n",a[i]);
	}
	
	// Calculate diagonal elements of tri-diagonal 'A' matrix
	for (long int i=0; i<num_points; i++) {
		diag1[i] = h[i];
//		fprintf(stdout, "diag1 = %f\n",diag1[i]);
	}
	diag1[num_points-2] = 0.0;
	
	for (long int i=0; i<num_points; i++) {
		diag2[i] = 2.0 * ( h[i-1] + h[i] );
//		fprintf(stdout, "diag2 = %f\n",diag2[i]);		
	}
	diag2[0] = 1.0;
	diag2[num_points-1] = 1.0;

	for (long int i=0; i<num_points; i++) {
		diag3[i] = h[i];
//		fprintf(stdout, "diag3 = %f\n",diag3[i]);
	}
	diag3[0] = 0.0;
	
	
	// Calculate elements of 'Q' matrix
	for (long int i=1; i<num_points-1; i++) {
		Q[i] = 3.0 / h[i] * (a[i+1]-a[i]) + 3.0 / h[i-1] * (a[i-1]-a[i]);
//		fprintf(stdout, "Q = %f\n",Q[i]);
	}
	Q[0] = 0.0;
	Q[num_points-1] = 0.0;
	
	// Solve for 'c' coefficients using Gaussian Elimination and then Backward substitution
	for (long int i=1; i<num_points-1; i++) {
//		fprintf(stdout, "%f %f %f %f %f\n", Q[i], diag2[i], diag3[i-1], diag1[i-1], diag2[i-1]);
		diag2[i] = diag2[i] - diag3[i-1] * diag1[i-1] / diag2[i-1];
		Q[i] = Q[i] - Q[i-1] * diag1[i-1] / diag2[i-1];
//		fprintf(stdout, "Q = %f\n", Q[i]);
	}
	
	c[num_points-1] = Q[num_points-1] / diag2[num_points-1];
	for(long int i=num_points-1 ; i>=0 ; i--) {
		c[i] = (Q[i] - diag3[i] * c[i+1]) / diag2[i];
//		fprintf(stdout, "c = %f\n",c[i]);
	}
	
	// Back out 'b' and 'd' coefficients
	for (int i=0; i<num_points-1; i++) {
		b[i] = ( a[i+1] - a[i] ) / h[i] - h[i] / 3.0 * (2 * c[i] + c[i+1]);
		d[i] = ( c[i+1] - c[i] ) / 3.0 / h[i];
//		fprintf(stdout, "b = %f d = %f\n",b[i], d[i]);
	}	
	
	// Deallocate memory
	free_1d_array_double(h);
	free_1d_array_double(diag1);
	free_1d_array_double(diag2);
	free_1d_array_double(diag3);
	free_1d_array_double(Q);
}

//----------------------------------------------------------------------------
double Cubic_Spline_t::evaluate_function(double interp_point)
{
	double f=0;
	for (long int i=0; i<num_points; i++) {
		if ( ( interp_point >= points[i] ) && ( interp_point <= points[i+1] ) ) {
			f = a[i] + b[i] * ( interp_point - points[i] )
						+ c[i] * SQ( interp_point - points[i])
						+ d[i] * CU( interp_point - points[i]);
//			fprintf(stdout, "within class %f %f %f %f\n", f, interp_point, points[i], points[i+1]);
		}
	}
	return f;
}

//----------------------------------------------------------------------------
double Cubic_Spline_t::evaluate_first_derivative(double interp_point)
{
	double df=0;
	for (long int i=0; i<num_points; i++) {
		if ( ( interp_point >= points[i] ) && ( interp_point <= points[i+1] ) ) {
			df = b[i] + 2.0 * c[i] * ( interp_point - points[i] )
					+ 3.0 * d[i] * SQ( interp_point - points[i] );
//			fprintf(stdout, "%f\n",df);
		}
	}
	return df;
}

//----------------------------------------------------------------------------
double Cubic_Spline_t::evaluate_second_derivative(double interp_point)
{
	double d2f=0;	
	for (long int i=0; i<num_points; i++) {
		if ( ( interp_point >= points[i] ) && ( interp_point <= points[i+1] ) ) {
			d2f = 2.0 * c[i] + 6.0 * d[i] * ( interp_point - points[i] );
//			fprintf(stdout, "%f\n",d2f);
		}
	}
	return d2f;
}

//----------------------------------------------------------------------------
