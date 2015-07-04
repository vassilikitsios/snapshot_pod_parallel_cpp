
#include <math.h>
#include <iostream>

#include "Utils.h"
#include "Snapshot_Data_t.h"
#include "Inner_Product_Calculator_t.h"

#define COUNTER_NOTIFY_INTERVAL		1000

//----------------------------------------------------------------------------
Inner_Product_Calculator_t::Inner_Product_Calculator_t()
{
}

//----------------------------------------------------------------------------
Inner_Product_Calculator_t::~Inner_Product_Calculator_t()
{	
}

//----------------------------------------------------------------------------
void Inner_Product_Calculator_t::set_problem_parameters(char *temporal_correlation_type, char *normalisation, bool two_dimensional, long int N_points)
{
	this->temporal_correlation_type = temporal_correlation_type;
	this->normalisation = normalisation;
	this->two_dimensional = two_dimensional;
	this->N_points = N_points;	
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product(Snapshot_Data_t **a, Snapshot_Data_t **b, 
														   Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average,
														   double *cell_volumes, long int i, long int j)
{	
	double correlation = 0.0;
		
	if (strcmp(temporal_correlation_type,"u")==0) {
		if (strcmp(normalisation,"none")==0) {
			correlation = calculate_inner_product_for_u(a, b, cell_volumes, i, j);
		} else if (strcmp(normalisation,"local")==0) {
			correlation = calculate_inner_product_for_u_local_norm(a, b, rms_field, rms_spatial_average, cell_volumes, i, j);
		}
		
	} else if (strcmp(temporal_correlation_type,"v")==0) {
		if (strcmp(normalisation,"none")==0) {
			correlation = calculate_inner_product_for_v(a, b, cell_volumes, i, j);
		} else if (strcmp(normalisation,"local")==0) {
			correlation = calculate_inner_product_for_v_local_norm(a, b, rms_field, rms_spatial_average, cell_volumes, i, j);
		}
		
	} else if ( (strcmp(temporal_correlation_type,"w")==0) && (!two_dimensional) ) {
		if (strcmp(normalisation,"none")==0) {
			correlation = calculate_inner_product_for_w(a, b, cell_volumes, i, j);
		} else if (strcmp(normalisation,"local")==0) {
			correlation = calculate_inner_product_for_w_local_norm(a, b, rms_field, rms_spatial_average, cell_volumes, i, j);
		}
		
	} else if (strcmp(temporal_correlation_type,"p")==0) {
		if (strcmp(normalisation,"none")==0) {
			correlation = calculate_inner_product_for_p(a, b, cell_volumes, i, j);
		} else if (strcmp(normalisation,"local")==0) {
			correlation = calculate_inner_product_for_p_local_norm(a, b, rms_field, rms_spatial_average, cell_volumes, i, j);
		}
		
	} else if (strcmp(temporal_correlation_type,"ke")==0) {
		if (strcmp(normalisation,"none")==0) {
			correlation = calculate_inner_product_for_ke(a, b, cell_volumes, i, j);
		} else if (strcmp(normalisation,"local")==0) {
			correlation = calculate_inner_product_for_ke_local_norm(a, b, rms_field, rms_spatial_average, cell_volumes, i, j);
		} else if (strcmp(normalisation,"global")==0) {
			correlation = calculate_inner_product_for_ke_global_norm(a, b, rms_spatial_average, cell_volumes, i, j);
		}
		
	} else if (strcmp(temporal_correlation_type,"q")==0) {
		if (strcmp(normalisation,"none")==0) {
			correlation = calculate_inner_product_for_q(a, b, cell_volumes, i, j);
		} else if (strcmp(normalisation,"local")==0) {
			correlation = calculate_inner_product_for_q_local_norm(a, b, rms_field, rms_spatial_average, cell_volumes, i, j);
		} else if (strcmp(normalisation,"global")==0) {
			correlation = calculate_inner_product_for_q_global_norm(a, b, rms_spatial_average, cell_volumes, i, j);
		}
	}
	
	return correlation;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_u(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) inner_product = inner_product + ( a[k][i].u * b[k][j].u ) * cell_volumes[k];
	return inner_product;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_u_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, 
																			Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average, 
																			double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) inner_product = inner_product + ( a[k][i].u * b[k][j].u ) * cell_volumes[k]/ ( SQ(rms_field[k].u) + SQ(rms_spatial_average.u));
	return inner_product;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_v(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) inner_product = inner_product + ( a[k][i].v * b[k][j].v ) * cell_volumes[k];
	return inner_product;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_v_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, 
																			Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average, 
																			double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) inner_product = inner_product + ( a[k][i].v * b[k][j].v ) * cell_volumes[k] / ( SQ(rms_field[k].v) + SQ(rms_spatial_average.v));
	return inner_product;
	
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_w(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) inner_product = inner_product + ( a[k][i].w * b[k][j].w ) * cell_volumes[k];
	return inner_product;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_w_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, 
																			Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average, 
																			double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) inner_product = inner_product + ( a[k][i].w * b[k][j].w ) * cell_volumes[k] / ( SQ(rms_field[k].w) + SQ(rms_spatial_average.w));
	return inner_product;
	
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_p(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) inner_product = inner_product + ( a[k][i].p * b[k][j].p ) * cell_volumes[k];
	return inner_product;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_p_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, 
																			Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average, 
																			double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) inner_product = inner_product + ( a[k][i].p * b[k][j].p ) * cell_volumes[k]
		/ ( SQ(rms_field[k].p) + SQ(rms_spatial_average.p));
	return inner_product;
	
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_ke(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j)
{
	double inner_product_u, inner_product_v, inner_product_w;
	inner_product_u = calculate_inner_product_for_u(a, b, cell_volumes, i, j);
	inner_product_v = calculate_inner_product_for_v(a, b, cell_volumes, i, j);
	inner_product_w = calculate_inner_product_for_w(a, b, cell_volumes, i, j);
	return inner_product_u + inner_product_v + inner_product_w;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_ke_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, 
																			 Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average, 
																			 double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) {
//		inner_product = inner_product + ( a[k][i].u * b[k][j].u + a[k][i].v * b[k][j].v + a[k][i].w * b[k][j].w) * cell_volumes[k] / ( SQ(rms_field[k].u) + SQ(rms_field[k].v) + SQ(rms_field[k].w) + SQ(rms_spatial_average.u) + SQ(rms_spatial_average.v) + SQ(rms_spatial_average.w) );
		
		inner_product = inner_product + ( a[k][i].u * b[k][j].u / ( SQ(rms_field[k].u) + SQ(rms_spatial_average.u) )
										 + a[k][i].v * b[k][j].v / ( SQ(rms_field[k].v) + SQ(rms_spatial_average.v) )
										 + a[k][i].w * b[k][j].w / ( SQ(rms_field[k].w) + SQ(rms_spatial_average.w) ) ) * cell_volumes[k];
	}
	return inner_product;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_ke_global_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t rms_spatial_average,
																			 double *cell_volumes, long int i, long int j)
{
	double inner_product_u, inner_product_v, inner_product_w;
	inner_product_u = calculate_inner_product_for_u(a, b, cell_volumes, i, j);
	inner_product_v = calculate_inner_product_for_v(a, b, cell_volumes, i, j);
	inner_product_w = calculate_inner_product_for_w(a, b, cell_volumes, i, j);
	return inner_product_u/SQ(rms_spatial_average.u) + inner_product_v/SQ(rms_spatial_average.v) + inner_product_w/SQ(rms_spatial_average.w);
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_q(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j)
{
	double inner_product_ke, inner_product_p;	
	inner_product_ke = calculate_inner_product_for_ke(a, b, cell_volumes, i, j);
	inner_product_p = calculate_inner_product_for_p(a, b, cell_volumes, i, j);
	return inner_product_ke + inner_product_p;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_q_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, 
																			Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average, 
																			double *cell_volumes, long int i, long int j)
{
	double inner_product = 0.0;
	for (int k=0; k<N_points; k++) 
		inner_product = inner_product + ( a[k][i].u * b[k][j].u + a[k][i].v * b[k][j].v + a[k][i].w * b[k][j].w) * SQ(rms_field[k].p) * cell_volumes[k]
		+ a[k][i].p * b[k][j].p * ( SQ(rms_field[k].u) + SQ(rms_field[k].v) + SQ(rms_field[k].w) ) * cell_volumes[k];
	return inner_product;
}

//----------------------------------------------------------------------------
double Inner_Product_Calculator_t::calculate_inner_product_for_q_global_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t rms_spatial_average,
																			 double *cell_volumes, long int i, long int j)
{
	double inner_product_u, inner_product_v, inner_product_w, inner_product_p;
	inner_product_u = calculate_inner_product_for_u(a, b, cell_volumes, i, j);
	inner_product_v = calculate_inner_product_for_v(a, b, cell_volumes, i, j);
	inner_product_w = calculate_inner_product_for_w(a, b, cell_volumes, i, j);
	inner_product_p = calculate_inner_product_for_p(a, b, cell_volumes, i, j);
	return inner_product_u/SQ(rms_spatial_average.u) 
	+ inner_product_v/SQ(rms_spatial_average.v) 
	+ inner_product_w/SQ(rms_spatial_average.w) 
	+ inner_product_p/(rms_spatial_average.p);
//	return ( inner_product_u + inner_product_v + inner_product_w ) / ( SQ(rms_spatial_average.u) + SQ(rms_spatial_average.v) + SQ(rms_spatial_average.w) ) + inner_product_p / SQ(rms_spatial_average.p);
}

//----------------------------------------------------------------------------
