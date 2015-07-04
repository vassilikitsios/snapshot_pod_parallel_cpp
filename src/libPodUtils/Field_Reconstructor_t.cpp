
#include <math.h>
#include <iostream>

#include "Utils.h"
#include "Snapshot_Data_t.h"
#include "real_sym_eig_solver.h"
#include "Field_Reconstructor_t.h"

#define COUNTER_NOTIFY_INTERVAL		1000

//----------------------------------------------------------------------------
Field_Reconstructor_t::Field_Reconstructor_t()
{
}

//----------------------------------------------------------------------------
Field_Reconstructor_t::~Field_Reconstructor_t()
{	
}

//----------------------------------------------------------------------------
void Field_Reconstructor_t::set_problem_parameters(char *temporal_correlation_type, long int N_points)
{
	this->temporal_correlation_type = temporal_correlation_type;
	this->N_points = N_points;
}

//----------------------------------------------------------------------------
void Field_Reconstructor_t::reconstruct_field(long int num_modes, double *inner_products, Snapshot_Data_t **spatial_modes, 
											  Snapshot_Data_t *pod_reconstruction)
{	
	// Assuming the spatial modes are zero for the appropriate temporal correlation types this code should work for all cases
	
	for (long int node=0; node<N_points; node++) {
		pod_reconstruction[node].u = 0.0;
		pod_reconstruction[node].v = 0.0;
		pod_reconstruction[node].w = 0.0;
		pod_reconstruction[node].p = 0.0;
		
		for (long int mode=0; mode<=num_modes; mode++) {
			pod_reconstruction[node].u = pod_reconstruction[node].u + inner_products[mode] * spatial_modes[node][mode].u;
			pod_reconstruction[node].v = pod_reconstruction[node].v + inner_products[mode] * spatial_modes[node][mode].v;
			pod_reconstruction[node].w = pod_reconstruction[node].w + inner_products[mode] * spatial_modes[node][mode].w;
			pod_reconstruction[node].p = pod_reconstruction[node].p + inner_products[mode] * spatial_modes[node][mode].p;
		}
	}
}

//----------------------------------------------------------------------------
void Field_Reconstructor_t::calculate_spatial_L2_error(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, Snapshot_Data_t *pod_reconstruction,
													   double *cell_volumes, Snapshot_Data_t **spatial_L2_error)
{	
	if (strcmp(temporal_correlation_type,"ke")==0) {
		calculate_spatial_L2_error_for_ke(num_modes, snapshot_num, snapshots, pod_reconstruction, cell_volumes, spatial_L2_error);
		
	} else if (strcmp(temporal_correlation_type,"p")==0) {
		calculate_spatial_L2_error_for_p(num_modes, snapshot_num, snapshots, pod_reconstruction, cell_volumes, spatial_L2_error);
		
	} else if (strcmp(temporal_correlation_type,"q")==0) {
		calculate_spatial_L2_error_for_q(num_modes, snapshot_num, snapshots, pod_reconstruction, cell_volumes, spatial_L2_error);
		
	} else if (strcmp(temporal_correlation_type,"u")==0) {
		calculate_spatial_L2_error_for_u(num_modes, snapshot_num, snapshots, pod_reconstruction, cell_volumes, spatial_L2_error);
		
	} else if (strcmp(temporal_correlation_type,"v")==0) {
		calculate_spatial_L2_error_for_v(num_modes, snapshot_num, snapshots, pod_reconstruction, cell_volumes, spatial_L2_error);
		
	} else if (strcmp(temporal_correlation_type,"w")==0) {
		calculate_spatial_L2_error_for_w(num_modes, snapshot_num, snapshots, pod_reconstruction, cell_volumes, spatial_L2_error);
	}
}

//----------------------------------------------------------------------------
void Field_Reconstructor_t::calculate_spatial_L2_error_for_u(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
															 Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error)
{	
	for (long int node=0; node<N_points; node++) {
		spatial_L2_error[num_modes][snapshot_num].u = spatial_L2_error[num_modes][snapshot_num].u 
		+ ( SQ(snapshots[node][0].u - pod_reconstruction[node].u) ) * cell_volumes[node];
		spatial_L2_error[num_modes][snapshot_num].v = 0.0;
		spatial_L2_error[num_modes][snapshot_num].w = 0.0;
		spatial_L2_error[num_modes][snapshot_num].p = 0.0;
	}
}

//----------------------------------------------------------------------------
void Field_Reconstructor_t::calculate_spatial_L2_error_for_v(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
															 Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error)
{	
	for (long int node=0; node<N_points; node++) {
		spatial_L2_error[num_modes][snapshot_num].u = 0.0;
		spatial_L2_error[num_modes][snapshot_num].v = spatial_L2_error[num_modes][snapshot_num].v 
		+ ( SQ(snapshots[node][0].v - pod_reconstruction[node].v) ) * cell_volumes[node];   
		spatial_L2_error[num_modes][snapshot_num].w = 0.0;
		spatial_L2_error[num_modes][snapshot_num].p = 0.0;
	}
}

//----------------------------------------------------------------------------
void Field_Reconstructor_t::calculate_spatial_L2_error_for_w(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
															 Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error)
{	
	for (long int node=0; node<N_points; node++) {
		spatial_L2_error[num_modes][snapshot_num].u = 0.0;
		spatial_L2_error[num_modes][snapshot_num].v = 0.0;   
		spatial_L2_error[num_modes][snapshot_num].w = spatial_L2_error[num_modes][snapshot_num].w 
		+ ( SQ(snapshots[node][0].w - pod_reconstruction[node].w) ) * cell_volumes[node];
		spatial_L2_error[num_modes][snapshot_num].p = 0.0;
	}
}

//----------------------------------------------------------------------------
void Field_Reconstructor_t::calculate_spatial_L2_error_for_p(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
															 Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error)
{	
	for (long int node=0; node<N_points; node++) {
		spatial_L2_error[num_modes][snapshot_num].u = 0.0;
		spatial_L2_error[num_modes][snapshot_num].v = 0.0;   
		spatial_L2_error[num_modes][snapshot_num].w = 0.0;
		spatial_L2_error[num_modes][snapshot_num].p = spatial_L2_error[num_modes][snapshot_num].p 
		+ ( SQ(snapshots[node][0].p - pod_reconstruction[node].p) ) * cell_volumes[node];
	}
}

//----------------------------------------------------------------------------
void Field_Reconstructor_t::calculate_spatial_L2_error_for_ke(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
															  Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error)
{	
	for (long int node=0; node<N_points; node++) {
		
		spatial_L2_error[num_modes][snapshot_num].u = spatial_L2_error[num_modes][snapshot_num].u 
		+ ( SQ(snapshots[node][0].u - pod_reconstruction[node].u) ) * cell_volumes[node];
		
		spatial_L2_error[num_modes][snapshot_num].v = spatial_L2_error[num_modes][snapshot_num].v 
		+ ( SQ(snapshots[node][0].v - pod_reconstruction[node].v) ) * cell_volumes[node];   
		
		spatial_L2_error[num_modes][snapshot_num].w = spatial_L2_error[num_modes][snapshot_num].w 
		+ ( SQ(snapshots[node][0].w - pod_reconstruction[node].w) ) * cell_volumes[node];

		spatial_L2_error[num_modes][snapshot_num].p = 0.0;
	}
}

//----------------------------------------------------------------------------
void Field_Reconstructor_t::calculate_spatial_L2_error_for_q(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
															 Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error)
{	
	for (long int node=0; node<N_points; node++) {
		spatial_L2_error[num_modes][snapshot_num].u = spatial_L2_error[num_modes][snapshot_num].u 
		+ ( SQ(snapshots[node][0].u - pod_reconstruction[node].u) ) * cell_volumes[node];
		spatial_L2_error[num_modes][snapshot_num].v = spatial_L2_error[num_modes][snapshot_num].v 
		+ ( SQ(snapshots[node][0].v - pod_reconstruction[node].v) ) * cell_volumes[node];   
		spatial_L2_error[num_modes][snapshot_num].w = spatial_L2_error[num_modes][snapshot_num].w 
		+ ( SQ(snapshots[node][0].w - pod_reconstruction[node].w) ) * cell_volumes[node];
		spatial_L2_error[num_modes][snapshot_num].p = spatial_L2_error[num_modes][snapshot_num].p 
		+ ( SQ(snapshots[node][0].p - pod_reconstruction[node].p) ) * cell_volumes[node];
	}
}

//----------------------------------------------------------------------------
