
#include <math.h>
#include <iostream>

#include "Utils.h"
#include "Snapshot_Data_t.h"
#include "real_sym_eig_solver.h"
#include "Spatial_Mode_Calculator_t.h"

#define COUNTER_NOTIFY_INTERVAL		1000

//----------------------------------------------------------------------------
Spatial_Mode_Calculator_t::Spatial_Mode_Calculator_t()
{
}

//----------------------------------------------------------------------------
Spatial_Mode_Calculator_t::~Spatial_Mode_Calculator_t()
{	
	free_1d_array_Snapshot_Data_t(spatial_mode);
}

//----------------------------------------------------------------------------
void Spatial_Mode_Calculator_t::set_problem_parameters(char *temporal_correlation_type, char *normalisation, long int N_points, long int N_snapshots)
{
	this->temporal_correlation_type = temporal_correlation_type;
	this->normalisation = normalisation;
	this->N_points = N_points;
	this->N_snapshots = N_snapshots;
	spatial_mode = malloc_1d_array_Snapshot_Data_t(N_points);
}

//----------------------------------------------------------------------------
void Spatial_Mode_Calculator_t::calculate_spatial_mode(long int mode_num, double **A, Snapshot_Data_t **snapshots, double *eigenvalues) 
{
	if (strcmp(temporal_correlation_type,"ke")==0) {
		calculate_spatial_mode_ke(mode_num, A, snapshots);
		
	} else if (strcmp(temporal_correlation_type,"p")==0) {
		calculate_spatial_mode_p(mode_num, A, snapshots);
		
	} else if (strcmp(temporal_correlation_type,"q")==0) {
		calculate_spatial_mode_q(mode_num, A, snapshots);
		
	} else if (strcmp(temporal_correlation_type,"u")==0) {
		calculate_spatial_mode_u(mode_num, A, snapshots);
		
	} else if (strcmp(temporal_correlation_type,"v")==0) {
		calculate_spatial_mode_v(mode_num, A, snapshots);
		
	} else if (strcmp(temporal_correlation_type,"w")==0) {
		calculate_spatial_mode_w(mode_num, A, snapshots);
	}
	
	for (long int k=0; k<N_points; k++)	{
		spatial_mode[k].u = spatial_mode[k].u / N_snapshots / eigenvalues[mode_num];
		spatial_mode[k].v = spatial_mode[k].v / N_snapshots / eigenvalues[mode_num];
		spatial_mode[k].w = spatial_mode[k].w / N_snapshots / eigenvalues[mode_num];
		spatial_mode[k].p = spatial_mode[k].p / N_snapshots / eigenvalues[mode_num];
	}
}

//----------------------------------------------------------------------------
void Spatial_Mode_Calculator_t::calculate_spatial_mode_ke(long int mode_num, double **A, Snapshot_Data_t **snapshots) 
{	
	for (long int k=0; k<N_points; k++)	{
		spatial_mode[k].u = 0.0;
		spatial_mode[k].v = 0.0;
		spatial_mode[k].w = 0.0;
		spatial_mode[k].p = 0.0;
		
		for (int j=0; j<N_snapshots; j++)	{						// For each snapshot
			spatial_mode[k].u = spatial_mode[k].u + A[j][mode_num] * snapshots[k][j].u;
			spatial_mode[k].v = spatial_mode[k].v + A[j][mode_num] * snapshots[k][j].v;
			spatial_mode[k].w = spatial_mode[k].w + A[j][mode_num] * snapshots[k][j].w;
		}	
	}
}

//----------------------------------------------------------------------------
void Spatial_Mode_Calculator_t::calculate_spatial_mode_u(long int mode_num, double **A, Snapshot_Data_t **snapshots) 
{	
	for (long int k=0; k<N_points; k++)	{
		spatial_mode[k].u = 0.0;
		spatial_mode[k].v = 0.0;
		spatial_mode[k].w = 0.0;
		spatial_mode[k].p = 0.0;
		
		for (int j=0; j<N_snapshots; j++)
			spatial_mode[k].u = spatial_mode[k].u + A[j][mode_num] * snapshots[k][j].u;
	}
}

//----------------------------------------------------------------------------
void Spatial_Mode_Calculator_t::calculate_spatial_mode_v(long int mode_num, double **A, Snapshot_Data_t **snapshots) 
{	
	for (long int k=0; k<N_points; k++)	{
		spatial_mode[k].u = 0.0;
		spatial_mode[k].v = 0.0;
		spatial_mode[k].w = 0.0;
		spatial_mode[k].p = 0.0;
		
		for (int j=0; j<N_snapshots; j++)
			spatial_mode[k].v = spatial_mode[k].v + A[j][mode_num] * snapshots[k][j].v;
	}
}

//----------------------------------------------------------------------------
void Spatial_Mode_Calculator_t::calculate_spatial_mode_w(long int mode_num, double **A, Snapshot_Data_t **snapshots) 
{	
	for (long int k=0; k<N_points; k++)	{
		spatial_mode[k].u = 0.0;
		spatial_mode[k].v = 0.0;
		spatial_mode[k].w = 0.0;
		spatial_mode[k].p = 0.0;
		
		for (int j=0; j<N_snapshots; j++)
			spatial_mode[k].w = spatial_mode[k].w + A[j][mode_num] * snapshots[k][j].w;
	}
}

//----------------------------------------------------------------------------
void Spatial_Mode_Calculator_t::calculate_spatial_mode_p(long int mode_num, double **A, Snapshot_Data_t **snapshots) 
{		
	for (long int k=0; k<N_points; k++)	{
		spatial_mode[k].u = 0.0;
		spatial_mode[k].v = 0.0;
		spatial_mode[k].w = 0.0;
		spatial_mode[k].p = 0.0;
		
		for (int j=0; j<N_snapshots; j++)	{						// For each snapshot
			spatial_mode[k].u = 0.0;
			spatial_mode[k].v = 0.0;
			spatial_mode[k].w = 0.0;
			spatial_mode[k].p = spatial_mode[k].p + A[j][mode_num] * snapshots[k][j].p;
		}	
	}
}

//----------------------------------------------------------------------------
void Spatial_Mode_Calculator_t::calculate_spatial_mode_q(long int mode_num, double **A, Snapshot_Data_t **snapshots) 
{			
	for (long int k=0; k<N_points; k++)	{
		spatial_mode[k].u = 0.0;
		spatial_mode[k].v = 0.0;
		spatial_mode[k].w = 0.0;
		spatial_mode[k].p = 0.0;
		
		for (int j=0; j<N_snapshots; j++)	{						// For each snapshot
			spatial_mode[k].u = spatial_mode[k].u + A[j][mode_num] * snapshots[k][j].u;
			spatial_mode[k].v = spatial_mode[k].v + A[j][mode_num] * snapshots[k][j].v;
			spatial_mode[k].w = spatial_mode[k].w + A[j][mode_num] * snapshots[k][j].w;
			spatial_mode[k].p = spatial_mode[k].p + A[j][mode_num] * snapshots[k][j].p;
		}	
	}
}

//----------------------------------------------------------------------------
double Spatial_Mode_Calculator_t::test_spatial_orthogonality(Snapshot_Data_t **spatial_modes, double *cell_volumes, 
															 Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average,
															 long int i, long int j)
{
	double norm;
	if (strcmp(normalisation,"local")==0) {
		norm = test_spatial_orthogonality_norm_local(spatial_modes, cell_volumes, rms_field, rms_spatial_average, i, j);
		
	} else if (strcmp(normalisation,"global")==0) {
		norm = test_spatial_orthogonality_norm_global(spatial_modes, cell_volumes, rms_spatial_average, i, j);
		
	} else if (strcmp(normalisation,"none")==0) {
		norm = test_spatial_orthogonality(spatial_modes, cell_volumes, i, j);
	}
	return norm;
}

//----------------------------------------------------------------------------
double Spatial_Mode_Calculator_t::test_spatial_orthogonality_norm_global(Snapshot_Data_t **spatial_modes, double *cell_volumes,
																		 Snapshot_Data_t rms_spatial_average, long int i, long int j) 
{	
	double norm = 0.0, norm_u = 0.0, norm_v = 0.0, norm_w = 0.0, norm_p = 0.0;
	for (long int k=0; k<N_points; k++) {
		norm_u = norm_u + spatial_modes[k][i].u * spatial_modes[k][j].u * cell_volumes[k] / SQ(rms_spatial_average.u);
		norm_v = norm_v + spatial_modes[k][i].v * spatial_modes[k][j].v * cell_volumes[k] / SQ(rms_spatial_average.v);
		norm_w = norm_w + spatial_modes[k][i].w * spatial_modes[k][j].w * cell_volumes[k] / SQ(rms_spatial_average.w);
		norm_p = norm_p + spatial_modes[k][i].p * spatial_modes[k][j].p * cell_volumes[k] / SQ(rms_spatial_average.p);
	}
	norm = norm_u + norm_v + norm_w + norm_p;
	return norm;
}

//----------------------------------------------------------------------------
double Spatial_Mode_Calculator_t::test_spatial_orthogonality_norm_local(Snapshot_Data_t **spatial_modes, double *cell_volumes,
																		Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average, 
																		long int i, long int j) 
{	
	double norm = 0.0, norm_u = 0.0, norm_v = 0.0, norm_w = 0.0, norm_p = 0.0;
	for (long int k=0; k<N_points; k++) {
		norm_u = norm_u + spatial_modes[k][i].u * spatial_modes[k][j].u * cell_volumes[k] / ( SQ(rms_field[k].u) + SQ(rms_spatial_average.u) );
		norm_v = norm_v + spatial_modes[k][i].v * spatial_modes[k][j].v * cell_volumes[k] / ( SQ(rms_field[k].v) + SQ(rms_spatial_average.v) );
		norm_w = norm_w + spatial_modes[k][i].w * spatial_modes[k][j].w * cell_volumes[k] / ( SQ(rms_field[k].w) + SQ(rms_spatial_average.w) );
		norm_p = norm_p + spatial_modes[k][i].p * spatial_modes[k][j].p * cell_volumes[k] / ( SQ(rms_field[k].p) + SQ(rms_spatial_average.p) );
	}
	norm = norm_u + norm_v + norm_w + norm_p;
	return norm;
}

//----------------------------------------------------------------------------
double Spatial_Mode_Calculator_t::test_spatial_orthogonality(Snapshot_Data_t **spatial_modes, double *cell_volumes, long int i, long int j) 
{	
	double norm = 0.0, norm_u = 0.0, norm_v = 0.0, norm_w = 0.0, norm_p = 0.0;
	for (long int k=0; k<N_points; k++) {
		norm_u = norm_u + spatial_modes[k][i].u * spatial_modes[k][j].u * cell_volumes[k];
		norm_v = norm_v + spatial_modes[k][i].v * spatial_modes[k][j].v * cell_volumes[k];
		norm_w = norm_w + spatial_modes[k][i].w * spatial_modes[k][j].w * cell_volumes[k];
		norm_p = norm_p + spatial_modes[k][i].p * spatial_modes[k][j].p * cell_volumes[k];
	}
	norm = norm_u + norm_v + norm_w + norm_p;
	return norm;
}

//----------------------------------------------------------------------------
void Spatial_Mode_Calculator_t::norm_check(long int i, long int j, double norm)
{
	if (i==j) {
		if (sqrt(SQ(norm - 1.0)) > ZERO*1.0e3) 
			fprintf(stdout,"WARNING norm check failed for spatial mode pair (%ld, %ld), should be %g, resulted in %g\n", i, j, 1.0, norm);
	} else {
		if (sqrt(SQ(norm)) > ZERO*1.0e3) 
			fprintf(stdout,"WARNING norm check failed for spatial mode pair (%ld, %ld), should be %g, resulted in %g\n", i, j, 0.0, norm);
	}
}

//----------------------------------------------------------------------------

