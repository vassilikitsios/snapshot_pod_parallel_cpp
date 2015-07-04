/*=========================================================================
 pod_reconstruction_mpi
 =========================================================================
 
 This program performs a POD reconstruction of the snapshots used to create
 the basis.
 
 =========================================================================*/

#include "mpi.h"

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

#include "Utils.h"
#include "PodUtils.h"
#include "input_data.h"

#define COUNTER_NOTIFY_INTERVAL		10000

//----------------------------------------------------------------------------
// mpi helper functions
//----------------------------------------------------------------------------
double report_time(double start_time, string message);

void test_spatial_orthogonality_MPI(POD_t *pod);
void calculate_L2_reconstruction_error_MPI(POD_t *pod, bool two_dimensional, int myrank);

int nprocs, myrank;

//----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	//............................................................................
	// Initialising MPI
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) quit_error_mpi(myrank, nprocs, (char*)"MPI didn't initialise...");
	if (MPI_Comm_size(MPI_COMM_WORLD, &nprocs) != MPI_SUCCESS) quit_error_mpi(myrank, nprocs, (char*)"Couldn't get nprocs...");
	if (MPI_Comm_rank(MPI_COMM_WORLD, &myrank) != MPI_SUCCESS) quit_error_mpi(myrank, nprocs, (char*)"Couldn't get myrank...");
	double start_time = MPI_Wtime();
	
	//............................................................................
	if (myrank==0) fprintf(stdout, "Cleaning results directory...\n");
	string results_directory("./results");	
	if (myrank==0) rmdir(results_directory.c_str());
	if (myrank==0) mkdir(results_directory.c_str(), 0777);	
	
	//............................................................................
	if (myrank==0) fprintf(stdout, "Reading input deck...\n");
	InputData_t input_deck;
	
	//............................................................................
	POD_t pod(input_deck.temporal_correlation_type, input_deck.normalisation, input_deck.dt, input_deck.two_dimensional,
			  input_deck.snapshot_list_filename, input_deck.snapshot_dir, input_deck.file_name_with_cell_volumes,
			  input_deck.x_min, input_deck.x_max, input_deck.y_min, input_deck.y_max, input_deck.z_min, input_deck.z_max,
			  input_deck.spatial_modes_dir);
	
	long int N_nodes_per_cpu = (long int) round(pod.get_total_number_of_points()/nprocs)+1;
	long int first_node_number = myrank*N_nodes_per_cpu;
	long int last_node_number = (myrank + 1) * N_nodes_per_cpu - 1;
	if (myrank==nprocs-1) last_node_number = pod.get_total_number_of_points() - 1;
	pod.set_node_number_range(first_node_number, last_node_number);
	
	pod.allocate_memory();
	pod.report_memory_usage(myrank);
	if (myrank==0) start_time = report_time(start_time,"memory allocation");
	
	pod.read_spatial_modes(myrank);
	pod.read_mean_field(input_deck.spatial_modes_dir, myrank);
	pod.read_rms_field(input_deck.spatial_modes_dir, myrank);
	pod.read_spatial_rms_history(input_deck.spatial_modes_dir);
	if (myrank==0)pod.write_spatial_rms_history();
	
	if (input_deck.test_orthogonality) {
		test_spatial_orthogonality_MPI(&pod);
		if (myrank==0) start_time = report_time(start_time,"spatial orthogonality tests");
	}
	
	if (input_deck.test_reconstruction) {	
		calculate_L2_reconstruction_error_MPI(&pod, input_deck.two_dimensional, myrank);
		if (myrank==0) {
			pod.write_reconstruction_L2_error();
			pod.calculate_L_infinite_reconstruction_error();
			pod.write_reconstruction_L_infinite_error();
			start_time = report_time(start_time,"reconstruction tests");
		}
	}
		
	//............................................................................
	// finalising MPI
	if (MPI_Finalize() != MPI_SUCCESS) quit_error_mpi(myrank, nprocs, (char*)"MPI didn't finalize...");

	//............................................................................
	return 0;
}

//----------------------------------------------------------------------------
double report_time(double start_time, string message)
{
	double end_time = MPI_Wtime(); 
	double elapsed_time = end_time - start_time;
	fprintf(stdout,"DURATION %s: sec = %f, min = %f, hours = %f\n", message.c_str(), elapsed_time, elapsed_time/60.0, elapsed_time/3600.0);
	return end_time;
}

//----------------------------------------------------------------------------
void test_spatial_orthogonality_MPI(POD_t *pod)
{	
	double norm_this_cpu, norm_all_cpus;
	int buffer_count = 1;
	for (long int i=0; i<pod->get_number_of_modes(); i++) {
		for (long int j=0; j<pod->get_number_of_modes(); j++) {			
			norm_this_cpu = pod->test_spatial_orthogonality(i, j);
			MPI_Allreduce(&norm_this_cpu, &norm_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			if (myrank==0) pod->spatial_mode_calculator.norm_check(i, j, norm_all_cpus);
		}
	}
}

//----------------------------------------------------------------------------
void calculate_L2_reconstruction_error_MPI(POD_t *pod, bool two_dimensional, int myrank)
{
	double inner_product_this_cpu, p_this_cpu, u_this_cpu, v_this_cpu, w_this_cpu;
	double inner_product_all_cpus, p_all_cpus, u_all_cpus, v_all_cpus, w_all_cpus;
	double ke;
	int buffer_count = 1;
	long int first_node_number = pod->get_first_point_number();
	long int last_node_number = pod->get_last_point_number();
	long int num_snapshots = pod->get_number_of_snapshots();
	long int num_modes = pod->get_number_of_modes();
	long int num_points = pod->get_number_of_points();

	double *inner_products = malloc_1d_array_double(num_points);
	Snapshot_Data_t **snapshot = malloc_2d_array_Snapshot_Data_t(num_points,1);
	Snapshot_Data_t *pod_reconstruction = malloc_1d_array_Snapshot_Data_t(num_points);
	vtkUnstructuredGridReader *reader;
	
	for (long int snapshot_num=0; snapshot_num < num_snapshots; snapshot_num++) {

		// Read snapshot
		fprintf(stdout, "Reading file %ld of %ld name: %s on processor %d for data\n", 
				snapshot_num+1, num_snapshots, pod->list_of_snapshots[snapshot_num].c_str(), myrank);
		reader = vtkUnstructuredGridReader::New();
		reader->SetFileName(pod->list_of_snapshots[snapshot_num].c_str());
		reader->ReadAllScalarsOn();		
		reader->ReadAllVectorsOn();
		reader->Update();
		for(long int node = first_node_number; node <= last_node_number; node++) {
			snapshot[node-first_node_number][0].p = reader->GetOutput()->GetPointData()->GetScalars("p")->GetComponent(node,0) 
			- pod->mean_field[node-first_node_number].p;
			snapshot[node-first_node_number][0].u = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,0)
			- pod->mean_field[node-first_node_number].u;
			snapshot[node-first_node_number][0].v = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,1)
			- pod->mean_field[node-first_node_number].v;
			snapshot[node-first_node_number][0].w = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,2)
			- pod->mean_field[node-first_node_number].w;
			if(two_dimensional) snapshot[node-first_node_number][0].w = 0.0;
		}
		reader->Delete();
		
		// Reconstruct and scale each mode
		ke = SQ(pod->spatial_rms_history[snapshot_num].u) + SQ(pod->spatial_rms_history[snapshot_num].v) + SQ(pod->spatial_rms_history[snapshot_num].w);
		for (long int mode_num=0; mode_num < num_modes; mode_num++) {
			inner_product_this_cpu = pod->inner_product_calculator.calculate_inner_product(snapshot, pod->spatial_modes, pod->rms_field, 
																						   pod->rms_spatial_average, pod->cell_volumes, 0, mode_num);
			MPI_Allreduce(&inner_product_this_cpu, &inner_product_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			inner_products[mode_num] = inner_product_all_cpus;
			
			pod->field_reconstructor.reconstruct_field(mode_num, inner_products, pod->spatial_modes, pod_reconstruction);
			pod->field_reconstructor.calculate_spatial_L2_error(mode_num, snapshot_num, snapshot, pod_reconstruction, pod->cell_volumes, pod->spatial_L2_error);
			
			u_this_cpu = pod->spatial_L2_error[mode_num][snapshot_num].u;
			v_this_cpu = pod->spatial_L2_error[mode_num][snapshot_num].v;
			w_this_cpu = pod->spatial_L2_error[mode_num][snapshot_num].w;
			p_this_cpu = pod->spatial_L2_error[mode_num][snapshot_num].p;
			
			MPI_Allreduce(&u_this_cpu, &u_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&v_this_cpu, &v_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&w_this_cpu, &w_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&p_this_cpu, &p_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
			pod->spatial_L2_error[mode_num][snapshot_num].u = u_all_cpus / ke;
			pod->spatial_L2_error[mode_num][snapshot_num].v = v_all_cpus / ke;
			pod->spatial_L2_error[mode_num][snapshot_num].w = w_all_cpus / ke;
			if (pod->two_dimensional) pod->spatial_L2_error[mode_num][snapshot_num].w = 0.0;
			pod->spatial_L2_error[mode_num][snapshot_num].p = p_all_cpus / ke;
		}
	}

	free_2d_array_Snapshot_Data_t(snapshot, num_points);
	free_1d_array_Snapshot_Data_t(pod_reconstruction);
	free_1d_array_double(inner_products);
}

//----------------------------------------------------------------------------
