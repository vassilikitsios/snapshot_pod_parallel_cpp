/*=========================================================================
 snapshot_pod_mpi
 =========================================================================
 
 This program performs the snapshot POD analysis over a series of cpus decomposed in space. 
 The eigenvalue problem is solved on all processors to avoid additional MPI calls.
 The snapshots are binary vtk files and the output eigenmodes are written as binary vtk files.
 
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

void calculate_spatial_average_of_temporal_rms_MPI(POD_t *pod);
void calculate_spatial_rms_history_MPI(POD_t *pod);
void calculate_and_write_temporal_correlations_MPI(POD_t *pod, int myrank);

void write_first_snapshot_unsteady_component_MPI(POD_t *pod);
void write_mean_field_MPI(POD_t *pod);
void write_rms_field_MPI(POD_t *pod);
void calculate_and_write_spatial_modes_MPI(POD_t *pod);

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
			  "none");		// Needed for reconstruction code
	
	long int N_nodes_per_cpu = (long int) round(pod.get_total_number_of_points()/nprocs)+1;
	long int first_node_number = myrank*N_nodes_per_cpu;
	long int last_node_number = (myrank + 1) * N_nodes_per_cpu - 1;
	if (myrank==nprocs-1) last_node_number = pod.get_total_number_of_points() - 1;
	pod.set_node_number_range(first_node_number, last_node_number);
	
	pod.allocate_memory();
	pod.report_memory_usage(myrank);
	if (myrank==0) start_time = report_time(start_time,"memory allocation");

	pod.read_temporal_correlations(input_deck.restart);
	pod.read_snapshots(myrank);
	if (myrank==0) start_time = report_time(start_time,"read files");
	
	pod.calculate_temporal_mean();
	pod.subtract_temporal_mean_from_snapshots();
	pod.calculate_temporal_rms();
	calculate_spatial_rms_history_MPI(&pod);
	if (myrank==0) pod.write_spatial_rms_history();
	
	calculate_spatial_average_of_temporal_rms_MPI(&pod);
	if (myrank==0) pod.write_rms_spatial_average();
	
	if (myrank==0) start_time = report_time(start_time,"calculate and write statstical measures, and substract mean from snapshots");
	
	calculate_and_write_temporal_correlations_MPI(&pod, myrank);
	if (myrank==0) start_time = report_time(start_time,"calculate and write temporal correlations");
	
	pod.calculate_and_scale_eigenvalues(input_deck.check_eigensolution);
	pod.calculate_number_of_modes_for_threshhold_energy(input_deck.threshhold_energy, 
														input_deck.min_number_of_eigenvalues, input_deck.max_number_of_eigenvalues);
	if (myrank==0) pod.write_eigenvalues();
	if (myrank==0) start_time = report_time(start_time,"calculate and write eigenvalues");
	
	pod.calculate_and_scale_temporal_modes();
	if (myrank==0) pod.write_temporal_modes();
	if (myrank==0) start_time = report_time(start_time,"calculate and write temporal_modes");
	
	write_mean_field_MPI(&pod);
	calculate_and_write_spatial_modes_MPI(&pod);
	write_first_snapshot_unsteady_component_MPI(&pod);
	write_rms_field_MPI(&pod);
	if (myrank==0) start_time = report_time(start_time,"calculate and write spatial_modes");

	if (input_deck.test_orthogonality) {
		pod.test_temporal_orthogonality();
		if (myrank==0) start_time = report_time(start_time,"temporal orthogonality tests");
	}
	
	if (myrank==0) pod.write_summary();
	
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
void calculate_spatial_average_of_temporal_rms_MPI(POD_t *pod)
{		
	double p_this_cpu, u_this_cpu, v_this_cpu, w_this_cpu;
	double p_all_cpus, u_all_cpus, v_all_cpus, w_all_cpus;
	int buffer_count = 1;

	pod->calculate_spatial_average_of_temporal_rms();
	
	u_this_cpu = pod->rms_spatial_average.u;
	v_this_cpu = pod->rms_spatial_average.v;
	w_this_cpu = pod->rms_spatial_average.w;
	p_this_cpu = pod->rms_spatial_average.p;
	
	MPI_Allreduce(&u_this_cpu, &u_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&v_this_cpu, &v_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&w_this_cpu, &w_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&p_this_cpu, &p_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	pod->rms_spatial_average.u = u_all_cpus;
	pod->rms_spatial_average.v = v_all_cpus;
	pod->rms_spatial_average.w = w_all_cpus;
	pod->rms_spatial_average.p = p_all_cpus;
}

//----------------------------------------------------------------------------
void calculate_spatial_rms_history_MPI(POD_t *pod)
{
    double sum_p2_this_cpu, sum_u2_this_cpu, sum_v2_this_cpu, sum_w2_this_cpu;
    double sum_p2_all_cpus, sum_u2_all_cpus, sum_v2_all_cpus, sum_w2_all_cpus;
	int buffer_count = 1;
	
	pod->calculate_spatial_rms_history();
	
	for (long int j=0; j<pod->get_number_of_snapshots(); j++) {
		sum_u2_this_cpu = SQ(pod->spatial_rms_history[j].u);
		sum_v2_this_cpu = SQ(pod->spatial_rms_history[j].v);
		sum_w2_this_cpu = SQ(pod->spatial_rms_history[j].w);
		sum_p2_this_cpu = SQ(pod->spatial_rms_history[j].p);
		
		MPI_Allreduce(&sum_u2_this_cpu, &sum_u2_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&sum_v2_this_cpu, &sum_v2_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&sum_w2_this_cpu, &sum_w2_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&sum_p2_this_cpu, &sum_p2_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		pod->spatial_rms_history[j].u = sqrt(sum_u2_all_cpus);
		pod->spatial_rms_history[j].v = sqrt(sum_v2_all_cpus);
		pod->spatial_rms_history[j].w = sqrt(sum_w2_all_cpus);
		pod->spatial_rms_history[j].p = sqrt(sum_p2_all_cpus);
	}
}

//----------------------------------------------------------------------------
void calculate_and_write_temporal_correlations_MPI(POD_t *pod, int myrank)
{
	double correlation_this_cpu, correlation_all_cpus;
	int buffer_count = 1;
	long int N_snapshots = pod->get_number_of_snapshots();
	
	string output_filename("./results/temporal_correlations.dat");
	FILE *fout;
	if (myrank==0) {																			// Write previously calculated correlations to file
		fout = fopen(output_filename.c_str(),"w");
		fprintf(fout,"#\n# snapshot i, snapshot j, correlation value\n#\n");
		for (long int i=0; i<pod->get_last_restart_snapshot_i(); i++)
			for (long int j=0; j<pod->get_last_restart_snapshot_j(); j++)
				fprintf(fout,"%ld %ld %16.10e\n", i+1, j+1, pod->A[i][j]);
	}
	
	fprintf(stdout,"Calculating temporal correlations... \n");
	for (long int i=pod->get_last_restart_snapshot_i(); i<N_snapshots; i++) {					// Rows of temporal_correlation
		for (long int j=pod->get_last_restart_snapshot_j(); j<N_snapshots; j++) {				// Columns of temporal_correlation
			
			if (j>=i) {
				correlation_this_cpu = pod->calculate_inner_product(i, j);
				MPI_Allreduce(&correlation_this_cpu, &correlation_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				pod->A[i][j] = correlation_all_cpus;
			} else {
				pod->A[i][j] = pod->A[j][i];		// Symmetric matrix
			}
			if ( (i*N_snapshots+j)%COUNTER_NOTIFY_INTERVAL == 0)
				cout << "Correlation completed for " << i*N_snapshots+j << " of " << N_snapshots*N_snapshots << endl;
		}
		if (myrank==0)																			// Write newly calculated correlations to file
			for (long int j=pod->get_last_restart_snapshot_j(); j<N_snapshots; j++)
				fprintf(fout,"%ld %ld %16.10e\n", i+1, j+1, pod->A[i][j]);
	}	
	if ( (N_snapshots*N_snapshots)%COUNTER_NOTIFY_INTERVAL != 0)
		cout << "Correlation completed for " << N_snapshots*N_snapshots << " of " << N_snapshots*N_snapshots << endl << endl;
	
	if (myrank==0) fclose(fout);	
}

//----------------------------------------------------------------------------
void write_mean_field_MPI(POD_t *pod)
{
	if (myrank==0) fprintf(stdout, "Writing mean field...\n");	
	
	vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
	writer->SetInput(pod->first_snapshot);
	writer->SetFileTypeToBinary();
	
	double p_this_cpu, u_this_cpu, v_this_cpu, w_this_cpu;
	double p_all_cpus, u_all_cpus, v_all_cpus, w_all_cpus;
	int buffer_count = 1;
	
	for(long int node=0; node < pod->get_total_number_of_points(); node++) {
		
		if ( (node < pod->get_first_point_number()) || (node > pod->get_last_point_number() ) ) {
			p_this_cpu = 0.0;
			u_this_cpu = 0.0;
			v_this_cpu = 0.0;
			w_this_cpu = 0.0;
		} else {
			p_this_cpu = pod->mean_field[node - pod->get_first_point_number()].p;
			u_this_cpu = pod->mean_field[node - pod->get_first_point_number()].u;
			v_this_cpu = pod->mean_field[node - pod->get_first_point_number()].v;
			w_this_cpu = pod->mean_field[node - pod->get_first_point_number()].w;
		}
		
		MPI_Allreduce(&u_this_cpu, &u_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&v_this_cpu, &v_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&w_this_cpu, &w_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&p_this_cpu, &p_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		pod->first_snapshot->GetPointData()->GetScalars("p")->SetComponent(node, 0, p_all_cpus);
		pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 0, u_all_cpus);
		pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 1, v_all_cpus);
		pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 2, w_all_cpus);
	}
	writer->SetFileName("./results/mean_field.vtk");
	if (myrank==0) writer->Write();
	writer->Delete();	
}

//----------------------------------------------------------------------------
void write_first_snapshot_unsteady_component_MPI(POD_t *pod)
{
	if (myrank==0) fprintf(stdout, "Writing unsteady component of first snapshot... \n");
	
	vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
	writer->SetInput(pod->first_snapshot);
	writer->SetFileTypeToBinary();
	
	double p_this_cpu, u_this_cpu, v_this_cpu, w_this_cpu;
	double p_all_cpus, u_all_cpus, v_all_cpus, w_all_cpus;
	double snapshot_p, snapshot_u, snapshot_v, snapshot_w;
	int buffer_count = 1;
	
	for(int node=0; node < pod->get_total_number_of_points(); node++) {
		
		if ( (node < pod->get_first_point_number()) || (node > pod->get_last_point_number() ) ) {
			p_this_cpu = 0.0;
			u_this_cpu = 0.0;
			v_this_cpu = 0.0;
			w_this_cpu = 0.0;
		} else {
			p_this_cpu = pod->mean_field[node - pod->get_first_point_number()].p;
			u_this_cpu = pod->mean_field[node - pod->get_first_point_number()].u;
			v_this_cpu = pod->mean_field[node - pod->get_first_point_number()].v;
			w_this_cpu = pod->mean_field[node - pod->get_first_point_number()].w;
		}
		
		MPI_Allreduce(&u_this_cpu, &u_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&v_this_cpu, &v_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&w_this_cpu, &w_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&p_this_cpu, &p_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		snapshot_p = pod->first_snapshot->GetPointData()->GetScalars("p")->GetComponent(node, 0);
		snapshot_u = pod->first_snapshot->GetPointData()->GetVectors("u")->GetComponent(node, 0);
		snapshot_v = pod->first_snapshot->GetPointData()->GetVectors("u")->GetComponent(node, 1);
		snapshot_w = pod->first_snapshot->GetPointData()->GetVectors("u")->GetComponent(node, 2);
		
		pod->first_snapshot->GetPointData()->GetScalars("p")->SetComponent(node, 0, snapshot_p - p_all_cpus);
		pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 0, snapshot_u - u_all_cpus);
		pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 1, snapshot_v - v_all_cpus);
		pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 2, snapshot_w - w_all_cpus);
	}
	
	writer->SetFileName("./results/first_snapshot_unsteady_component.vtk");
	if (myrank==0) writer->Write();
	writer->Delete();
}


//----------------------------------------------------------------------------
void write_rms_field_MPI(POD_t *pod)
{
	if (myrank==0) fprintf(stdout, "Writing rms field... \n");	
	
	vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
	writer->SetInput(pod->first_snapshot);
	writer->SetFileTypeToBinary();
	
	double p_this_cpu, u_this_cpu, v_this_cpu, w_this_cpu;
	double p_all_cpus, u_all_cpus, v_all_cpus, w_all_cpus;
	int buffer_count = 1;
	
	for(int node=0; node < pod->get_total_number_of_points(); node++) {
		
		if ( (node < pod->get_first_point_number()) || (node > pod->get_last_point_number() ) ) {
			p_this_cpu = 0.0;
			u_this_cpu = 0.0;
			v_this_cpu = 0.0;
			w_this_cpu = 0.0;
		} else {
			p_this_cpu = pod->rms_field[node - pod->get_first_point_number()].p;
			u_this_cpu = pod->rms_field[node - pod->get_first_point_number()].u;
			v_this_cpu = pod->rms_field[node - pod->get_first_point_number()].v;
			w_this_cpu = pod->rms_field[node - pod->get_first_point_number()].w;
		}
		
		MPI_Allreduce(&u_this_cpu, &u_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&v_this_cpu, &v_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&w_this_cpu, &w_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&p_this_cpu, &p_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		pod->first_snapshot->GetPointData()->GetScalars("p")->SetComponent(node, 0, p_all_cpus);
		pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 0, u_all_cpus);
		pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 1, v_all_cpus);
		pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 2, w_all_cpus);
	}
	writer->SetFileName("./results/rms_field.vtk");
    if (myrank==0) writer->Write();
    writer->Delete();
}

//----------------------------------------------------------------------------
void calculate_and_write_spatial_modes_MPI(POD_t *pod)
{
	if (myrank==0) fprintf(stdout, "Writing spatial modes ... \n");

	FILE *fout;
	if (myrank==0) fout = fopen("./results/spatial_modes.list","w");
	
	string output_filename_string;
	string output_location("./results/");
	string name("spatial_mode_");
	string pad_zeros;
	string mode_number_string;
	string extension(".vtk");
	
	double p_this_cpu, u_this_cpu, v_this_cpu, w_this_cpu;
	double p_all_cpus, u_all_cpus, v_all_cpus, w_all_cpus;
	int buffer_count = 1;
	
	vtkUnstructuredGridWriter *writer;
	
	for (int mode=0; mode<pod->get_number_of_modes(); mode++) {						// For each mode

		// Set filename
		mode_number_string = convert_int_to_string(mode+1);
		if (mode+1>9999) {
			quit_error((char*)"Too many modes to output, increase digits in output file names");
		} else if ( (mode+1<9999) && (mode+1>=1000) ) {
			pad_zeros = "";
		} else if ( (mode+1<1000) && (mode+1>=100) ) {
			pad_zeros = "0";
		} else if ( (mode+1<100) && (mode+1>=10) ) {
			pad_zeros = "00";
		} else if (mode+1<10) {
			pad_zeros = "000";
		}
		output_filename_string = output_location + name + pad_zeros + mode_number_string + extension;
		
		// Set the data
		pod->calculate_and_scale_mode_i(mode);
		for(long int node=0; node < pod->get_total_number_of_points(); node++) {

			if ( (node < pod->get_first_point_number()) || (node > pod->get_last_point_number() ) ) {
				p_this_cpu = 0.0;
				u_this_cpu = 0.0;
				v_this_cpu = 0.0;
				w_this_cpu = 0.0;
			} else {
				p_this_cpu = pod->spatial_mode_calculator.spatial_mode[node - pod->get_first_point_number()].p;
				u_this_cpu = pod->spatial_mode_calculator.spatial_mode[node - pod->get_first_point_number()].u;
				v_this_cpu = pod->spatial_mode_calculator.spatial_mode[node - pod->get_first_point_number()].v;
				w_this_cpu = pod->spatial_mode_calculator.spatial_mode[node - pod->get_first_point_number()].w;
			}
			
			MPI_Allreduce(&u_this_cpu, &u_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&v_this_cpu, &v_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&w_this_cpu, &w_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&p_this_cpu, &p_all_cpus, buffer_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
			pod->first_snapshot->GetPointData()->GetScalars("p")->SetComponent(node, 0, p_all_cpus);
			pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 0, u_all_cpus);
			pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 1, v_all_cpus);
			pod->first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 2, w_all_cpus);
		}

		// Write the file
		writer = vtkUnstructuredGridWriter::New();
		writer->SetInput(pod->first_snapshot);
		writer->SetFileTypeToBinary();
		writer->SetFileName(output_filename_string.c_str());
		if (myrank==0) writer->Write();
		writer->Delete();
		output_filename_string = name + pad_zeros + mode_number_string + extension;
		if (myrank==0) fprintf(fout,"%s\n",output_filename_string.c_str());
	}
	
	if (myrank==0) fclose(fout);
}

//----------------------------------------------------------------------------
