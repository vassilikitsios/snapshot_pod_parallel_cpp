/*=========================================================================
 snapshot_pod
 =========================================================================
 
 This program performs the snapshot POD analysis. The snapshots are 
 binary vtk files and the output eigenmodes are written as binary vtk files
 
 =========================================================================*/

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

#include "Utils.h"
#include "PodUtils.h"
#include "input_data.h"

//----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	//............................................................................
	cout << "Cleaning results directory..." << endl;
	string results_directory("./results");	rmdir(results_directory.c_str());		mkdir(results_directory.c_str(), 0777);	
	
	//............................................................................
	cout << "Reading input deck..." << endl;	
	InputData_t input_deck;

	//............................................................................
	POD_t pod(input_deck.temporal_correlation_type, input_deck.normalisation, input_deck.dt, input_deck.two_dimensional,
			  input_deck.snapshot_list_filename, input_deck.snapshot_dir, input_deck.file_name_with_cell_volumes,
			  input_deck.x_min, input_deck.x_max, input_deck.y_min, input_deck.y_max, input_deck.z_min, input_deck.z_max,
			  "none");		// Needed for reconstruction code
	
	pod.allocate_memory();
	pod.report_memory_usage();
	
	pod.read_temporal_correlations(input_deck.restart);
	pod.read_snapshots();

	pod.calculate_temporal_mean();
	pod.write_mean_field();
	
	pod.subtract_temporal_mean_from_snapshots();
	pod.write_first_snapshot_unsteady_component();
	
	pod.calculate_temporal_rms();
	pod.write_rms_field();
	
	pod.calculate_spatial_rms_history();
	pod.write_spatial_rms_history();
	
	pod.calculate_spatial_average_of_temporal_rms();
	pod.write_rms_spatial_average();
	
	pod.calculate_and_write_temporal_correlations();
	
	pod.calculate_and_scale_eigenvalues(input_deck.check_eigensolution);
	pod.calculate_number_of_modes_for_threshhold_energy(input_deck.threshhold_energy,
		input_deck.min_number_of_eigenvalues, input_deck.max_number_of_eigenvalues);
	pod.write_eigenvalues();
	
	pod.calculate_and_scale_temporal_modes();
	if (input_deck.test_orthogonality) pod.test_temporal_orthogonality();
	pod.write_temporal_modes();
	
	pod.calculate_and_write_spatial_modes();
	
	pod.write_summary();	

	//............................................................................
	return 0;
}

//----------------------------------------------------------------------------
