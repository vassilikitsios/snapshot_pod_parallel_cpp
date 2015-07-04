/*=========================================================================
 pod_reconstruction
 =========================================================================
 
 This program performs a POD reconstruction of the snapshots used to create
 the basis.
 
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
			  input_deck.spatial_modes_dir);
	
	pod.allocate_memory();
	pod.report_memory_usage();
	
	pod.read_spatial_modes();
	pod.read_mean_field(input_deck.spatial_modes_dir);
	pod.read_rms_field(input_deck.spatial_modes_dir);
	pod.read_spatial_rms_history(input_deck.spatial_modes_dir);
	pod.write_spatial_rms_history();
	
	if (input_deck.test_orthogonality) 
		pod.test_spatial_orthogonality();
	
	if (input_deck.test_reconstruction) {	
		pod.calculate_L2_reconstruction_error();
		pod.write_reconstruction_L2_error();
		pod.calculate_L_infinite_reconstruction_error();
		pod.write_reconstruction_L_infinite_error();
	 }

	//............................................................................
	return 0;
}

//----------------------------------------------------------------------------
