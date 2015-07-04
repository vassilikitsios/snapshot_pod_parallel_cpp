
#include "Utils.h"
#include "input_data.h"

//-------------------------------------------------------------------------
InputData_t::InputData_t()
{
	ParamListHead_t param_list_head;
	
	param_list_head.get_params_from_file((char*)"snapshot_pod.in");
	param_list_head.print_params();

	temporal_correlation_type = param_list_head.get_string_param((char*)"temporal_correlation_type");
	normalisation = param_list_head.get_string_param((char*)"normalisation");
	
	snapshot_list_filename = param_list_head.get_string_param((char*)"snapshot_list_filename");
	snapshot_dir = param_list_head.get_string_param((char*)"snapshot_dir");
	file_name_with_cell_volumes = param_list_head.get_string_param((char*)"file_name_with_cell_volumes");
	
	dt = param_list_head.get_double_param((char*)"dt");	
	
	char *test_orthogonality_char;
	test_orthogonality_char = param_list_head.get_string_param((char*)"test_orthogonality");
	if (strcmp(test_orthogonality_char,"true")==0) {
		test_orthogonality = true;
	} else if (strcmp(test_orthogonality_char,"false")==0) {
		test_orthogonality = false;
	} else {
		quit_error((char*)"<test_reconstruction> = 'true' or 'false'\n");
	}
	
	char *check_eigensolution_char;
	check_eigensolution_char = param_list_head.get_string_param((char*)"check_eigensolution");
	if (strcmp(check_eigensolution_char,"true")==0) {
		check_eigensolution = true;
	} else if (strcmp(check_eigensolution_char,"false")==0) {
		check_eigensolution = false;
	} else {
		quit_error((char*)"<check_eigensolution> = 'true' or 'false'\n");
	}
	
	restart = param_list_head.get_string_param((char*)"restart");
	
	threshhold_energy = param_list_head.get_double_param((char*)"threshhold_energy");
	min_number_of_eigenvalues = param_list_head.get_int_param((char*)"min_number_of_eigenvalues");
	max_number_of_eigenvalues = param_list_head.get_int_param((char*)"max_number_of_eigenvalues");
	
	two_dimensional = param_list_head.get_string_param((char*)"two_dimensional");		
	
	x_min = param_list_head.get_double_param((char*) "x_min");	
	x_max = param_list_head.get_double_param((char*) "x_max");	
	y_min = param_list_head.get_double_param((char*) "y_min");	
	y_max = param_list_head.get_double_param((char*) "y_max");	
	z_min = param_list_head.get_double_param((char*) "z_min");	
	z_max = param_list_head.get_double_param((char*) "z_max");
}

//-------------------------------------------------------------------------
InputData_t::~InputData_t()
{
}
