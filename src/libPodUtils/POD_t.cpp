
#include "Utils.h"
#include "Snapshot_Data_t.h"
#include "real_sym_eig_solver.h"

#include "Inner_Product_Calculator_t.h"
#include "Spatial_Mode_Calculator_t.h"
#include "Field_Reconstructor_t.h"
#include "POD_t.h"

#define COUNTER_NOTIFY_INTERVAL		1000

//----------------------------------------------------------------------------
POD_t::POD_t(char *temporal_correlation_type, char *normalisation, double dt, char *two_dimensional_char,
			 char *snapshot_list_filename, char *snapshot_dir, char *file_name_with_cell_volumes,
			 double x_min, double x_max, double y_min, double y_max, double z_min, double z_max,
			 char *spatial_modes_dir)
{
	if ( (strcmp(temporal_correlation_type,"u")!=0) && 
			(strcmp(temporal_correlation_type,"v")!=0) &&
			(strcmp(temporal_correlation_type,"w")!=0) &&
			(strcmp(temporal_correlation_type,"ke")!=0) && 
			(strcmp(temporal_correlation_type,"p")!=0) && 
		(strcmp(temporal_correlation_type,"q")!=0) ) {
		quit_error((char*) "<temporal_correlation_type> = 'u', 'v', 'w','ke', 'p', 'q'\n");
	} else {
		this->temporal_correlation_type = temporal_correlation_type;
	}

	if ( (strcmp(normalisation,"none")!=0) && 
		(strcmp(normalisation,"local")!=0) &&
		(strcmp(normalisation,"global")!=0) ) {
		quit_error((char*) "<temporal_correlation_type> = 'none', 'local', 'global'\n");
	} else {
		this->normalisation = normalisation;
	}
	
	if (strcmp(two_dimensional_char,"true")==0) {
		two_dimensional = true;
	} else if (strcmp(two_dimensional_char,"false")==0) {
		two_dimensional = false;
	} else {
		quit_error((char*) "<two_dimensional> = 'true' or 'false'\n");
	}

	this->dt = dt;
	this->x_min = x_min;
	this->x_max = x_max;
	this->y_min = y_min;
	this->y_max = y_max;
	this->z_min = z_min;
	this->z_max = z_max;
	
	read_list_of_snapshots(snapshot_list_filename, snapshot_dir);
	this->N_snapshots = list_of_snapshots.size();
	first_snapshot_number = 0;
	last_snapshot_number = N_snapshots - 1;
	
	if (strcmp(spatial_modes_dir,"none")!=0) {
		read_list_of_spatial_modes(spatial_modes_dir);
		this->N_modes = list_of_spatial_modes.size();	
	} else {
		this->N_modes = 0;
	}
	
	first_snapshot = vtkUnstructuredGrid::New();
	read_and_save_first_snapshot(file_name_with_cell_volumes);
	this->N_points = first_snapshot->GetNumberOfPoints();
	first_node_number = 0;
	last_node_number = N_points - 1;
	last_restart_snapshot_i = 0;
	last_restart_snapshot_j = 0;
	
	snapshots_allocated = false;
	spatial_modes_allocated = false;
}

//----------------------------------------------------------------------------
void POD_t::allocate_memory() 
{
	fprintf(stdout, "Allocating memory ... \n");
	
	inner_product_calculator.set_problem_parameters(temporal_correlation_type, normalisation, two_dimensional, N_points);
	spatial_mode_calculator.set_problem_parameters(temporal_correlation_type, normalisation, N_points, N_snapshots);
	field_reconstructor.set_problem_parameters(temporal_correlation_type, N_points);

	if (N_modes==0) {
		snapshots = malloc_2d_array_Snapshot_Data_t(N_points, N_snapshots);
		A = malloc_2d_array_double(N_snapshots, N_snapshots);
		eigenvalues = malloc_1d_array_double(N_snapshots);
		scaled_eigenvalues = malloc_1d_array_double(N_snapshots);
		snapshots_allocated = true;
	} else {
		spatial_modes = malloc_2d_array_Snapshot_Data_t(N_points, N_modes);
		spatial_L2_error = malloc_2d_array_Snapshot_Data_t(N_snapshots, N_snapshots);
		max_spatial_L2_error = malloc_1d_array_Snapshot_Data_t(N_snapshots);	
		spatial_modes_allocated = true;
	}
	
	cell_volumes = malloc_1d_array_double(N_points);
	mean_field = malloc_1d_array_Snapshot_Data_t(N_points);
	rms_field = malloc_1d_array_Snapshot_Data_t(N_points);
	spatial_rms_history = malloc_1d_array_Snapshot_Data_t(N_snapshots);
}

//----------------------------------------------------------------------------
void POD_t::report_memory_usage(int myrank)
{
	if (snapshots_allocated) {
		fprintf(stdout, 
				"MEMORY USAGE on myrank = %d\n \
				N_points = %ld\n \
				N_snapshots = %ld\n \
				sizeof(double*) = %ld bytes\n \
				sizeof(double) = %ld bytes\n \
				sizeof(Snapshot_Data_t*) = %ld bytes\n \
				sizeof(Snapshot_Data_t) = %ld bytes\n \
				FOR SNAPSHOTS: sizeof(Snapshot_Data_t*) * N_points + sizeof(Snapshot_Data_t) * N_points * N_snapshots = %g Gb\n \
				FOR TEMPORAL MODES: sizeof(double*) * N_snapshots + sizeof(double) * N_snapshots * N_snapshots = %g Gb\n", 
				myrank, N_points, N_snapshots, sizeof(double*), sizeof(double), sizeof(Snapshot_Data_t*), sizeof(Snapshot_Data_t),
				1.0*(sizeof(Snapshot_Data_t*) * N_points + sizeof(Snapshot_Data_t) * N_points * N_snapshots) / 1.0e9,
				1.0*(sizeof(double*) * N_snapshots + sizeof(double) * N_snapshots * N_snapshots) / 1.0e9);

	} else if (spatial_modes_allocated) {
		fprintf(stdout, 
				"MEMORY USAGE on myrank = %d\n \
				N_points = %ld\n \
				N_snapshots = %ld\n \
				sizeof(double*) = %ld bytes\n \
				sizeof(double) = %ld bytes\n \
				sizeof(Snapshot_Data_t*) = %ld bytes\n \
				sizeof(Snapshot_Data_t) = %ld bytes\n \
				FOR SPATIAL MODES: sizeof(Snapshot_Data_t*) * N_points + sizeof(Snapshot_Data_t) * N_points * N_modes = %g Gb\n \
				FOR TEMPORAL MODES: sizeof(double*) * N_snapshots + sizeof(double) * N_snapshots * N_snapshots = %g Gb\n", 
				myrank, N_points, N_snapshots, sizeof(double*), sizeof(double), sizeof(Snapshot_Data_t*), sizeof(Snapshot_Data_t),
				1.0*(sizeof(Snapshot_Data_t*) * N_points + sizeof(Snapshot_Data_t) * N_points * N_modes) / 1.0e9,
				1.0*(sizeof(double*) * N_snapshots + sizeof(double) * N_snapshots * N_snapshots) / 1.0e9);

	} else {
		quit_error("Problem with the memory allocation");
	}
}

//----------------------------------------------------------------------------
void POD_t::report_memory_usage()
{
	report_memory_usage(0);
}

//----------------------------------------------------------------------------
POD_t::~POD_t()
{	
	cout << "Deallocating memory ... " << endl;
	first_snapshot->Delete();
	
	free_1d_array_double(cell_volumes);
	free_1d_array_Snapshot_Data_t(rms_field);
	free_1d_array_Snapshot_Data_t(spatial_rms_history);
	free_1d_array_Snapshot_Data_t(mean_field);
	
	if (snapshots_allocated) {
		free_2d_array_Snapshot_Data_t(snapshots, N_points);
		free_2d_array_double(A, N_snapshots);
		free_1d_array_double(eigenvalues);
		free_1d_array_double(scaled_eigenvalues);

	} else if (spatial_modes_allocated) {
		free_2d_array_Snapshot_Data_t(spatial_modes, N_points);
		free_1d_array_Snapshot_Data_t(max_spatial_L2_error);
		free_2d_array_Snapshot_Data_t(spatial_L2_error, N_snapshots);
	}
}

//----------------------------------------------------------------------------
void POD_t::read_list_of_snapshots(char *snapshot_list_filename, char *snapshot_dir)
{
	ifstream fin(snapshot_list_filename);
	string snapshot_filename; 
	while(getline(fin, snapshot_filename)) {
		snapshot_filename = snapshot_dir + snapshot_filename;
		list_of_snapshots.push_back(snapshot_filename);
	}
}

//----------------------------------------------------------------------------
void POD_t::read_list_of_spatial_modes(char *spatial_mode_dir)
{
	string spatial_mode_list_filename;
	spatial_mode_list_filename = spatial_mode_dir;
	spatial_mode_list_filename = spatial_mode_list_filename + "/spatial_modes.list";
	ifstream fin(spatial_mode_list_filename.c_str());
	string spatial_mode_filename; 
	while(getline(fin, spatial_mode_filename)) {
		spatial_mode_filename = spatial_mode_dir + spatial_mode_filename;
		list_of_spatial_modes.push_back(spatial_mode_filename);
	}
}

//----------------------------------------------------------------------------
void POD_t::read_and_save_first_snapshot(char *file_name_with_cell_volumes)
{
	cout << file_name_with_cell_volumes << endl;
	
	vtkUnstructuredGridReader *first_snapshot_reader = vtkUnstructuredGridReader::New();
	first_snapshot_reader->SetFileName(file_name_with_cell_volumes);
	first_snapshot_reader->ReadAllScalarsOn();	first_snapshot_reader->ReadAllVectorsOn();	first_snapshot_reader->Update();
	cout << "Number of scalars in file = " << first_snapshot_reader->GetNumberOfScalarsInFile() << endl;
	cout << "Number of vectors in file = " << first_snapshot_reader->GetNumberOfVectorsInFile() << endl;
	cout << "Number of cells = " << first_snapshot_reader->GetOutput()->GetNumberOfCells() << endl;
	cout << "Number of points = " << first_snapshot_reader->GetOutput()->GetNumberOfPoints() << endl;
	
/*	vtkExtractUnstructuredGrid *extract = vtkExtractUnstructuredGrid::New();
	extract->SetExtent(x_min, x_max, y_min, y_max, z_min, z_max);
	extract->SetInput(first_snapshot_reader->GetOutput());
	extract->MergingOn();	extract->Update();
	cout << "Number of cells in extracted file = " << extract->GetOutput()->GetNumberOfCells() << endl;
	cout << "Number of points in extracted file = " << extract->GetOutput()->GetNumberOfPoints() << endl;
*/	
	first_snapshot->DeepCopy(first_snapshot_reader->GetOutput());
//	first_snapshot->DeepCopy(extract->GetOutput());
//	extract->Delete();
	first_snapshot_reader->Delete();
}

//----------------------------------------------------------------------------
void POD_t::set_node_number_range(long int first_node_number, long int last_node_number)
{
	this->first_node_number = first_node_number;
	this->last_node_number = last_node_number;
	N_points = last_node_number - first_node_number + 1;
}

//----------------------------------------------------------------------------
long int POD_t::get_total_number_of_points()
{
	double total_N_points = first_snapshot->GetNumberOfPoints();
	return total_N_points;
}

//----------------------------------------------------------------------------
long int POD_t::get_number_of_points()
{
	return N_points;
}

//----------------------------------------------------------------------------
long int POD_t::get_first_point_number()
{
	return first_node_number;
}

//----------------------------------------------------------------------------
long int POD_t::get_last_point_number()
{
	return last_node_number;
}

//----------------------------------------------------------------------------
void POD_t::set_snapshot_number_range(long int first_snapshot_number, long int last_snapshot_number)
{
	this->first_snapshot_number = first_snapshot_number;
	this->last_snapshot_number = last_snapshot_number;
	N_snapshots = last_snapshot_number - first_snapshot_number + 1;
}

//----------------------------------------------------------------------------
long int POD_t::get_number_of_snapshots()
{
	return N_snapshots;	
}

//----------------------------------------------------------------------------
long int POD_t::get_number_of_modes()
{
	return N_modes;
}

//----------------------------------------------------------------------------
long int POD_t::get_last_restart_snapshot_i()
{
	return last_restart_snapshot_i;
}

//----------------------------------------------------------------------------
long int POD_t::get_last_restart_snapshot_j()
{
	return last_restart_snapshot_j;
}

//----------------------------------------------------------------------------
void POD_t::read_snapshots()
{	
	read_snapshots(0);
}

//----------------------------------------------------------------------------
void POD_t::read_snapshots(int myrank)
{		
	vtkUnstructuredGridReader *reader;
	
	for(long int snapshot_num = first_snapshot_number; snapshot_num <= last_snapshot_number; snapshot_num++) {
		fprintf(stdout, "Reading file %ld of %ld name: %s on processor %d for data\n", snapshot_num+1, N_snapshots, list_of_snapshots[snapshot_num].c_str(), myrank);
		reader = vtkUnstructuredGridReader::New();
		reader->SetFileName(list_of_snapshots[snapshot_num].c_str());
		reader->ReadAllScalarsOn();
		reader->ReadAllVectorsOn();
		reader->Update();

/*		vtkExtractUnstructuredGrid *extract = vtkExtractUnstructuredGrid::New();
		extract->SetExtent(x_min, x_max, y_min, y_max, z_min, z_max);
		extract->SetInput(reader->GetOutput());
		extract->MergingOn();	extract->Update();
*/		
		for(long int node = first_node_number; node <= last_node_number; node++) {
/*			snapshots[node-first_node_number][snapshot_num-first_snapshot_number].p = extract->GetOutput()->GetPointData()->GetScalars("p")->GetComponent(node,0);
			snapshots[node-first_node_number][snapshot_num-first_snapshot_number].u = extract->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,0);
			snapshots[node-first_node_number][snapshot_num-first_snapshot_number].v = extract->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,1);
			snapshots[node-first_node_number][snapshot_num-first_snapshot_number].w = extract->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,2);	
*/			snapshots[node-first_node_number][snapshot_num-first_snapshot_number].p = reader->GetOutput()->GetPointData()->GetScalars("p")->GetComponent(node,0);
			snapshots[node-first_node_number][snapshot_num-first_snapshot_number].u = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,0);
			snapshots[node-first_node_number][snapshot_num-first_snapshot_number].v = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,1);
			snapshots[node-first_node_number][snapshot_num-first_snapshot_number].w = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,2);
			if(two_dimensional) snapshots[node-first_node_number][snapshot_num-first_snapshot_number].w = 0.0;
		}
		reader->Delete();
//		extract->Delete();
	}	

	for(long int node = first_node_number; node <= last_node_number; node++)
		cell_volumes[node-first_node_number] = first_snapshot->GetPointData()->GetScalars("cell_volume")->GetComponent(node,0);
}

//----------------------------------------------------------------------------
void POD_t::read_spatial_modes()
{	
	read_spatial_modes(0);
}

//----------------------------------------------------------------------------
void POD_t::read_spatial_modes(int myrank)
{		
	vtkUnstructuredGridReader *reader;
	
	for(long int i = 0; i < N_modes; i++) {
		fprintf(stdout, "Reading file %ld of %ld name: %s on processor %d for data\n", i+1, N_modes, list_of_spatial_modes[i].c_str(), myrank);
		reader = vtkUnstructuredGridReader::New();
		reader->SetFileName(list_of_spatial_modes[i].c_str());
		reader->ReadAllScalarsOn();		
		reader->ReadAllVectorsOn();
		reader->Update();
		
		for(long int node = first_node_number; node <= last_node_number; node++) {
			spatial_modes[node-first_node_number][i].p = reader->GetOutput()->GetPointData()->GetScalars("p")->GetComponent(node,0);
			spatial_modes[node-first_node_number][i].u = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,0);
			spatial_modes[node-first_node_number][i].v = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,1);
			spatial_modes[node-first_node_number][i].w = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,2);
			if(two_dimensional) spatial_modes[node-first_node_number][i].w = 0.0;
		}
		reader->Delete();
	}
	
	for(long int node = first_node_number; node <= last_node_number; node++)
		cell_volumes[node-first_node_number] = first_snapshot->GetPointData()->GetScalars("cell_volume")->GetComponent(node,0);
}

//----------------------------------------------------------------------------
void POD_t::read_mean_field(char *restart_dir)
{		
	read_mean_field(restart_dir, 0);
}

//----------------------------------------------------------------------------
void POD_t::read_mean_field(char *restart_dir, int myrank)
{		
	string mean_field_file_name;
	mean_field_file_name = restart_dir;
	mean_field_file_name = mean_field_file_name + "/mean_field.vtk";
	
	fprintf(stdout, "Reading mean field file: %s on processor %d for data\n", mean_field_file_name.c_str(), myrank);
	vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
	reader->SetFileName(mean_field_file_name.c_str());
	reader->ReadAllScalarsOn();		
	reader->ReadAllVectorsOn();
	reader->Update();
	
/*	vtkExtractUnstructuredGrid *extract = vtkExtractUnstructuredGrid::New();
	extract->SetExtent(x_min, x_max, y_min, y_max, z_min, z_max);
	extract->SetInput(reader->GetOutput());
	extract->MergingOn();	extract->Update();
*/	
	for(long int node = first_node_number; node <= last_node_number; node++) {
/*		mean_field[node-first_node_number].p = extract->GetOutput()->GetPointData()->GetScalars("p")->GetComponent(node,0);
		mean_field[node-first_node_number].u = extract->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,0);
		mean_field[node-first_node_number].v = extract->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,1);
		mean_field[node-first_node_number].w = extract->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,2);		
*/		mean_field[node-first_node_number].p = reader->GetOutput()->GetPointData()->GetScalars("p")->GetComponent(node,0);
		mean_field[node-first_node_number].u = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,0);
		mean_field[node-first_node_number].v = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,1);
		mean_field[node-first_node_number].w = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,2);
		if(two_dimensional) mean_field[node-first_node_number].w = 0.0;
	}
//	extract->Delete();
	reader->Delete();
}


//----------------------------------------------------------------------------
void POD_t::read_rms_field(char *restart_dir)
{		
	read_rms_field(restart_dir, 0);
}

//----------------------------------------------------------------------------
void POD_t::read_rms_field(char *restart_dir, int myrank)
{		
	string rms_field_file_name;
	rms_field_file_name = restart_dir;
	rms_field_file_name = rms_field_file_name + "/rms_field.vtk";
	
	fprintf(stdout, "Reading rms field file: %s on processor %d for data\n", rms_field_file_name.c_str(), myrank);
	vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
	reader->SetFileName(rms_field_file_name.c_str());
	reader->ReadAllScalarsOn();		
	reader->ReadAllVectorsOn();
	reader->Update();
	
	vtkExtractUnstructuredGrid *extract = vtkExtractUnstructuredGrid::New();
	extract->SetExtent(x_min, x_max, y_min, y_max, z_min, z_max);
	extract->SetInput(reader->GetOutput());
	extract->MergingOn();	extract->Update();
	
	for(long int node = first_node_number; node <= last_node_number; node++) {
/*		rms_field[node-first_node_number].p = extract->GetOutput()->GetPointData()->GetScalars("p")->GetComponent(node,0);
		rms_field[node-first_node_number].u = extract->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,0);
		rms_field[node-first_node_number].v = extract->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,1);
		rms_field[node-first_node_number].w = extract->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,2);
*/		rms_field[node-first_node_number].p = reader->GetOutput()->GetPointData()->GetScalars("p")->GetComponent(node,0);
		rms_field[node-first_node_number].u = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,0);
		rms_field[node-first_node_number].v = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,1);
		rms_field[node-first_node_number].w = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,2);
		if(two_dimensional) rms_field[node-first_node_number].w = 0.0;
	}
//	extract->Delete();
	reader->Delete();
}

//----------------------------------------------------------------------------
void POD_t::calculate_temporal_mean()
{		
    double sum_p, sum_u, sum_v, sum_w;
	
	for (int i=0; i<N_points; i++) {
		sum_p = 0.0; sum_u = 0.0; sum_v = 0.0; sum_w = 0.0;
		
		for (int j=0; j<N_snapshots; j++) {
			sum_u = sum_u + snapshots[i][j].u;
			sum_v = sum_v + snapshots[i][j].v;
			sum_w = sum_w + snapshots[i][j].w;
			sum_p = sum_p + snapshots[i][j].p;
		}
		
		mean_field[i].u = sum_u / N_snapshots;
		mean_field[i].v = sum_v / N_snapshots;
		mean_field[i].w = sum_w / N_snapshots;
		mean_field[i].p = sum_p / N_snapshots;		
	}	
}

//----------------------------------------------------------------------------
void POD_t::subtract_temporal_mean_from_snapshots()
{		
	for (int i=0; i<N_points; i++) {
		for (int j=0; j<N_snapshots; j++) {
			snapshots[i][j].u = snapshots[i][j].u - mean_field[i].u;
			snapshots[i][j].v = snapshots[i][j].v - mean_field[i].v;
			snapshots[i][j].w = snapshots[i][j].w - mean_field[i].w;
			snapshots[i][j].p = snapshots[i][j].p - mean_field[i].p;
		}
	}	
}

//----------------------------------------------------------------------------
void POD_t::calculate_temporal_rms()
{		
    double sum_p2, sum_u2, sum_v2, sum_w2;
	for (int i=0; i<N_points; i++) {
		sum_p2 = 0.0; sum_u2 = 0.0; sum_v2 = 0.0; sum_w2 = 0.0;
		for (int j=0; j<N_snapshots; j++) {
			sum_u2 = sum_u2 + SQ(snapshots[i][j].u);
			sum_v2 = sum_v2 + SQ(snapshots[i][j].v);
			sum_w2 = sum_w2 + SQ(snapshots[i][j].w);
			sum_p2 = sum_p2 + SQ(snapshots[i][j].p);
		}
		rms_field[i].u = sqrt(sum_u2 / N_snapshots);
		rms_field[i].v = sqrt(sum_v2 / N_snapshots);
		rms_field[i].w = sqrt(sum_w2 / N_snapshots);
		rms_field[i].p = sqrt(sum_p2 / N_snapshots);
	}
}

//----------------------------------------------------------------------------
void POD_t::calculate_spatial_rms_history()
{
    double sum_p2, sum_u2, sum_v2, sum_w2;
	for (long int j=0; j<N_snapshots; j++) {
		sum_p2 = 0.0; sum_u2 = 0.0; sum_v2 = 0.0; sum_w2 = 0.0;
		for (long int i=0; i<N_points; i++) {
			sum_u2 = sum_u2 + SQ(snapshots[i][j].u) * cell_volumes[i];
			sum_v2 = sum_v2 + SQ(snapshots[i][j].v) * cell_volumes[i];
			sum_w2 = sum_w2 + SQ(snapshots[i][j].w) * cell_volumes[i];
			sum_p2 = sum_p2 + SQ(snapshots[i][j].p) * cell_volumes[i];
		}
		spatial_rms_history[j].u = sqrt(sum_u2);
		spatial_rms_history[j].v = sqrt(sum_v2);
		spatial_rms_history[j].w = sqrt(sum_w2);
		spatial_rms_history[j].p = sqrt(sum_p2);
	}
}

//----------------------------------------------------------------------------
void POD_t::calculate_spatial_average_of_temporal_rms()
{		
	for (int i=0; i<N_snapshots; i++) {
		rms_spatial_average.u = rms_spatial_average.u + SQ(spatial_rms_history[i].u);
		rms_spatial_average.v = rms_spatial_average.v + SQ(spatial_rms_history[i].v);
		rms_spatial_average.w = rms_spatial_average.w + SQ(spatial_rms_history[i].w);
		rms_spatial_average.p = rms_spatial_average.p + SQ(spatial_rms_history[i].p);
	}
	rms_spatial_average.u = sqrt(rms_spatial_average.u / N_snapshots);
	rms_spatial_average.v = sqrt(rms_spatial_average.v / N_snapshots);
	rms_spatial_average.w = sqrt(rms_spatial_average.w / N_snapshots);
	rms_spatial_average.p = sqrt(rms_spatial_average.p / N_snapshots);
}

//----------------------------------------------------------------------------
void POD_t::write_spatial_rms_history()
{
	cout << "Writing spatial rms history... ";	
	
	string output_filename("./results/spatial_rms_history.dat");
	
	FILE *fout = fopen(output_filename.c_str(),"w");
	fprintf(fout,"#\n# snapshot_num, u_rms, v_rms, w_rms, p_rms\n#\n");
	for (long int j=0; j<N_snapshots; j++)
		fprintf(fout,"%ld %16.10e %16.10e %16.10e %16.10e\n", 
				j+1, spatial_rms_history[j].u, spatial_rms_history[j].v, spatial_rms_history[j].w, spatial_rms_history[j].p);
	
	fclose(fout);	

	cout << "done." << endl;
}

//----------------------------------------------------------------------------
void POD_t::read_spatial_rms_history(char *restart_dir)
{
	string spatial_rms_history_filename;
	spatial_rms_history_filename = restart_dir;
	spatial_rms_history_filename = spatial_rms_history_filename + "/spatial_rms_history.dat";
	
	// Open file
	ifstream fin(spatial_rms_history_filename.c_str());
	if (fin.fail()) quit_error((char*)"file %s not found\n", spatial_rms_history_filename.c_str());
	
	// Remove initial comment lines
	char next_char;
	char comment_char;		comment_char = '#';	
	next_char = fin.peek();
	while(next_char == comment_char) {
		fin.ignore(INT_MAX, '\n');
		next_char = fin.peek();
	}
	
	// Get data
	long int snapshot_i_read;
	double u_rms_read, v_rms_read, w_rms_read, p_rms_read;
	while(fin){
		fin >> snapshot_i_read >> u_rms_read >> v_rms_read >> w_rms_read >> p_rms_read;
		spatial_rms_history[snapshot_i_read-1].u = u_rms_read;
		spatial_rms_history[snapshot_i_read-1].v = v_rms_read;
		spatial_rms_history[snapshot_i_read-1].w = w_rms_read;
		spatial_rms_history[snapshot_i_read-1].p = p_rms_read;
	}
	fin.close();
}

//----------------------------------------------------------------------------
void POD_t::calculate_and_write_temporal_correlations()
{	
	string output_filename("./results/temporal_correlations.dat");

	FILE *fout = fopen(output_filename.c_str(),"w");
	fprintf(fout,"#\n# snapshot i, snapshot j, correlation value\n#\n");
	for (long int i=0; i<last_restart_snapshot_i; i++)
		for (long int j=0; j<last_restart_snapshot_j; j++)
			fprintf(fout,"%ld %ld %16.10e\n", i+1, j+1, A[i][j]);
	
	fprintf(stdout,"Calculating temporal correlations... \n");
	for (long int i=last_restart_snapshot_i; i<N_snapshots; i++) {					// Rows of temporal_correlation
		for (long int j=last_restart_snapshot_j; j<N_snapshots; j++) {				// Columns of temporal_correlation

			if (j>=i) {
				A[i][j] = inner_product_calculator.calculate_inner_product(snapshots, snapshots, rms_field, rms_spatial_average, cell_volumes, i, j) / N_snapshots;
			} else {
				A[i][j] = A[j][i];		// Symmetric matrix
			}
			
			if ( (i*N_snapshots+j)%COUNTER_NOTIFY_INTERVAL == 0)
				cout << "Correlation completed for " << i*N_snapshots+j << " of " << N_snapshots*N_snapshots << endl;
		}
		for (long int j=0; j<N_snapshots; j++)
			fprintf(fout,"%ld %ld %16.10e\n", i+1, j+1, A[i][j]);
	}	
	if ( (N_snapshots*N_snapshots)%COUNTER_NOTIFY_INTERVAL != 0)
		cout << "Correlation completed for " << N_snapshots*N_snapshots << " of " << N_snapshots*N_snapshots << endl << endl;
	
	fclose(fout);
}

//----------------------------------------------------------------------------
void POD_t::read_temporal_correlations(char *restart_dir)
{
	string restart_dir_string, temporal_correlations_filename;
	restart_dir_string = restart_dir;
	temporal_correlations_filename = restart_dir_string + "/temporal_correlations.dat";
	
	if (restart_dir_string == "false") {
		last_restart_snapshot_i = 0;
		last_restart_snapshot_j = 0;		
		
	} else {
		fprintf(stdout, "Reading temporal correlations ... \n");
		
		// Open file
		ifstream fin(temporal_correlations_filename.c_str());
		if (fin.fail()) quit_error((char*)"file %s not found\n", temporal_correlations_filename.c_str());
		
		// Remove initial comment lines
		char next_char;
		char comment_char;		comment_char = '#';	
		next_char = fin.peek();
		while(next_char == comment_char) {
			fin.ignore(INT_MAX, '\n');
			next_char = fin.peek();
		}
		
		// Get data
		long int snapshot_i_read, snapshot_j_read;
		double temporal_correlation_read;
		while(fin){
			fin >> snapshot_i_read >> snapshot_j_read >> temporal_correlation_read;
			A[snapshot_i_read-1][snapshot_j_read-1] = temporal_correlation_read;
		}
		fin.close();
		
		last_restart_snapshot_i = snapshot_i_read;
		last_restart_snapshot_j = snapshot_j_read;
	}
}

//----------------------------------------------------------------------------
double POD_t::calculate_inner_product(long int i, long int j)
{
	// Used in MPI version of calculate_temporal_correlations
	double inner_product = inner_product_calculator.calculate_inner_product(snapshots, snapshots, rms_field, rms_spatial_average, cell_volumes, i, j) / N_snapshots;
	return inner_product;
}

//----------------------------------------------------------------------------
void POD_t::calculate_and_scale_eigenvalues(bool check_eigensolution) 
{
	calculate_eigensolution(N_snapshots, A, eigenvalues, check_eigensolution);
	scale_eigenvalues();
}

//----------------------------------------------------------------------------
void POD_t::scale_eigenvalues()
{
	// Find sum of eigenvalues and set those eigenvalues smaller than round off error to zero
	double sum_eigenvalues = 0.0;
	N_non_zero_eigenvalues = 0;
	
	for(int i=0; i<N_snapshots; i++) {
		if(eigenvalues[i]>ZERO) {
			sum_eigenvalues = sum_eigenvalues + eigenvalues[i];
			N_non_zero_eigenvalues++;
		} else {
			eigenvalues[i] = 0.0;
			scaled_eigenvalues[i] = 0.0;
		}
	}
	fprintf(stdout, "num_eigenvalues = %ld, num_non_zero_eigenvalues = %ld\n", N_snapshots, N_non_zero_eigenvalues);
	
	// Now scale eigenvalues to a percentage value
	if (sum_eigenvalues>ZERO) {
		for(int i=0; i<N_non_zero_eigenvalues; i++) scaled_eigenvalues[i] = eigenvalues[i] / sum_eigenvalues * 100.0;
	} else {
		quit_error((char*) "Something wrong with the eigenvalues, they sum to less than round off error\n");
	}
}

//----------------------------------------------------------------------------
void POD_t::calculate_number_of_modes_for_threshhold_energy(double threshhold_energy, int min_number_of_modes, int max_number_of_modes)
{
	double sum_eigenvalues = scaled_eigenvalues[0];
	
	N_modes = 0;
	while ( (sum_eigenvalues<threshhold_energy) && (N_modes<N_non_zero_eigenvalues) && (scaled_eigenvalues[N_modes]>0) ) {
		N_modes++;
		sum_eigenvalues = sum_eigenvalues + scaled_eigenvalues[N_modes];
	}
	fprintf(stdout, "The first %ld modes of %ld contain %f percent of the energy\n", N_modes, N_non_zero_eigenvalues, sum_eigenvalues);
	
	if (N_modes<min_number_of_modes) {
		N_modes = min_number_of_modes;
	} else if (N_modes>max_number_of_modes) {
		N_modes = max_number_of_modes;
	}
	fprintf(stdout, "The first %ld eigenmodes will be written to file.\n", N_modes);
}

//----------------------------------------------------------------------------
void POD_t::calculate_and_scale_temporal_modes()
{
	double norm; 
	for (long int i=0; i<N_non_zero_eigenvalues; i++) {														// For each eigensolution
		norm = 0.0;
		if (eigenvalues[i]<ZERO) {
			for (long int j=0; j<N_snapshots; j++) A[j][i] = 0.0;
		} else {
			for (long int j=0; j<N_snapshots; j++) norm = norm + A[j][i] * A[j][i];							// Square and sum the eigenvector
			norm = norm / eigenvalues[i] / N_snapshots;														// Calculate norm
			for (long int j=0; j<N_snapshots; j++) A[j][i] = A[j][i] / sqrt(norm);							// Scale the eigenvector
		}
	}
}

//----------------------------------------------------------------------------
void POD_t::test_temporal_orthogonality()
{	
	double norm;
	for (long int i=0; i<N_non_zero_eigenvalues; i++) {
		for (long int k=0; k<N_modes; k++) {
			norm = 0.0;
			for (long int j=0; j<N_snapshots; j++) norm = norm + A[j][i] * A[j][k];
			
			if (i==k) {
				if (sqrt(SQ(norm - eigenvalues[i] * N_snapshots)) > ZERO*1.0e3) 
					fprintf(stdout,"WARNING norm check failed for temporal mode pair (%ld, %ld), should be %g, resulted in %g\n", i, k, eigenvalues[i] * N_snapshots, norm);
			} else {
				if (sqrt(SQ(norm)) > ZERO*1.0e3) 
					fprintf(stdout,"WARNING norm check failed for temporal mode pair (%ld, %ld), should be %g, resulted in %g\n", i, k, 0.0, norm);
			}
		}
	}
}

//----------------------------------------------------------------------------
void POD_t::test_spatial_orthogonality() 
{		
	double norm;
	for (long int i=0; i<N_modes; i++) {
		for (long int j=0; j<N_modes; j++) {			
			norm = spatial_mode_calculator.test_spatial_orthogonality(spatial_modes, cell_volumes, rms_field, rms_spatial_average, i, j);
			spatial_mode_calculator.norm_check(i, j, norm);
		}
	}
}

//----------------------------------------------------------------------------
double POD_t::test_spatial_orthogonality(long int i, long int j) 
{		
	double norm = 0.0;
	norm = spatial_mode_calculator.test_spatial_orthogonality(spatial_modes, cell_volumes, rms_field, rms_spatial_average, i, j);
	return norm;
}

//----------------------------------------------------------------------------
void POD_t::calculate_L2_reconstruction_error()
{
	double *inner_products = malloc_1d_array_double(N_modes);
	Snapshot_Data_t **snapshot = malloc_2d_array_Snapshot_Data_t(N_points,1);
	Snapshot_Data_t *pod_reconstruction = malloc_1d_array_Snapshot_Data_t(N_points);
	vtkUnstructuredGridReader *reader;
	double ke;
	
	for (long int snapshot_num=0; snapshot_num<N_snapshots; snapshot_num++) {
		
		// Read snapshot
		fprintf(stdout, "Reading file %ld of %ld name: %s on processor %d for data\n", snapshot_num+1, N_snapshots, list_of_snapshots[snapshot_num].c_str(), 0);
		reader = vtkUnstructuredGridReader::New();
		reader->SetFileName(list_of_snapshots[snapshot_num].c_str());
		reader->ReadAllScalarsOn();		
		reader->ReadAllVectorsOn();
		reader->Update();
		for(long int node = first_node_number; node <= last_node_number; node++) {
			snapshot[node-first_node_number][0].p = reader->GetOutput()->GetPointData()->GetScalars("p")->GetComponent(node,0) 
			- mean_field[node-first_node_number].p;
			snapshot[node-first_node_number][0].u = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,0)
			- mean_field[node-first_node_number].u;
			snapshot[node-first_node_number][0].v = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,1)
			- mean_field[node-first_node_number].v;
			snapshot[node-first_node_number][0].w = reader->GetOutput()->GetPointData()->GetVectors("u")->GetComponent(node,2)
			- mean_field[node-first_node_number].w;
			if(two_dimensional) snapshot[node-first_node_number][0].w = 0.0;
		}
		reader->Delete();
		
		// Reconstruct for each mode
		ke = SQ(spatial_rms_history[snapshot_num].u) + SQ(spatial_rms_history[snapshot_num].v) + SQ(spatial_rms_history[snapshot_num].w);
		for (long int mode_num=0; mode_num<N_modes; mode_num++) {
			inner_products[mode_num] = inner_product_calculator.calculate_inner_product(snapshot, spatial_modes, rms_field, rms_spatial_average, 
																						cell_volumes, 0, mode_num);
			field_reconstructor.reconstruct_field(mode_num, inner_products, spatial_modes, pod_reconstruction);
			field_reconstructor.calculate_spatial_L2_error(mode_num, snapshot_num, snapshot, pod_reconstruction, cell_volumes, spatial_L2_error);
			if (two_dimensional) spatial_L2_error[mode_num][snapshot_num].w = 0.0;
			
			spatial_L2_error[mode_num][snapshot_num].p = spatial_L2_error[mode_num][snapshot_num].p / ke;
			spatial_L2_error[mode_num][snapshot_num].u = spatial_L2_error[mode_num][snapshot_num].u / ke;
			spatial_L2_error[mode_num][snapshot_num].v = spatial_L2_error[mode_num][snapshot_num].v / ke;
			spatial_L2_error[mode_num][snapshot_num].w = spatial_L2_error[mode_num][snapshot_num].w / ke;
			if (two_dimensional) spatial_L2_error[mode_num][snapshot_num].w = 0.0;
		}		
	}
	
	free_2d_array_Snapshot_Data_t(snapshot, N_points);
	free_1d_array_Snapshot_Data_t(pod_reconstruction);
	free_1d_array_double(inner_products);
}

//----------------------------------------------------------------------------
void POD_t::calculate_L_infinite_reconstruction_error()
{
	double total_error, max_total_error;
	long int snapshot_num_with_max_error = 0;
	for (long int mode_num=0; mode_num<N_modes; mode_num++) {
		max_total_error = 0.0;
		for (long int snapshot_num=0; snapshot_num<N_snapshots; snapshot_num++) {
			total_error = spatial_L2_error[mode_num][snapshot_num].u + spatial_L2_error[mode_num][snapshot_num].v + spatial_L2_error[mode_num][snapshot_num].w  + spatial_L2_error[mode_num][snapshot_num].p;
			if (total_error > max_total_error) {
				max_total_error = total_error;
				snapshot_num_with_max_error = snapshot_num;
			}
		}
		max_spatial_L2_error[mode_num].u = spatial_L2_error[mode_num][snapshot_num_with_max_error].u;
		max_spatial_L2_error[mode_num].v = spatial_L2_error[mode_num][snapshot_num_with_max_error].v;
		max_spatial_L2_error[mode_num].w = spatial_L2_error[mode_num][snapshot_num_with_max_error].w;
		max_spatial_L2_error[mode_num].p = spatial_L2_error[mode_num][snapshot_num_with_max_error].p;
	}
}

//----------------------------------------------------------------------------
void POD_t::write_eigenvalues() 
{
	char filename_out[] = "./results/eigenvalues.dat";
	FILE *fout;
	fout=fopen(filename_out,"w");
	
	double sum_eigenvalues = 0;
	fprintf(fout,"#\n# mode, unscaled eigenvalues, scaled eigenvalues, scaled cumulative eigenvalues, scaled singular values\n#\n");
	for (long int mode=0; mode<N_non_zero_eigenvalues; mode++) {
		sum_eigenvalues = sum_eigenvalues + scaled_eigenvalues[mode];
		fprintf(fout,"%ld %16.10e %16.10e %16.10e %16.10e\n", 
				mode+1, eigenvalues[mode], scaled_eigenvalues[mode], sum_eigenvalues, sqrt(eigenvalues[mode]));
	}
	fclose(fout);
}

//----------------------------------------------------------------------------
void POD_t::write_temporal_modes()
{
	cout << "Writing temporal modes... ";
	FILE *fout;									
	
	string temporal_mode_filename;
	string output_location("./results/");
	string name("temporal_mode_");
	string pad_zeros;
	string mode_number_string;
	string extension(".dat");
	
	for (long int mode=0; mode<N_non_zero_eigenvalues; mode++) {
		mode_number_string = convert_int_to_string(mode+1);
		if (mode+1>9999) {
			quit_error((char*)"Too many modes to output, increase digits in output file names\n");
		} else if ( (mode+1<9999) && (mode+1>=1000) ) {
			pad_zeros = "";
		} else if ( (mode+1<1000) && (mode+1>=100) ) {
			pad_zeros = "0";
		} else if ( (mode+1<100) && (mode+1>=10) ) {
			pad_zeros = "00";
		} else if (mode+1<10) {
			pad_zeros = "000";
		}
		temporal_mode_filename = output_location + name + pad_zeros + mode_number_string + extension;
		
		fout=fopen(temporal_mode_filename.c_str(),"w");
		fprintf(fout,"#\n# time, fourier coefficient\n#\n");
		for (long int time_instant=0; time_instant<N_snapshots; time_instant++) fprintf(fout,"%16.10e %16.10e\n", dt*time_instant, A[time_instant][mode]);
		fclose(fout);
	}
	cout << "done." << endl << endl;
}

//----------------------------------------------------------------------------
void POD_t::write_mean_field()
{
	vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
	writer->SetInput(first_snapshot);
	writer->SetFileTypeToBinary();
	
	cout << "Writing mean field... ";	
	for(int node=0; node < first_snapshot->GetNumberOfPoints(); node++) {
		first_snapshot->GetPointData()->GetScalars("p")->SetComponent(node, 0, mean_field[node].p);
		first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 0, mean_field[node].u);
		first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 1, mean_field[node].v);
		first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 2, mean_field[node].w);
	}
	writer->SetFileName("./results/mean_field.vtk");
	writer->Write();
	writer->Delete();	
	cout << "done." << endl << endl;
}

//----------------------------------------------------------------------------
void POD_t::write_first_snapshot_unsteady_component()
{
	vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
	writer->SetInput(first_snapshot);
	writer->SetFileTypeToBinary();
	
	cout << "Writing unsteady component of first snapshot... ";	
	double snapshot_p, snapshot_u, snapshot_v, snapshot_w;
	for(int node=0; node < first_snapshot->GetNumberOfPoints(); node++) {
		snapshot_p = first_snapshot->GetPointData()->GetScalars("p")->GetComponent(node, 0);
		snapshot_u = first_snapshot->GetPointData()->GetVectors("u")->GetComponent(node, 0);
		snapshot_v = first_snapshot->GetPointData()->GetVectors("u")->GetComponent(node, 1);
		snapshot_w = first_snapshot->GetPointData()->GetVectors("u")->GetComponent(node, 2);
		first_snapshot->GetPointData()->GetScalars("p")->SetComponent(node, 0, snapshot_p - mean_field[node].p);
		first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 0, snapshot_u - mean_field[node].u);
		first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 1, snapshot_v - mean_field[node].v);
		first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 2, snapshot_w - mean_field[node].w);
	}
	writer->SetFileName("./results/first_snapshot_unsteady_component.vtk");
	writer->Write();
	writer->Delete();
	cout << "done." << endl << endl;
}


//----------------------------------------------------------------------------
void POD_t::write_rms_field()
{
	vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
	writer->SetInput(first_snapshot);
	writer->SetFileTypeToBinary();
	
	cout << "Writing rms field... ";	
	for(int node=0; node < first_snapshot->GetNumberOfPoints(); node++) {
		first_snapshot->GetPointData()->GetScalars("p")->SetComponent(node, 0, rms_field[node].p);
		first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 0, rms_field[node].u);
		first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 1, rms_field[node].v);
		first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 2, rms_field[node].w);
	}
	writer->SetFileName("./results/rms_field.vtk");
	writer->Write();
	writer->Delete();	
	cout << "done." << endl << endl;
}

//----------------------------------------------------------------------------
void POD_t::write_rms_spatial_average()
{
	cout << "Writing rms spatial average ... ";	
	
	string output_filename("./results/rms_spatial_average.dat");
	
	FILE *fout;
	fout=fopen(output_filename.c_str(),"w");
	
	fprintf(fout,"#\n# u, v, w, p\n#\n");
	fprintf(fout,"%16.10e %16.10e %16.10e %16.10e\n", rms_spatial_average.u, rms_spatial_average.v, rms_spatial_average.w, rms_spatial_average.p);
	fclose(fout);

	cout << "done." << endl;
}

//----------------------------------------------------------------------------
void POD_t::calculate_and_scale_mode_i(long int mode_i) 
{	
	spatial_mode_calculator.calculate_spatial_mode(mode_i, A, snapshots, eigenvalues);
}

//----------------------------------------------------------------------------
void POD_t::calculate_and_write_spatial_modes()
{
	FILE *fout = fopen("./results/spatial_modes.list","w");
	
	string output_filename_string;
	string output_location("./results/");
	string name("spatial_mode_");
	string pad_zeros;
	string mode_number_string;
	string extension(".vtk");
	
	vtkUnstructuredGridWriter *writer;
	
	for (int mode=0; mode<N_modes; mode++) {						// For each mode
		
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
		calculate_and_scale_mode_i(mode); 
		for(int node=0; node < first_snapshot->GetNumberOfPoints(); node++) {
			first_snapshot->GetPointData()->GetScalars("p")->SetComponent(node, 0, spatial_mode_calculator.spatial_mode[node].p);
			first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 0, spatial_mode_calculator.spatial_mode[node].u);
			first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 1, spatial_mode_calculator.spatial_mode[node].v);
			first_snapshot->GetPointData()->GetVectors("u")->SetComponent(node, 2, spatial_mode_calculator.spatial_mode[node].w);
		}
		
		// Write the file
		writer = vtkUnstructuredGridWriter::New();
		writer->SetInput(first_snapshot);
		writer->SetFileTypeToBinary();
		writer->SetFileName(output_filename_string.c_str());
		writer->Write();
		writer->Delete();
		output_filename_string = name + pad_zeros + mode_number_string + extension;
		fprintf(fout,"%s\n",output_filename_string.c_str());
	}
	
	fclose(fout);
}

//----------------------------------------------------------------------------
void POD_t::write_reconstruction_L2_error()
{
	cout << "Writing reconstruction L2 error... ";
	
	string output_filename;
	string output_location("./results/");
	string name("spatial_reconstruction_L2_error_");
	string pad_zeros;
	string mode_number_string;
	string extension(".dat");
	
	FILE *fout;
	double total_error;
	
	for (int mode=0; mode<N_modes; mode++) {						// For each mode
		
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
		output_filename = output_location + name + pad_zeros + mode_number_string + extension;
		fout=fopen(output_filename.c_str(),"w");
		
		// Write the data
		fprintf(fout,"#\n# time, u L2 error, v L2 error, w L2 error, p L2 error, total L2 error\n#\n");
		for (long int snapshot=0; snapshot<N_snapshots; snapshot++) {
			total_error = spatial_L2_error[mode][snapshot].p 
			+ spatial_L2_error[mode][snapshot].u + spatial_L2_error[mode][snapshot].v + spatial_L2_error[mode][snapshot].w;

			fprintf(fout,"%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e\n", 
					(snapshot+1)*dt, 
					spatial_L2_error[mode][snapshot].u,
					spatial_L2_error[mode][snapshot].v,
					spatial_L2_error[mode][snapshot].w,
					spatial_L2_error[mode][snapshot].p,
					total_error);
		}
		
		fclose(fout);
		
	}		
	cout << "done." << endl << endl;
}

//----------------------------------------------------------------------------
void POD_t::write_reconstruction_L_infinite_error()
{
	cout << "Writing max spatial reconstruction L2 error ... ";	

	double total_error;
	string output_filename("./results/max_spatial_reconstruction_L2_error.dat");
	
	FILE *fout;
	fout=fopen(output_filename.c_str(),"w");
	
	fprintf(fout,"#\n# number of modes used in reconstruction, u L2 error, v L2 error, w L2 error, p L2 error, total L2 error\n#\n");
	for (long int mode=0; mode<N_modes; mode++) {
		total_error = max_spatial_L2_error[mode].u + max_spatial_L2_error[mode].v + max_spatial_L2_error[mode].w + max_spatial_L2_error[mode].p;
		
		fprintf(fout,"%ld %16.10e %16.10e %16.10e %16.10e %16.10e\n", 
				mode+1, 
				max_spatial_L2_error[mode].u, 
				max_spatial_L2_error[mode].v, 
				max_spatial_L2_error[mode].w, 
				max_spatial_L2_error[mode].p,
				total_error);
	}
	
	fclose(fout);
	
	cout << "done." << endl << endl;
	
}

//----------------------------------------------------------------------------
void POD_t::write_summary()
{
	fprintf(stdout,"-----------------------------------------------\n");
	fprintf(stdout,"                  Summary\n");
	fprintf(stdout,"-----------------------------------------------\n");
	fprintf(stdout,"Number of nodes used in analysis = %ld\n", N_points);
	fprintf(stdout,"Number of snapshots used in analysis = %ld\n", N_snapshots);
	fprintf(stdout,"Number of modes output = %ld\n", N_modes);		
	fprintf(stdout,"-----------------------------------------------\n\n");
}

//----------------------------------------------------------------------------
