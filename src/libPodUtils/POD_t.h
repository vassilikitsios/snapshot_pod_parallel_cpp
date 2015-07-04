
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkPointData.h"
#include "vtkExtractUnstructuredGrid.h"

#include <vector>

//-------------------------------------------------------------------------
// POD_t
//-------------------------------------------------------------------------
class POD_t
	{
	private:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		long int N_snapshots, N_points, N_non_zero_eigenvalues, N_modes;
		char *temporal_correlation_type, *normalisation;
		double dt;
		long int first_node_number, last_node_number;
		long int first_snapshot_number, last_snapshot_number;		// for parallelisation - but not used as spatial parallelisaiton is more efficient
		long int last_restart_snapshot_i, last_restart_snapshot_j;	// from reading temporal correlations file
		bool snapshots_allocated, spatial_modes_allocated;
		//-------------------------------------------------------------------------		
		void read_list_of_snapshots(char *snapshot_list_filename, char *snapshot_dir);
		void read_list_of_spatial_modes(char *spatial_modes_dir);
		void read_and_save_first_snapshot(char *first_snapshot_filename);
		//-------------------------------------------------------------------------
		void scale_eigenvalues();
		//-------------------------------------------------------------------------
	public:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		vector<string> list_of_snapshots, list_of_spatial_modes;
		
		vtkUnstructuredGrid *first_snapshot;			// Actually the snapshot with the cell volumes stored, usually the first snapshot but not necessarily
		Snapshot_Data_t **snapshots;
		double *cell_volumes;
		Snapshot_Data_t *mean_field;					// temporal mean of the snapshots
		Snapshot_Data_t *rms_field;						// rms of the unsteady components i.e. the mean has already been subtracted
		Snapshot_Data_t rms_spatial_average;			// spatial average of the rms fields
		Snapshot_Data_t *spatial_rms_history;			// spatial rms of each snapshot
		
		double **A;										// Initially stores the temporal correlation, after eigensolution is complete it stores temporal modes
		double *eigenvalues;
		double *scaled_eigenvalues;
		double x_min, x_max, y_min, y_max, z_min, z_max;	// Domain extents
		Snapshot_Data_t **spatial_modes;
		//-------------------------------------------------------------------------
		bool two_dimensional;
		Inner_Product_Calculator_t inner_product_calculator;
		Spatial_Mode_Calculator_t spatial_mode_calculator;
		Field_Reconstructor_t field_reconstructor;
		//-------------------------------------------------------------------------
		Snapshot_Data_t **spatial_L2_error;
		Snapshot_Data_t *max_spatial_L2_error;	
		//-------------------------------------------------------------------------
		// constructor
		//-------------------------------------------------------------------------
		POD_t(char *temporal_correlation_type, char *normalisation, double dt, char *two_dimensional_char,
			  char *snapshot_list_filename, char *snapshot_dir, char *file_name_with_cell_volumes,
			  double x_min, double x_max, double y_min, double y_max, double z_min, double z_max,
			  char *spatial_modes_dir);
		//-------------------------------------------------------------------------
		// applying test functions
		//-------------------------------------------------------------------------
		void allocate_memory();
		void report_memory_usage(int myrank);
		void report_memory_usage();
		
		void set_node_number_range(long int first_node_number, long int last_node_number);
		long int get_number_of_points();
		long int get_first_point_number();
		long int get_last_point_number();
		long int get_total_number_of_points();

		void set_snapshot_number_range(long int first_snapshot_number, long int last_snapshot_number);
		long int get_number_of_snapshots();
		long int get_number_of_modes();

		void read_snapshots();
		void read_snapshots(int myrank);
		
		void read_spatial_modes();
		void read_spatial_modes(int myrank);
		
		void read_temporal_correlations(char *temporal_correlations_filename);
		long int get_last_restart_snapshot_i();
		long int get_last_restart_snapshot_j();
		
		void calculate_temporal_mean();
		void write_mean_field();
		void read_mean_field(char *restart_dir, int myrank);
		void read_mean_field(char *restart_dir);
		
		void subtract_temporal_mean_from_snapshots();
		void write_first_snapshot_unsteady_component();
		
		void calculate_temporal_rms();
		void write_rms_field();
		void read_rms_field(char *restart_dir, int myrank);
		void read_rms_field(char *restart_dir);
		
		void calculate_spatial_rms_history();
		void write_spatial_rms_history();
		void read_spatial_rms_history(char *restart_dir);
		
		void calculate_spatial_average_of_temporal_rms();
		void write_rms_spatial_average();
		
		void calculate_and_write_temporal_correlations();
		double calculate_inner_product(long int i, long int j);				// Used for MPI driver
		
		void calculate_and_scale_eigenvalues(bool check_eigensolution);
		void calculate_number_of_modes_for_threshhold_energy(double threshhold_energy, int min_number_of_eigenvalues, int max_number_of_eigenvalues);
		void write_eigenvalues();
		
		void calculate_and_scale_temporal_modes();
		void test_temporal_orthogonality();
		void write_temporal_modes();
		void read_temporal_pod_modes(char *file_location);
		
		void calculate_and_scale_mode_i(long int mode_i);
		void calculate_and_write_spatial_modes();
		
		void test_spatial_orthogonality();
		double test_spatial_orthogonality(long int i, long int j);
		
		void calculate_L2_reconstruction_error();
		void write_reconstruction_L2_error();
		
		void calculate_L_infinite_reconstruction_error();
		void write_reconstruction_L_infinite_error();
		
		void write_summary();
		//-------------------------------------------------------------------------
		// destructor
		//-------------------------------------------------------------------------
		~POD_t();
		//-------------------------------------------------------------------------
	};

