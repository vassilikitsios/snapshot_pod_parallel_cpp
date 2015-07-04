
//-------------------------------------------------------------------------
// Field_Reconstructor_t
//-------------------------------------------------------------------------
class Field_Reconstructor_t
	{
	private:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		long int N_points;
		char *temporal_correlation_type;
		//-------------------------------------------------------------------------
		void calculate_spatial_L2_error_for_u(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
											  Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error);
		void calculate_spatial_L2_error_for_v(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
											  Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error);
		void calculate_spatial_L2_error_for_w(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
											  Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error);
		void calculate_spatial_L2_error_for_ke(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
											   Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error);
		void calculate_spatial_L2_error_for_p(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
											  Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error);
		void calculate_spatial_L2_error_for_q(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, 
											  Snapshot_Data_t *pod_reconstruction, double *cell_volumes, Snapshot_Data_t **spatial_L2_error);
		//-------------------------------------------------------------------------
	public:
		//-------------------------------------------------------------------------
		// constructor
		//-------------------------------------------------------------------------
		Field_Reconstructor_t();
		//-------------------------------------------------------------------------
		// applying test functions
		//-------------------------------------------------------------------------		
		void set_problem_parameters(char *temporal_correlation_type, long int N_points);
		void reconstruct_field(long int num_modes, double *inner_products, Snapshot_Data_t **spatial_modes, Snapshot_Data_t *pod_reconstruction);
		void calculate_spatial_L2_error(long int num_modes, long int snapshot_num, Snapshot_Data_t **snapshots, Snapshot_Data_t *pod_reconstruction,
										double *cell_volumes, Snapshot_Data_t **spatial_L2_error);
		//-------------------------------------------------------------------------
		// destructor
		//-------------------------------------------------------------------------
		~Field_Reconstructor_t();
		//-------------------------------------------------------------------------
	};

