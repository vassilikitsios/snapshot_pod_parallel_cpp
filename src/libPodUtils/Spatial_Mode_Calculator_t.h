
//-------------------------------------------------------------------------
// Spatial_Mode_Calculator_t
//-------------------------------------------------------------------------
class Spatial_Mode_Calculator_t
	{
	private:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		long int N_points, N_snapshots, N_modes;
		char *temporal_correlation_type, *normalisation;
		//-------------------------------------------------------------------------		
		void calculate_spatial_mode_u(long int mode_num, double **A, Snapshot_Data_t **snapshots);
		void calculate_spatial_mode_v(long int mode_num, double **A, Snapshot_Data_t **snapshots);
		void calculate_spatial_mode_w(long int mode_num, double **A, Snapshot_Data_t **snapshots);
		void calculate_spatial_mode_ke(long int mode_num, double **A, Snapshot_Data_t **snapshots);
		void calculate_spatial_mode_p(long int mode_num, double **A, Snapshot_Data_t **snapshots);
		void calculate_spatial_mode_q(long int mode_num, double **A, Snapshot_Data_t **snapshots);
		//-------------------------------------------------------------------------
		double test_spatial_orthogonality(Snapshot_Data_t **spatial_modes, double *cell_volumes, long int i, long int j);
		double test_spatial_orthogonality_norm_local(Snapshot_Data_t **spatial_modes, double *cell_volumes,
													 Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average, long int i, long int j);
		double test_spatial_orthogonality_norm_global(Snapshot_Data_t **spatial_modes, double *cell_volumes,
													  Snapshot_Data_t rms_spatial_average, long int i, long int j);
		//-------------------------------------------------------------------------
	public:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		Snapshot_Data_t *spatial_mode;
		//-------------------------------------------------------------------------
		// constructor
		//-------------------------------------------------------------------------
		Spatial_Mode_Calculator_t();
		//-------------------------------------------------------------------------
		// applying test functions
		//-------------------------------------------------------------------------		
		void set_problem_parameters(char *temporal_correlation_type, char *normalisation, long int N_points, long int N_snapshots);
		void calculate_spatial_mode(long int mode_num, double **A, Snapshot_Data_t **snapshots, double *eigenvalues);

		double test_spatial_orthogonality(Snapshot_Data_t **spatial_modes, double *cell_volumes, 
										  Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average,
										  long int i, long int j);
		void norm_check(long int i, long int j, double norm);
		//-------------------------------------------------------------------------
		// destructor
		//-------------------------------------------------------------------------
		~Spatial_Mode_Calculator_t();
		//-------------------------------------------------------------------------
	};

