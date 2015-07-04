
//============================================================================
// InputData_t
//============================================================================
class InputData_t
{
private:
	//-------------------------------------------------------------------------
public:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	double x_min;							// Domain of interest
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;

	char *snapshot_list_filename;			// file containing list of snapshot files
	char *snapshot_dir;						// location of snapshots
	char *file_name_with_cell_volumes;		// file number with cell volumes, first file number is 0
	
	char *two_dimensional;					// if data is 2D "true" or "false"
	char *temporal_correlation_type;		// u, v, w, p, ke, q
	char *normalisation;					// none local global
	char *restart;							// 'false' or location of the temporal correlations file
	double dt;								// time between each snapshot - assuming constant dt
	
	bool test_orthogonality;
	bool check_eigensolution;
	
	double threshhold_energy;
	int min_number_of_eigenvalues;
	int max_number_of_eigenvalues;
	
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------
	InputData_t();
	//-------------------------------------------------------------------------
	~InputData_t();
	//-------------------------------------------------------------------------
};

