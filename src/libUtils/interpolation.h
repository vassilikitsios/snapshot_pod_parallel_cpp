
//============================================================================
// Dist_t
//============================================================================
class Dist_t
	{
	public:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		double delta_x;
		double delta_y;
		double delta_z;
		double delta_mag;
		long int node_num;
		//-------------------------------------------------------------------------	
	};

//============================================================================
// Interpolation_t
//============================================================================
class Interpolation_t
{
private:	
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	int order;						// Polynomial Order
	int p;							// Number of terms in polynomial
	int max_redundancy;				// Maximum number of additional points required above that needed to calculate a polynomial surface of order "order"
	int m;							// Number of points used in stencil
	double epsilon;					// Constant used for Gaussian weights
	double singularity_tolerance;	// Tolerance for determining if a stencil is singular or not
	bool two_PI_periodic;			// Flag if mesh is 2PI periodic
	
	long int num_data_points;		// Number of points in the source mesh
	double *D;						// Multiply this matrix by the data points in the stencil to get the results of the interpolation
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------		
	void order_closest_nodes_on_extruded_2D_data(double **x_data, double x_interp_point, double y_interp_point);
	void order_closest_nodes_on_2D_data(double **x_data, double x_interp_point, double y_interp_point);
	
	void generate_stencil_and_test_stability(double x_interp_point, double y_interp_point);
	void generate_stencil(double x_interp_point, double y_interp_point);

	void generate_w_matrix(double *w);
	void generate_w_matrix_for_continuity(double *w);
	
	void generate_B_matrix(double *B);
	void generate_B_matrix_for_continuity(double *B);
	void generate_b_matrix(double *b);
	void generate_b_x_matrix(double *b_x);
	void generate_b_y_matrix(double *b_y);

	void generate_C_matrix(double *interpolation_data, double *C);
	void generate_C_matrix_for_continuity(double *u, double *v, double *C);	
	//-------------------------------------------------------------------------	
public:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	Dist_t *distance_to_points;		// Sorted list of distance to points used to determine the stencil of closest points to the desired interpolation point
	double f, df_dx, df_dy, d2f_dx2, d2f_dxdy, d2f_dy2;			// Results of interpolation
	double u, du_dx, du_dy, d2u_dx2, d2u_dxdy, d2u_dy2;			// Results of interpolation for continuity conservation
	double v, dv_dx, dv_dy, d2v_dx2, d2v_dxdy, d2v_dy2;
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Interpolation_t(long int num_data_points, int order, int max_redundancy, double epsilon, double singularity_tolerance, bool two_PI_periodic);
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------	
	void form_interpolation_function_on_extruded_2D_data(double **x_data, double x_interp_point, double y_interp_point);
	void form_interpolation_function_on_2D_data(double **x_data, double x_interp_point, double y_interp_point);

	void evaluate_interpolation_function(double *interpolation_data);
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Interpolation_t();
	//-------------------------------------------------------------------------
};

//============================================================================
// Non class function
//============================================================================
Dist_t* malloc_Dist_t(long int N);
void free_Dist_t(Dist_t *dist);
int qsort_distance_compare_function(Dist_t *dist1, Dist_t *dist2);
