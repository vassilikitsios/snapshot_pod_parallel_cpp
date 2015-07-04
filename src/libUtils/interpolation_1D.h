
//============================================================================
// Interpolation_1D_t
//============================================================================
class Interpolation_1D_t
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
	
	long int num_data_points;		// Number of points in the source mesh
	double *D;						// Multiply this matrix by the data points in the stencil to get the results of the interpolation
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------
	void generate_stencil(double x_interp_point);
	void order_closest_points(double *x_data, double x_interp_point);
	void generate_w_matrix(double *w);
	void generate_B_matrix(double *B);
	void generate_C_matrix(double *interpolation_data, double *C);
	//-------------------------------------------------------------------------	
public:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	Dist_t *distance_to_points;			// Sorted list of distance to points used to determine the stencil of closest points to the desired interpolation point
	double f, df_dx, d2f_dx2;			// Results of interpolation
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Interpolation_1D_t(long int num_data_points, int order, int max_redundancy, double epsilon, double singularity_tolerance);
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------	
	void form_interpolation_function(double *x_data, double x_interp_point);
	void evaluate_interpolation_function(double *interpolation_data);
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Interpolation_1D_t();
	//-------------------------------------------------------------------------
};
