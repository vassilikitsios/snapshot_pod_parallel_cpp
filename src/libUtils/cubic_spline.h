
//****************************************************************************
// Cubic_Spline_t
//****************************************************************************
class Cubic_Spline_t
{
private:	
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	double *a, *b, *c, *d, *points;
	//-------------------------------------------------------------------------	
public:
	//-------------------------------------------------------------------------
	// variables
	//-------------------------------------------------------------------------
	long int num_points;			// Numer of points in the spline
	//-------------------------------------------------------------------------
	// constructor
	//-------------------------------------------------------------------------
	Cubic_Spline_t(long int num_points);
	//-------------------------------------------------------------------------
	// functions
	//-------------------------------------------------------------------------	
	void calculate_interpolation_function(double *points, double *data);
	double evaluate_function(double interp_point);
	double evaluate_first_derivative(double interp_point);
	double evaluate_second_derivative(double interp_point);
	//-------------------------------------------------------------------------
	// destructor
	//-------------------------------------------------------------------------
	~Cubic_Spline_t();
	//-------------------------------------------------------------------------
};
