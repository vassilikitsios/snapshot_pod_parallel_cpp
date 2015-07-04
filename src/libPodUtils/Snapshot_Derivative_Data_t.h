
//-------------------------------------------------------------------------
// Snapshot_Derivative_Data_t
//-------------------------------------------------------------------------
class Snapshot_Derivative_Data_t
	{
	private:
		//-------------------------------------------------------------------------
	public:
		//-------------------------------------------------------------------------
		double cell_volume;
		double p;
		double u, du_dx, du_dy, du_dz;
		double v, dv_dx, dv_dy, dv_dz;
		double w, dw_dx, dw_dy, dw_dz;
		double d2u_dx2, d2v_dy2, d2w_dz2;	// to check the surface integrals with test functions
//		double p, dp_dx, dp_dy, dp_dz;
//		double u, du_dx, du_dy, du_dz, d2u_dx2, d2u_dx_dy, d2u_dy2, d2u_dy_dz, d2u_dz2;
//		double v, dv_dx, dv_dy, dv_dz, d2v_dx2, d2v_dx_dy, d2v_dy2, d2v_dy_dz, d2v_dz2;
//		double w, dw_dx, dw_dy, dw_dz, d2w_dx2, d2w_dx_dy, d2w_dy2, d2w_dy_dz, d2w_dz2;
		//-------------------------------------------------------------------------
		// constructor
		//-------------------------------------------------------------------------
		Snapshot_Derivative_Data_t();
		//-------------------------------------------------------------------------
		// applying test functions
		//-------------------------------------------------------------------------		
		void apply_constant_test_field();
		void apply_linear_test_field(double x, double y);
		void apply_quadratic_test_field(double x, double y);
		//-------------------------------------------------------------------------
		// destructor
		//-------------------------------------------------------------------------
		~Snapshot_Derivative_Data_t();
		//-------------------------------------------------------------------------
	};


//-------------------------------------------------------------------------
// Non class functions
//-------------------------------------------------------------------------

Snapshot_Derivative_Data_t **malloc_2d_array_Snapshot_Derivative_Data_t(long int num_points, long int num_modes);
Snapshot_Derivative_Data_t *malloc_1d_array_Snapshot_Derivative_Data_t(long int Nx);
void free_2d_array_Snapshot_Derivative_Data_t(Snapshot_Derivative_Data_t **f, long int Nx);
void free_1d_array_Snapshot_Derivative_Data_t(Snapshot_Derivative_Data_t *f);

