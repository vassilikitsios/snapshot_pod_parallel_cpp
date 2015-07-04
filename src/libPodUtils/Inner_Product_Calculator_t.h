
//-------------------------------------------------------------------------
// Inner_Product_Calculator_t
//-------------------------------------------------------------------------
class Inner_Product_Calculator_t
	{
	private:
		//-------------------------------------------------------------------------
		// variables
		//-------------------------------------------------------------------------
		long int N_points;
		char *temporal_correlation_type, *normalisation;
		bool two_dimensional;
		//-------------------------------------------------------------------------
		// Correlations
		//-------------------------------------------------------------------------
		double calculate_inner_product_for_u(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_u_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average,
														double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_v(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_v_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average,
														double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_w(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_w_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average,
														double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_p(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_p_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average,
														double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_ke(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_ke_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average,
														double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_ke_global_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t rms_spatial_average,
														  double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_q(Snapshot_Data_t **a, Snapshot_Data_t **b, double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_q_local_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t *rms_field,  Snapshot_Data_t rms_spatial_average,
														double *cell_volumes, long int i, long int j);
		double calculate_inner_product_for_q_global_norm(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t rms_spatial_average,
														 double *cell_volumes, long int i, long int j);
		//-------------------------------------------------------------------------
	public:
		//-------------------------------------------------------------------------
		// constructor
		//-------------------------------------------------------------------------
		Inner_Product_Calculator_t();
		//-------------------------------------------------------------------------
		// applying test functions
		//-------------------------------------------------------------------------		
		void set_problem_parameters(char *temporal_correlation_type, char *normalisation, bool two_dimensional, long int N_points);
		double calculate_inner_product(Snapshot_Data_t **a, Snapshot_Data_t **b, Snapshot_Data_t *rms_field, Snapshot_Data_t rms_spatial_average,
									   double *cell_volumes, long int i, long int j);
		//-------------------------------------------------------------------------
		// destructor
		//-------------------------------------------------------------------------
		~Inner_Product_Calculator_t();
		//-------------------------------------------------------------------------
	};

