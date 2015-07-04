
//-------------------------------------------------------------------------
// Snapshot_Data_t
//-------------------------------------------------------------------------
class Snapshot_Data_t
	{
	private:
		//-------------------------------------------------------------------------
	public:
		//-------------------------------------------------------------------------
		double u, v, w, p;
	};

//-------------------------------------------------------------------------
// Non class functions
//-------------------------------------------------------------------------
Snapshot_Data_t **malloc_2d_array_Snapshot_Data_t(long int num_points, long int num_modes);
Snapshot_Data_t *malloc_1d_array_Snapshot_Data_t(long int Nx);
void free_2d_array_Snapshot_Data_t(Snapshot_Data_t **f, long int Nx);
void free_1d_array_Snapshot_Data_t(Snapshot_Data_t *f);
