//========================
//memory_allocations.h
//========================

//****************************************************************************
// memory allocations
//****************************************************************************
// int
extern int    *malloc_1d_array_int(long int Nx);
extern int   **malloc_2d_array_int(long int Nx, long int Ny);
extern int  ***malloc_3d_array_int(long int Nx, long int Ny, long int Nz);
extern int ****malloc_4d_array_int(int num_dimensions, long int Nx, long int Ny, long int Nz);
// int - fortran indices
extern int    *malloc_1d_array_int_fortran_indices(long int i0, long int i1);
extern int   **malloc_2d_array_int_fortran_indices(long int i0, long int i1, long int j0, long int j1);
extern int  ***malloc_3d_array_int_fortran_indices(long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
extern int ****malloc_4d_array_int_fortran_indices(long int I0, long int I1, long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
// long int
extern long int    *malloc_1d_array_long_int(long int Nx);
extern long int   **malloc_2d_array_long_int(long int Nx, long int Ny);
// bool
extern bool *malloc_1d_array_bool(long int Nx);
// double
extern double    *malloc_1d_array_double(long int Nx);
extern double   **malloc_2d_array_double(long int Nx, long int Ny);
extern double  ***malloc_3d_array_double(long int Nx, long int Ny, long int Nz);
extern double ****malloc_4d_array_double(int num_dimensions, long int Nx, long int Ny, long int Nz);
// double - fortran indices
extern double    *malloc_1d_array_double_fortran_indices(long int i0, long int i1);
extern double   **malloc_2d_array_double_fortran_indices(long int i0, long int i1, long int j0, long int j1);
extern double  ***malloc_3d_array_double_fortran_indices(long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
extern double ****malloc_4d_array_double_fortran_indices(long int I0, long int I1, long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
// char
extern char *malloc_1d_array_char(long int Nx);
extern char **malloc_2d_array_char(long int Nx, long int Ny);

//****************************************************************************
// memory re-allocations
//****************************************************************************
// double
extern double *realloc_1d_array_double(double *f_old, long int Nx);

//****************************************************************************
// memory de-allocations
//****************************************************************************
// double
extern void free_1d_array_double(double *f);
extern void free_2d_array_double(double **f, long int Nx);
extern void free_3d_array_double(double ***f, long int Nx, long int Ny);
extern void free_4d_array_double(double ****f, int num_dimensions, long int Nx, long int Ny);
// double - fortran indices
extern void free_1d_array_double_fortran_indices(double *f,long int i0, long int i1);
extern void free_2d_array_double_fortran_indices(double **f,long int i0, long int i1, long int j0, long int j1);
extern void free_3d_array_double_fortran_indices(double ***f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
extern void free_4d_array_double_fortran_indices(double ****f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1, long int l0, long int l1);
// int
extern void free_1d_array_int(int *f);
extern void free_2d_array_int(int **f, long int Nx);
extern void free_3d_array_int(int ***f, long int Nx, long int Ny);
extern void free_4d_array_int(int ****f, int num_dimensions, long int Nx, long int Ny);
// int - fortran indices
extern void free_1d_array_int_fortran_indices(int *f,long int i0, long int i1);
extern void free_2d_array_int_fortran_indices(int **f,long int i0, long int i1, long int j0, long int j1);
extern void free_3d_array_int_fortran_indices(int ***f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1);
extern void free_4d_array_int_fortran_indices(int ****f,long int i0, long int i1, long int j0, long int j1, long int k0, long int k1, long int l0, long int l1);
// long int
extern void free_1d_array_long_int(long int *f);
extern void free_2d_array_long_int(long int **f, long int Nx);
// bool
extern void free_1d_array_bool(bool *f);
//char
extern void free_1d_array_char(char *f);
extern void free_2d_array_char(char **f, long int Nx);

