extern int matrix_invert(double *A, double *Ainv, long int N);
extern void matrix_display(double *A, long int R, long int C);
extern void matrix_multiply(double *A, long int Arows, long int Acols, double *B, long int Brows, long int Bcols, double *product);
extern double matrix_determinant(double *A, long int N);
extern void display_augmented(double *A, double *B, long int N);
extern void matrix_transpose(double *A, long int R, long int C);
extern int matrix_cholesky(double *A, double *L, long int N);
extern void matrix_copy(double *A, long int R, long int C, double *Cp);
extern int matrix_cholesky_invert(double *A, double *Ainv, long int N);
extern void matrix_constant_multiply(double *A, long int Arows, long int Acols, double constant, double *product);
extern void matrix_add(double *A, long int Arows, long int Acols, double *B, long int Brows, long int Bcols, double *sum);
extern long int find_null_space_vectors(double *A, long int rows, long int cols, double **P, double absolute_tolerance);
extern void svd(double *A, long int rows, long int cols, double *U, double *D, double *V);
extern int matrix_row_reduce(double *A, long int rows, long int cols, double tolerance);
extern long int matrix_nullity(double *A, long int rows, long int cols, double absolute_tolerance);
extern long int find_null_space_decomposition(double *A, long int rows, long int cols, double **D, double **V, double tolerance);
extern long int find_null_space_decomposition(double *A, long int rows, long int cols, double tolerance);

