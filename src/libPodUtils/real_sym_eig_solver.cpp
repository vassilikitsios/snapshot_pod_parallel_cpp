#include <iostream>
#include <math.h>

using namespace std;

#include "Utils.h"
#include "real_sym_eig_solver.h"

//----------------------------------------------------------------------------
extern "C" int dsyevr_(char *jobz, char *range, char *uplo, long int *n,
					   double *a, long int *lda, double *vl, double *vu, long int *
					   il, long int *iu, double *abstol, long int *m, double *w, 
					   double *z__, long int *ldz, long int *isuppz, double *work, 
					   long int *lwork, long int *iwork, long int *liwork, long int *info);

/* 
DSYEVR computes selected eigenvalues and, optionally, eigenvectors   
of a real symmetric matrix T.  Eigenvalues and eigenvectors can be   
selected by specifying either a range of values or a range of   
indices for the desired eigenvalues.   

Whenever possible, DSYEVR calls DSTEGR to compute the   
eigenspectrum using Relatively Robust Representations.  DSTEGR   
computes eigenvalues by the dqds algorithm, while orthogonal   
eigenvectors are computed from various "good" L D L^T representations   
(also known as Relatively Robust Representations). Gram-Schmidt   
orthogonalization is avoided as far as possible. More specifically,   
the various steps of the algorithm are as follows. For the i-th   
unreduced block of T,   
(a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T   
is a relatively robust representation,   
(b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high   
relative accuracy by the dqds algorithm,   
(c) If there is a cluster of close eigenvalues, "choose" sigma_i   
close to the cluster, and go to step (a),   
(d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,   
compute the corresponding eigenvector by forming a   
rank-revealing twisted factorization.   
The desired accuracy of the output can be specified by the input   
parameter ABSTOL.   

For more details, see "A new O(n^2) algorithm for the symmetric   
tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,   
Computer Science Division Technical Report No. UCB//CSD-97-971,   
UC Berkeley, May 1997.   


Note 1 : DSYEVR calls DSTEGR when the full spectrum is requested   
on machines which conform to the ieee-754 floating point standard.   
DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and   
when partial spectrum requests are made.   

Normal execution of DSTEGR may create NaNs and infinities and   
hence may abort due to a floating point exception in environments   
which do not handle NaNs and infinities in the ieee standard default   
manner.   

Arguments   
=========   

JOBZ    (input) CHARACTER*1   
= 'N':  Compute eigenvalues only;   
= 'V':  Compute eigenvalues and eigenvectors.   

RANGE   (input) CHARACTER*1   
= 'A': all eigenvalues will be found.   
= 'V': all eigenvalues in the half-open interval (VL,VU]   
will be found.   
= 'I': the IL-th through IU-th eigenvalues will be found.   
********* For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and   
********* DSTEIN are called   

UPLO    (input) CHARACTER*1   
= 'U':  Upper triangle of A is stored;   
= 'L':  Lower triangle of A is stored.   

N       (input) INTEGER   
The order of the matrix A.  N >= 0.   

A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
On entry, the symmetric matrix A.  If UPLO = 'U', the   
leading N-by-N upper triangular part of A contains the   
upper triangular part of the matrix A.  If UPLO = 'L',   
the leading N-by-N lower triangular part of A contains   
the lower triangular part of the matrix A.   
On exit, the lower triangle (if UPLO='L') or the upper   
triangle (if UPLO='U') of A, including the diagonal, is   
destroyed.   

LDA     (input) INTEGER   
The leading dimension of the array A.  LDA >= max(1,N).   

VL      (input) DOUBLE PRECISION   
VU      (input) DOUBLE PRECISION   
If RANGE='V', the lower and upper bounds of the interval to   
be searched for eigenvalues. VL < VU.   
Not referenced if RANGE = 'A' or 'I'.   

IL      (input) INTEGER   
IU      (input) INTEGER   
If RANGE='I', the indices (in ascending order) of the   
smallest and largest eigenvalues to be returned.   
1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.   
Not referenced if RANGE = 'A' or 'V'.   

ABSTOL  (input) DOUBLE PRECISION   
The absolute error tolerance for the eigenvalues.   
An approximate eigenvalue is accepted as converged   
when it is determined to lie in an interval [a,b]   
of width less than or equal to   

ABSTOL + EPS *   max( |a|,|b| ) ,   

where EPS is the machine precision.  If ABSTOL is less than   
or equal to zero, then  EPS*|T|  will be used in its place,   
where |T| is the 1-norm of the tridiagonal matrix obtained   
by reducing A to tridiagonal form.   

See "Computing Small Singular Values of Bidiagonal Matrices   
with Guaranteed High Relative Accuracy," by Demmel and   
Kahan, LAPACK Working Note #3.   

If high relative accuracy is important, set ABSTOL to   
DLAMCH( 'Safe minimum' ).  Doing so will guarantee that   
eigenvalues are computed to high relative accuracy when   
possible in future releases.  The current code does not   
make any guarantees about high relative accuracy, but   
furutre releases will. See J. Barlow and J. Demmel,   
"Computing Accurate Eigensystems of Scaled Diagonally   
Dominant Matrices", LAPACK Working Note #7, for a discussion   
of which matrices define their eigenvalues to high relative   
accuracy.   

M       (output) INTEGER   
The total number of eigenvalues found.  0 <= M <= N.   
If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.   

W       (output) DOUBLE PRECISION array, dimension (N)   
The first M elements contain the selected eigenvalues in   
ascending order.   

Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))   
If JOBZ = 'V', then if INFO = 0, the first M columns of Z   
contain the orthonormal eigenvectors of the matrix A   
corresponding to the selected eigenvalues, with the i-th   
column of Z holding the eigenvector associated with W(i).   
If JOBZ = 'N', then Z is not referenced.   
Note: the user must ensure that at least max(1,M) columns are   
supplied in the array Z; if RANGE = 'V', the exact value of M   
is not known in advance and an upper bound must be used.   

LDZ     (input) INTEGER   
The leading dimension of the array Z.  LDZ >= 1, and if   
JOBZ = 'V', LDZ >= max(1,N).   

ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )   
The support of the eigenvectors in Z, i.e., the indices   
indicating the nonzero elements in Z. The i-th eigenvector   
is nonzero only in elements ISUPPZ( 2*i-1 ) through   
ISUPPZ( 2*i ).   
********* Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1   

WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

LWORK   (input) INTEGER   
The dimension of the array WORK.  LWORK >= max(1,26*N).   
For optimal efficiency, LWORK >= (NB+6)*N,   
where NB is the max of the blocksize for DSYTRD and DORMTR   
returned by ILAENV.   

If LWORK = -1, then a workspace query is assumed; the routine   
only calculates the optimal size of the WORK array, returns   
this value as the first entry of the WORK array, and no error   
message related to LWORK is issued by XERBLA.   

IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.   

LIWORK  (input) INTEGER   
The dimension of the array IWORK.  LIWORK >= max(1,10*N).   

If LIWORK = -1, then a workspace query is assumed; the   
routine only calculates the optimal size of the IWORK array,   
returns this value as the first entry of the IWORK array, and   
no error message related to LIWORK is issued by XERBLA.   

INFO    (output) INTEGER   
= 0:  successful exit   
< 0:  if INFO = -i, the i-th argument had an illegal value   
> 0:  Internal error   

Further Details   
===============   

Based on contributions by   
Inderjit Dhillon, IBM Almaden, USA   
Osni Marques, LBNL/NERSC, USA   
Ken Stanley, Computer Science Division, University of   
California at Berkeley, USA   

=====================================================================*/


//----------------------------------------------------------------------------
// helper function declaration
//----------------------------------------------------------------------------
void copy_matrix_A_to_vector_a(long int n, double **A, double *a);

void evaluate_error_of_all_eigensolutions(long int n, long int num_eigenvalues, double **A, 
										  double *reverse_order_eigenvectors, double *reverse_order_eigenvalues);
double evaluate_eigensolution_error(long int n, double **A, double *eigenvector, double eigenvalue);

void copy_and_swap_order_of_solution_output(long int n, long int num_eigenvalues, 
											double **A, double *reverse_order_eigenvectors, double *eigenvalues, double *reverse_order_eigenvalues);


//----------------------------------------------------------------------------
// function code
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
long int calculate_real_sym_eigensolution(long int n, double **A, double *eigenvalues)
{
	long int num_eigenvalues = calculate_eigensolution(n, A, eigenvalues, true);
	return num_eigenvalues;
}

//----------------------------------------------------------------------------
long int calculate_eigensolution(long int n, double **A, double *eigenvalues)
{
	long int num_eigenvalues = calculate_eigensolution(n, A, eigenvalues, true);
	return num_eigenvalues;
}

//----------------------------------------------------------------------------
long int calculate_eigensolution(long int n, double **A, double *eigenvalues, bool check_solution)
{
	cout << "Calculating eigenvalues... " << endl;
	
	
	// Set options for clapack routine and vector and matrices sizes
	char jobz = 'V';							// Find eigenvalues and eigenvectors
	char range = 'A';							// Find all eigensolutions
	char uplo = 'U';							// Upper triangular stored in 'a'

	long int lda = n;							// Leading dimension of matrix 'a'
	long int ldz = n;							// Leading dimension of matrix 'z'
	
	long int info = 1;						// Flag set to 0 if calculation successful
	long int num_eigenvalues = 0;				// Total number of eigenvalues found
	
	double vl = 0.0;							// Not referenced because range = 'A';
	double vu = 0.0;							// Not referenced because range = 'A';
	long int il = 0;							// Not referenced because range = 'A';
	long int iu = 0;							// Not referenced because range = 'A';
	
	double abstol = 1.0e-15;					// Error tolerance

	long int lwork, liwork;
	

	// Memory allocation
	long int *iwork, *isuppz;
	double *work, *a;
	double *reverse_order_eigenvalues, *reverse_order_eigenvectors;

	work = malloc_1d_array_double(1);		work[0] = 0.0;
	iwork = malloc_1d_array_long_int(1);	iwork[0] = 0;
	
	reverse_order_eigenvectors = malloc_1d_array_double(n*n);
	for (int i=0; i<n*n; i++) reverse_order_eigenvectors[i] = 0.0; 
	
	reverse_order_eigenvalues = malloc_1d_array_double(n);
	for (int i=0; i<n; i++) reverse_order_eigenvalues[i] = 0.0; 
	
	a = malloc_1d_array_double(n*lda);
	for (int i=0; i<n*lda; i++) a[i] = 0.0; 
	
	isuppz = malloc_1d_array_long_int(2*n);
	for (int i=0; i<2*n; i++) isuppz[i] = 0; 
	
		
	// Perform eigenvalue calcualtion
	copy_matrix_A_to_vector_a(n, A, a);
	
	lwork = 26*n;								// Initial guess of size of vector 'work'
	liwork = 10*n;								// Initial guess of size of vector 'iwork'
	cout << "Initial estimate of work array sizes: lwork = " << lwork << " liwork = " << liwork << endl;	
	lwork = -1;									// Flag to query the size of vector 'work'
	liwork = -1;								// Flag to query the size of vector 'iwork'
	dsyevr_(&jobz, &range, &uplo, &n,
			a, &lda, &vl, &vu, &il, &iu, &abstol, &num_eigenvalues, reverse_order_eigenvalues, 
			reverse_order_eigenvectors, &ldz, isuppz, work, 
			&lwork, iwork, &liwork, &info);		// Query optimal size of 'work' and 'work'
	
	lwork = (long int) work[0];							// Optimal size of vector 'work'
	liwork = (long int) iwork[0];							// Optimal size of vector 'iwork'				
	free_1d_array_double(work);
	free_1d_array_long_int(iwork);	
	cout << "Optimal size of work arrays: lwork = " << lwork << " liwork = " << liwork << endl;
		
	work = malloc_1d_array_double(lwork);
	for (int i=0; i<lwork; i++) work[i] = 0.0;

	iwork = malloc_1d_array_long_int(liwork);
	for (int i=0; i<liwork; i++) iwork[i] = 0;
	
	dsyevr_(&jobz, &range, &uplo, &n,
			a, &lda, &vl, &vu, &il, &iu, &abstol, &num_eigenvalues, reverse_order_eigenvalues, 
			reverse_order_eigenvectors, &ldz, isuppz, work, 
			&lwork, iwork, &liwork, &info);	
	cout << num_eigenvalues << " eigenvalues found out of a possible " << n << endl << endl;		
	
	// Checking eigenvalue solution
	if(check_solution) evaluate_error_of_all_eigensolutions(n, num_eigenvalues, A, reverse_order_eigenvectors, reverse_order_eigenvalues);
	
	// Swap order of solution output
	copy_and_swap_order_of_solution_output(n, num_eigenvalues, A, reverse_order_eigenvectors, eigenvalues, reverse_order_eigenvalues);
		
	// Free allocated memory
	free_1d_array_double(reverse_order_eigenvectors);
	free_1d_array_double(reverse_order_eigenvalues);
	free_1d_array_double(a);
	free_1d_array_double(work);
	free_1d_array_long_int(iwork);	
	free_1d_array_long_int(isuppz);
		
	return num_eigenvalues;
}

//----------------------------------------------------------------------------
void copy_matrix_A_to_vector_a(long int n, double **A, double *a) {
	int pos=0;
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			a[pos] = A[j][i];					//Copy matrix 'A' across to vector 'a'
			pos++;
		}
	}
}

//----------------------------------------------------------------------------
void evaluate_error_of_all_eigensolutions(long int n, long int num_eigenvalues, double **A, 
										  double *reverse_order_eigenvectors, double *reverse_order_eigenvalues)
{
	double *test_eigenvector;
	double error;

	cout << "Checking eigensolutions " << endl;
	test_eigenvector = malloc_1d_array_double(n);
	
	for (long int i=num_eigenvalues-1; i>=0; i--) {
		for (long int j=0; j<n; j++) test_eigenvector[j] = reverse_order_eigenvectors[i*num_eigenvalues+j];
		error = evaluate_eigensolution_error(n,A,test_eigenvector,reverse_order_eigenvalues[i]);
		if (error>ZERO) cout << "Root mean squared error for eigensolution " << num_eigenvalues - i << " = " << error << endl;		
	}
	cout << "Finished checking eigensolutions " << endl << endl;	
 	free_1d_array_double(test_eigenvector);	
}

//----------------------------------------------------------------------------
double evaluate_eigensolution_error(long int n, double **A, double *eigenvector, double eigenvalue)
{
	double product_row_of_A_and_eigenvector, error=0.0;
	
	for (long int i=0; i<n; i++) {
		product_row_of_A_and_eigenvector = 0.0;
		for (long int j=0; j<n; j++) product_row_of_A_and_eigenvector = product_row_of_A_and_eigenvector + A[i][j] * eigenvector[j];
		error = error + pow(product_row_of_A_and_eigenvector - eigenvalue * eigenvector[i],2.0);
	}
	error = sqrt(error/n);
	return error;
}

//----------------------------------------------------------------------------
void copy_and_swap_order_of_solution_output(long int n, long int num_eigenvalues, 
											double **A, double *reverse_order_eigenvectors, double *eigenvalues, double *reverse_order_eigenvalues)
{
	// Copy eigenvectors over to matrix 'A' swap so that the largest eigenvalues are first
	for (long int i=num_eigenvalues-1; i>=0; i--) {
		for (long int j=0; j<n; j++) {
			A[j][num_eigenvalues-i-1] = reverse_order_eigenvectors[i*num_eigenvalues+j];
		}
	}
	
	// Swap so that the largest eigenvalues are first
	for (long int i=num_eigenvalues-1; i>=0; i--) eigenvalues[num_eigenvalues-i-1] = reverse_order_eigenvalues[i];	
}

//----------------------------------------------------------------------------

