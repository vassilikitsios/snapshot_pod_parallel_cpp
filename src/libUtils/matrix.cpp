#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>

#include "matrix.h"
#include "misc.h"
#include "memory_allocations.h"

// ****************************************************************************
// Swaps rows row1 and row2 of both matricies A and B, each of size N x N
// ****************************************************************************
void swap_row(long int row1, long int row2, double *A, double *B, long int N)
{
	if (A==NULL || B==NULL || N<1 || row1<0 || row1>=N || row2<0 || row2>=N) quit_error((char*)"Bad input data to swap_row");
	double temp;
	for (long int col=0; col<N; col++) {
		temp=A[row1*N+col];
		A[row1*N+col]=A[row2*N+col];
		A[row2*N+col]=temp;
		
		temp=B[row1*N+col];
		B[row1*N+col]=B[row2*N+col];
		B[row2*N+col]=temp;
	}
}

// ****************************************************************************
// Multplies one row of both matricies A and B (each of size N x N) 
// by a certain factor
// ****************************************************************************
void multiply_row(long int row, double factor, double *A, double *B, long int N)
{
	if (A==NULL || B==NULL || N<1 || row<0 || row>=N || fabs(factor)<DBL_EPSILON) quit_error((char*)"Bad input data to multiply_row");
	for (long int col=0; col<N; col++) {
		A[row*N+col]*=factor;
		B[row*N+col]*=factor;		
	}
}

// ****************************************************************************
// Adds factor times row srow to row drow, in both matricies A and B 
// (each of size N x N)
// ****************************************************************************
void add_multiple_row(double factor, long int srow, long int drow, double *A, double *B, long int N)
{
	if (A==NULL || B==NULL || N<1 || srow<0 || srow>=N || drow<0 || drow>=N) quit_error((char*)"Bad input data to add_multiple_row");
	for (long int col=0; col<N; col++) {
		A[drow*N+col]+=factor*A[srow*N+col];
		B[drow*N+col]+=factor*B[srow*N+col];		
	}
}

// ****************************************************************************
// Displays the augmented matrix A|B, where A and B are N x N matricies.
// ****************************************************************************
void display_augmented(double *A, double *B, long int N)
{
	for (long int row=0; row<N; row++) {
		for (long int col=0; col<N; col++) {
			printf("%10f", A[row*N+col]);
			if (col<N-1) printf(", ");
		}
		printf(" | ");
		for (long int col=0; col<N; col++) {
			printf("%10f", B[row*N+col]);
			if (col<N-1) printf(", ");
		}
		printf("\n");
	}
	printf("\n");
}

// ****************************************************************************
// Calculates the inverse of N x N matrix A and stores it in matrix Ainv. 
// It is assumed that the memory required for this has already been allocated. 
// The data in matrix A is not destroyed, unless the same address is supplied for Ainv.
// (This case should be successfully handled also.)
// The function returns -1 if the matrix is not invertible.
// ****************************************************************************
int matrix_invert(double *A, double *Ainv, long int N)
{
	double *backup = malloc_1d_array_double(N*N);
	double pivot, factor;
	long int pivot_row;

	if (A==NULL || Ainv==NULL || N<1) quit_error((char*)"Bad input data to matrix_invert");
	
	// Backup the A matrix so that manipulations can be performed here
	// without damaging the original matrix.
	for (long int row=0; row<N; row++)
		for (long int col=0; col<N; col++)
			backup[row*N+col]=A[row*N+col];
	
	// First, fill Ainv with the identity matrix
	for (long int row=0; row<N; row++) {
		for (long int col=0; col<N; col++) {
			if (row==col) Ainv[row*N+col]=1;
			else Ainv[row*N+col]=0;
		}
	}

	// Calculate the inverse using Gauss-Jordan elimination
	for (long int col=0; col<N; col++) {
		//display_augmented(backup, Ainv, N);
		pivot_row=-1;
		pivot=0.0;
		for (long int row=col; row<N; row++) {
			if (fabs(backup[row*N+col])>fabs(pivot)) {
				pivot_row=row;
				pivot=backup[row*N+col];
			}
		}
		//printf("Pivot: %lf\n\n", pivot);

		if (pivot_row<0 || pivot_row>=N || fabs(pivot)<DBL_EPSILON) {
			free_1d_array_double(backup);
			return -1;
		} else {
			swap_row(col, pivot_row, backup, Ainv, N);
			multiply_row(col, 1.0/pivot, backup, Ainv, N);
			for (long int elim_row=0; elim_row<N; elim_row++) {
				if (elim_row!=col) {
					factor=-backup[elim_row*N+col];
					add_multiple_row(factor, col, elim_row, backup, Ainv, N);
				}
			}
		}
	}
	free_1d_array_double(backup);
	
	return 0;
}

// ****************************************************************************
// Displays the R x C matrix A
// ****************************************************************************
void matrix_display(double *A, long int R, long int C)
{
	for (long int rows=0; rows<R; rows++) {
		for (long int cols=0; cols<C; cols++) {
			printf("%10f", A[rows*C+cols]);
			if (cols!=C-1) printf(", ");
		}
		printf("\n");
	}
}

// ****************************************************************************
// Multiplies matrix A (Arows x Acols) by matrix B (Brows x Bcols) and writes the result to product. 
// It is assumed that the required amount of memory will already be allocated.
// If the matricies are not of the correct dimensions for multiplication then an error is thrown. 
// product is not allowed to point to the same area of memory as A or B
// ****************************************************************************
void matrix_multiply(double *A, long int Arows, long int Acols, double *B, long int Brows, long int Bcols, double *product)
{
	if (A==NULL || B==NULL || product==NULL) quit_error((char*)"Bad input data to matrix_multiply");

	if (Acols!=Brows) quit_error((char*)"Matricies are incompatible for multiplication");

	double *temp = malloc_1d_array_double(Arows*Bcols);
	
	for (long int row=0; row<Arows; row++) {
		for (long int col=0; col<Bcols; col++) {
			temp[row*Bcols+col]=0;
			for (long int sum_point=0; sum_point<Acols; sum_point++) temp[row*Bcols+col]+=A[row*Acols+sum_point]*B[sum_point*Bcols+col];
		}
	}

	for (long int row=0; row<Arows; row++)
		for (long int col=0; col<Bcols; col++)
			product[row*Bcols+col]=temp[row*Bcols+col];

	free_1d_array_double(temp);
}

// ****************************************************************************
// Swaps rows row1 and row2 of matrix A (of size N x N)
// ****************************************************************************
void single_swap_row(long int row1, long int row2, double *A, long int N)
{
	if (A==NULL || N<1 || row1<0 || row1>=N || row2<0 || row2>=N) quit_error((char*)"Bad input data to single_swap_row");
	double temp;
	for (long int col=0; col<N; col++) {
		temp=A[row1*N+col];
		A[row1*N+col]=A[row2*N+col];
		A[row2*N+col]=temp;
	}
}

// ****************************************************************************
// Multplies one row of matrix A (of size N x N) 
// by a certain factor
// ****************************************************************************
void single_multiply_row(long int row, double factor, double *A, long int N)
{
	if (A==NULL || N<1 || row<0 || row>=N || fabs(factor)<DBL_EPSILON) quit_error((char*)"Bad input data to single_multiply_row");
	for (long int col=0; col<N; col++) A[row*N+col]*=factor;
}

// ****************************************************************************
// Adds factor times row srow to row drow, in both matricies A and B 
// (each of size N x N)
// ****************************************************************************
void single_add_multiple_row(double factor, long int srow, long int drow, double *A, long int N)
{
	if (A==NULL || N<1 || srow<0 || srow>=N || drow<0 || drow>=N) quit_error((char*)"Bad input data to single_add_multiple_row");
	for (long int col=0; col<N; col++) A[drow*N+col]+=factor*A[srow*N+col];
}

// ****************************************************************************
// Calculate the determinant of a square matrix
// ****************************************************************************
double matrix_determinant(double *A, long int N)
{
	double *backup;
	double result;
	long int max_pivot;
	double max_value;

	if (N<1) quit_error((char*)"Attempt to calculate determinant of invalid matrix");
	if (A==NULL) quit_error((char*)"Attempt to calculate determinant of invalid matrix");

	// Copy the data to a place where it can safely be destroyed
	backup=(double *)malloc(sizeof(double)*N*N);
	if (backup==NULL) quit_error((char*)"Out of memory");

	for (long int row=0; row<N; row++)
		for (long int col=0; col<N; col++)
			backup[row*N+col]=A[row*N+col];
	
	result=1.0;

	// Go through the matrix, diagonalising to form an upper triangular matrix.
	// As part of this process, the rows are normalised to make the values along the diagonal
	// equal to one. As factors are pulled out, they are multiplied onto result to keep track
	// of the value of the determinant.
	
	for (long int i=0; i<N; i++) {
		// Find the element in this column which has the largest absolute value
		max_pivot=-1;
		max_value=0.0;
		for (long int row=i; row<N; row++) {
			if (fabs(backup[row*N+i])>fabs(max_value)) {
				max_pivot=row;
				max_value=backup[row*N+i];
			}
		}
		if (max_pivot<0) return 0.0;
		if (max_pivot!=i) {
			single_swap_row(max_pivot, i, backup, N);
			result*=-1.0;
		}
		result*=backup[i*N+i];
		single_multiply_row(i, 1.0/backup[i*N+i], backup, N);
		for (long int row=i+1; row<N; row++) single_add_multiple_row(-backup[row*N+i], i, row, backup, N);
	}

	free(backup);

	return result;
}

// ****************************************************************************
// Copies the contents of matrix A (R rows by C columns) onto the
// block of memory pointed to by Cp. It is assumed that sufficient 
// memory has been allocated to do this. 
// ****************************************************************************
void matrix_copy(double *A, long int R, long int C, double *Cp)
{
	if (A==NULL || Cp==NULL || R<1 || C<1) quit_error((char*)"Bad input to matrix_copy");
	
	for (long int row=0; row<R; row++)
		for (long int col=0; col<C; col++)
			Cp[row*C+col]=A[row*C+col];
}

// ****************************************************************************
// Finds the transpose of matrix A (R rows by C columns). 
// This is done in place. 
// The transposed matrix has C rows by R columns. 
// ****************************************************************************
void matrix_transpose(double *A, long int R, long int C)
{
	double *backup;

	if (A==NULL || R<1 || C<1) quit_error((char*)"Bad input to matrix_transpose");

	// Backup the A matrix so that manipulations can be performed here
	// without damaging the original matrix.
	backup=(double *)malloc(sizeof(double)*R*C);
	if (backup==NULL) quit_error((char*)"Out of memory");

	for (long int row=0; row<R; row++)
		for (long int col=0; col<C; col++)
			backup[row*C+col]=A[row*C+col];

	// Calculate the transpose
	for (long int row=0; row<C; row++)
		for (long int col=0; col<R; col++)
			A[row*R+col]=backup[col*C+row];
	
	free(backup);

	return;
}

// ****************************************************************************
// Internal function for calculating the Cholesky decomposition of a matrix.
// This code is from "Numerical Recipes in C : The Art of Scientific Computing" by Numerical Recipes Software, 
// Cambridge University Press, 1992. There have been some slight modifications made to it.
// Below is the original comment. Some of the array references are different in this version of the code:
// Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
// decomposition, A = L ∑ LT . On input, only the upper triangle of a need be given; it is not
// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
// elements which are returned in p[1..n].
// ****************************************************************************
int choldc(double *a, long int n, double p[])
{
	long int i,j,k;
	double sum;
	
	for (i=0; i<n; i++) 
	{
		for (j=i; j<n; j++)
		{
			for (sum=a[i*n+j],k=i-1;k>=0;k--) sum -= a[i*n+k]*a[j*n+k];
			if (i == j) 
			{
				// Check if a, with rounding errors, is not positive definite.
				if (sum <= 0.0)	return -1;
				p[i]=sqrt(sum);
			} else a[j*n+i]=sum/p[i];
		}
	}

	return 0;
}

// ****************************************************************************
// Calculates the Cholesky decomposition of N x N matrix A, such that A = L x transpose(L). 
// This is only possible if matrix A is symmetric and positive definite (all eigenvalues
// positive). 
// It is assumed that the memory required for this has already been allocated. 
// The data in matrix A is not destroyed, unless the same address is supplied for L.
// (This case should be successfully handled also.)
// The function returns -1 if the matrix is not symmetric or not positive definite. 
// ****************************************************************************
int matrix_cholesky(double *A, double *L, long int N)
{
	long int row, col;
	double *p, *backup;
	int result;
	
	if (A==NULL || L==NULL || N<1) quit_error((char*)"Bad input data to matrix_cholesky");
	
	// Check that the A matrix is at least symmetric. (It also needs to be positive definite
	// but that cannot be determined yet.)
	for (row=0; row<N; row++)
	{
		for (col=row; col<N; col++)
		{
			if (fabs(A[row*N+col]-A[col*N+row])>DBL_EPSILON) return -1;
		}
	}

	// Backup the A matrix so that manipulations can be performed here
	// without damaging the original matrix.
	backup=(double *)malloc(sizeof(double)*N*N);
	if (backup==NULL) quit_error((char*)"Out of memory");

	for (row=0; row<N; row++)
	{
		for (col=0; col<N; col++)
		{
			backup[row*N+col]=A[row*N+col];
		}
	}

	// Create the vector for the results along the main diagonal
	p=(double *)malloc(sizeof(double)*N);
	if (p==NULL) quit_error((char*)"Out of memory");

	// Calculate the Cholesky decomposition
	result=choldc(backup, N, p);

	// Assemble the answer in a nicer format
	for (row=0; row<N; row++)
	{
		for (col=0; col<N; col++)
		{
			if (row<col) // Upper
			{
				L[row*N+col]=0.0;
			}
			else if (row>col) // Lower
			{
				L[row*N+col]=backup[row*N+col];
			}
			else L[row*N+col]=p[row];  // Diagonal
		}
	}
	
	free(p);
	free(backup);

	return result;
}

// ****************************************************************************
// Invert the matrix A (size N by N) and store the result in Ainv, using the Cholesky
// decomposition. This is more numerically stable than standard inversion methods
// like Gauss-Jordan elimination but will only work if the matrix is symmetric
// and positive definite. 
// ****************************************************************************
int matrix_cholesky_invert(double *A, double *Ainv, long int N)
{
	double *L, *y, *z;
	int result;
	long int col, row, sum_count;
	double sum;

	L=(double *)malloc(sizeof(double)*N*N);
	if (L==NULL) quit_error((char*)"Out of memory");

	y=(double *)malloc(sizeof(double)*N);
	if (y==NULL) quit_error((char*)"Out of memory");

	z=(double *)malloc(sizeof(double)*N);
	if (z==NULL) quit_error((char*)"Out of memory");

	result=matrix_cholesky(A, L, N);

	// Generate the inverse matrix one column at a time
	for (col=0; col<N; col++)
	{
		// Create the y vector
		for (row=0; row<N; row++)
		{
			if (row==col) y[row]=1.0;
			else y[row]=0.0;
		}

		// Solve for the z vector
		for (row=0; row<N; row++)
		{
			sum=0.0;
			for (sum_count=0; sum_count<row; sum_count++)
			{
				sum+=L[row*N+sum_count]*z[sum_count];
			}
			z[row]=(y[row]-sum)/L[row*N+row];
		}

		// Solve for the solution vector - which is the next column
		// of Ainv. 
		for (row=N-1; row>=0; row--)
		{
			sum=0.0;
			for (sum_count=row+1; sum_count<N; sum_count++)
			{
				sum+=L[sum_count*N+row]*Ainv[sum_count*N+col];
			}
			Ainv[row*N+col]=(z[row]-sum)/L[row*N+row];
		}
	}

	free(z);
	free(y);
	free(L);
	
	return result;
}

// ****************************************************************************
// Multiply a matrix A (size Arows rows by Acols columns) by a constant, and store the result in product.
// product must point to an allocated area of memory of the correct size. It is allowed to point to the same area of 
// memory as A, and the contents of A will only be changed if this is the case. 
// ****************************************************************************
void matrix_constant_multiply(double *A, long int Arows, long int Acols, double constant, double *product)
{
	double *answer;
	long int row, col;

	if (A==NULL || Arows<1 || Acols <1 || product==NULL) quit_error((char*)"Illegal arguments to matrix_constant_multiply");
	
	answer=(double *)malloc(sizeof(double)*Arows*Acols);
	if (answer==NULL) quit_error((char*)"Out of memory");

	for (row=0; row<Arows; row++)
	{
		for (col=0; col<Acols; col++)
		{
			answer[row*Acols+col]=constant*A[row*Acols+col];
		}
	}
	
	for (row=0; row<Arows; row++)
	{
		for (col=0; col<Acols; col++)
		{
			product[row*Acols+col]=answer[row*Acols+col];
		}
	}

	free(answer);
	
	return;
}

// ****************************************************************************
// Adds matrix A (Arows by Acols) to matrix B (Brows by Bcols) and stores the result in sum. 
// sum must point to an allocated area of memory of the correct size. It is allowed to point to the same area of 
// memory as A and/or B, and the contents of A and/or B will only be changed if this is the case. 
// ****************************************************************************
void matrix_add(double *A, long int Arows, long int Acols, double *B, long int Brows, long int Bcols, double *sum)
{
	double *answer;
	long int row, col;

	if (A==NULL || Arows<1 || Acols <1 || B==NULL || Brows<1 || Bcols <1 || sum==NULL || Arows!=Brows || Acols!=Bcols) quit_error((char*)"Illegal arguments to matrix_add");

	answer=(double *)malloc(sizeof(double)*Arows*Acols);
	if (answer==NULL) quit_error((char*)"Out of memory");

	for (row=0; row<Arows; row++)
	{
		for (col=0; col<Acols; col++)
		{
			answer[row*Acols+col]=A[row*Acols+col]+B[row*Acols+col];
		}
	}
	
	for (row=0; row<Arows; row++)
	{
		for (col=0; col<Acols; col++)
		{
			sum[row*Acols+col]=answer[row*Acols+col];
		}
	}

	free(answer);

	return;
}


// ****************************************************************************
// The following code is concerned with calculating the singular value decomposition of a matrix.
// This code is from "Numerical Recipes in C : The Art of Scientific Computing" by Numerical Recipes Software, 
// Cambridge University Press, 1992. There have been some slight modifications made to it, to make it work in this
// new context. A wrapper function has been written to make the function easier to call, using the standard 
// conventions of this matrix module. 
// ****************************************************************************

#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) quit_error((char*)"allocation failure in dvector()");
	return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) quit_error((char*)"allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) quit_error((char*)"allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

double dpythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQ(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQ(absa/absb)));
}

/*
Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, 
A = U.W.transpose(V). The matrix U replaces a on output. The diagonal matrix of singular values
W is output as a vector w[1..n]. The matrix V (not its transpose) is output as v[1..n][1..n]. 
The matrix a must exist prior to calling this function (obviously), and so too must w and v. 
*/
void dsvdcmp(double **a, int m, int n, double w[], double **v)
{
	int flag,i,its,j,jj,k,l=0,nm=0;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=dvector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) quit_error((char*)"no convergence in 30 dsvdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_dvector(rv1,1,n);

	return;
}

// ****************************************************************************
// Compute the (thin) singular value decomposition of a matrix A (size: rows x cols).
// I.e. it finds U, D and V such that:
// U.D.transpose(V) = A
// V.transpose(V) = I
// U.transpose(U) = I (if rows = cols)
// The matricies U, D, V must already be allocated prior to calling this function.
// Their required sizes are: U (rows x cols), D (cols x cols), V(cols x cols)
// The contents of A are not changed. 
// ****************************************************************************
void svd(double *A, long int rows, long int cols, double *U, double *D, double *V)
{
	double **a, **v, *d;
	long int r, c;

	// Allocate the data structures which the core function expects
	a=dmatrix(1, rows, 1, cols);
	v=dmatrix(1, cols, 1, cols);
	d=dvector(1, cols);
	
	// Copy the input data into the core data structures
	for (r=1; r<=rows; r++)
	{
		for (c=1; c<=cols; c++)
		{
			a[r][c]=A[(r-1)*cols+(c-1)];
		}
	}

	// Call the core function
	dsvdcmp(a, rows, cols, d, v);
	
	// Copy the output data into the output data structures
	for (r=1; r<=rows; r++)
	{
		for (c=1; c<=cols; c++)
		{
			U[(r-1)*cols+(c-1)]=a[r][c];
		}
	}

	for (r=1; r<=cols; r++)
	{
		for (c=1; c<=cols; c++)
		{
			V[(r-1)*cols+(c-1)]=v[r][c];
		}
	}

	for (r=1; r<=cols; r++)
	{
		for (c=1; c<=cols; c++)
		{
			if (r==c) D[(r-1)*cols+(c-1)]=d[c];
			else D[(r-1)*cols+(c-1)]=0.0;
		}
	}
	
	// Free the core data structures
	free_dvector(d, 1, cols);
	free_dmatrix(v, 1, cols, 1, cols);
	free_dmatrix(a, 1, rows, 1, cols);
	
	return;
}

// ****************************************************************************
// Row reduce matrix A (size: rows x cols). 
// This function WILL NOT handle the case when rows > cols. 
// false is returned if the matrix does not have the maximum possible rank, true otherwise. 
// tolerance is the largest number which can be considered to be equal to zero, for the purposes of pivoting.
// ****************************************************************************
int matrix_row_reduce(double *A, long int rows, long int cols, double tolerance)
{
	long int r, i, pivot_row, col_offset;
	double pivot_val;
	
	if (A==NULL || rows<1 || cols<1) quit_error((char*)"Bad arguments to matrix_row_reduce");
	if (rows>cols) quit_error((char*)"matrix_row_reduce is limited to cases where rows<=cols");
	
	col_offset=0;
	for (i=0; i<rows; i++)
	{
		// Find the row with the best pivot and swap it
		pivot_val=fabs(tolerance);
		pivot_row=-1;
		while (pivot_row<0)
		{
			for (r=i; r<rows; r++)
			{
				if (fabs(A[r*cols+i+col_offset])>fabs(pivot_val))
				{
					pivot_val=A[r*cols+i+col_offset];
					pivot_row=r;
				}
			}
			if (pivot_row<0)
			{
				col_offset++;
				if (i+col_offset>=cols) return false;
			}
		}
		if (pivot_row>=0 && pivot_row<rows)
		{

			if (pivot_row!=i) single_swap_row(pivot_row, i, A, cols);
			single_multiply_row(i, 1.0/pivot_val, A, cols);

			for (r=0; r<rows; r++)
			{
				if (r!=i)
				{
					single_add_multiple_row(-A[r*cols+i+col_offset], i, r, A, cols);
				}
			}
		}

	}
	
	return true;
}

// ****************************************************************************
// This function takes a matrix A (size: rows x cols) and calculates useful information for solving the problem A . p = 0, where
// p is a vector of size cols x 1. 
// In fact, it finds a diagonal matrix D and square matrix V such that D . V . p = 0 and the solution space for p will be
// the same as that of A . p = 0. 
// D is a diagonal matrix of size cols x cols, with the singular values down the main diagonal (sorted from largest magnitude to smallest).
// V is a square matrix of size cols x cols, with the rows of the matrix being the corresponding singular unit vectors (which are mutually orthogonal).
// The nullity of A is returned, based on the specified tolerance. 
// ****************************************************************************
long int find_null_space_decomposition(double *A, long int rows, long int cols, double **D, double **V, double tolerance)
{
	double *Utemp, *Dtemp, *Vtemp, max_val, temp_val;
	long int nullity=cols, i, ii, max_loc=0;

	Utemp=(double *)malloc(sizeof(double)*rows*cols);
	if (Utemp==NULL) quit_error((char*)"Out of memory");
	Dtemp=(double *)malloc(sizeof(double)*cols*cols);
	if (Dtemp==NULL) quit_error((char*)"Out of memory");
	Vtemp=(double *)malloc(sizeof(double)*cols*cols);
	if (Vtemp==NULL) quit_error((char*)"Out of memory");
	
	// Calculate the Singular Value Decomposition
	svd(A, rows, cols, Utemp, Dtemp, Vtemp);

	// Transpose Vtemp
	matrix_transpose(Vtemp, cols, cols);

	// Sort the singular values into the correct order
	for (i=0; i<cols; i++)
	{
		// Find the maximum singular value in the remainder of the list
		max_val=-1.0;
		for (ii=i; ii<cols; ii++)
		{
			if (fabs(Dtemp[ii+cols*ii])>max_val)
			{
				max_val=fabs(Dtemp[ii+cols*ii]);
				max_loc=ii;
			}
		}

		// Update the nullity
		if (max_val>fabs(tolerance)) nullity--;

		// Swap the singular values
		temp_val=Dtemp[max_loc+cols*max_loc];
		Dtemp[max_loc+cols*max_loc]=Dtemp[i+cols*i];
		Dtemp[i+cols*i]=temp_val;
		
		// Swap the corresponding singular vectors
		single_swap_row(i, max_loc, Vtemp, cols);
	}
		
	// Return the stuff that is required and throw away the rest
	free(Utemp);
	*D=Dtemp;
	*V=Vtemp;

	return nullity;
}


// ****************************************************************************
// Same as the previous function, but does not return D and V
// ****************************************************************************
long int find_null_space_decomposition(double *A, long int rows, long int cols, double tolerance)
{
	double *Utemp, *Dtemp, *Vtemp, max_val, temp_val;
	long int nullity=cols, i, ii, max_loc=0;
	
	Utemp=(double *)malloc(sizeof(double)*rows*cols);
	if (Utemp==NULL) quit_error((char*)"Out of memory");
	Dtemp=(double *)malloc(sizeof(double)*cols*cols);
	if (Dtemp==NULL) quit_error((char*)"Out of memory");
	Vtemp=(double *)malloc(sizeof(double)*cols*cols);
	if (Vtemp==NULL) quit_error((char*)"Out of memory");
	
	// Calculate the Singular Value Decomposition
	svd(A, rows, cols, Utemp, Dtemp, Vtemp);
	
	// Transpose Vtemp
	matrix_transpose(Vtemp, cols, cols);
	
	// Sort the singular values into the correct order
	for (i=0; i<cols; i++)
	{
		// Find the maximum singular value in the remainder of the list
		max_val=-1.0;
		for (ii=i; ii<cols; ii++)
		{
			if (fabs(Dtemp[ii+cols*ii])>max_val)
			{
				max_val=fabs(Dtemp[ii+cols*ii]);
				max_loc=ii;
			}
		}
		
		// Update the nullity
		if (max_val>fabs(tolerance)) nullity--;
		
		// Swap the singular values
		temp_val=Dtemp[max_loc+cols*max_loc];
		Dtemp[max_loc+cols*max_loc]=Dtemp[i+cols*i];
		Dtemp[i+cols*i]=temp_val;
		
		// Swap the corresponding singular vectors
		single_swap_row(i, max_loc, Vtemp, cols);
	}
	
	// Return the stuff that is required and throw away the rest
	free(Utemp);
	free(Dtemp);
	free(Vtemp);
	
	return nullity;
}


// ****************************************************************************
// This function finds a set of vectors (P) such that, for each vector p:
// A.p = 0
// Note that A is a matrix of size rows x cols, P is a matrix of size cols x cols whose ROWS are the p vectors.
// P must NOT point to any allocated memory prior to calling this function. 
// The rank deficit (nullity) is returned (0 if matrix is full rank, 1 if matrix has one singular value, etc).
// tol defines how small a singular value must be in order to be considered a "zero". 
// ****************************************************************************
long int find_null_space_vectors(double *A, long int rows, long int cols, double **P, double tolerance)
{
	double *U, *D, *V, *R, max_val, min_val;
	long int i, rank_deficit, c, r;

	U=(double *)malloc(sizeof(double)*rows*cols);
	if (U==NULL) quit_error((char*)"Out of memory");
	D=(double *)malloc(sizeof(double)*cols*cols);
	if (D==NULL) quit_error((char*)"Out of memory");
	V=(double *)malloc(sizeof(double)*cols*cols);
	if (V==NULL) quit_error((char*)"Out of memory");
	
	// Calculate the Singular Value Decomposition
	svd(A, rows, cols, U, D, V);

	// Find the maximum and minimum magnitude singular values
	max_val=0.0;
	min_val=DBL_MAX;
	for (i=0; i<cols; i++)
	{
		if (fabs(D[i*cols+i])<fabs(min_val)) min_val=fabs(D[i*cols+i]);
		if (fabs(D[i*cols+i])>fabs(max_val)) max_val=fabs(D[i*cols+i]);
	}

	// Count the rank deficit of the matrix, and extract
	// relevant rows of the transpose(V) matrix. 
	rank_deficit=0;
	for (i=0; i<cols; i++)
	{
		if (fabs(D[i*cols+i])<tolerance) rank_deficit++;
	}
	
	// Transpose V
	matrix_transpose(V, cols, cols);

	// Generate a reduced version of transpose(V)
	R=(double *)malloc(sizeof(double)*cols*(rank_deficit));
	if (R==NULL) quit_error((char*)"Out of memory");
	r=0;
	for (i=0; i<cols; i++)
	{
		if (fabs(D[i*cols+i])<fabs(tolerance))
		{
			for (c=0; c<cols; c++)
			{
				R[r*cols+c]=V[i*cols+c];
			}
			r++;
		}
	}

	// Row reduce R
	matrix_row_reduce(R, rank_deficit, cols, tolerance);

	*P=R;
	
	free(V);
	free(D);
	free(U);
	
	return rank_deficit;
}

// ****************************************************************************
// This function finds the nullity of matrix A.
// Note that A is a matrix of size rows x cols 
// The rank deficit (nullity) is returned (0 if matrix is full rank, 1 if matrix has one singular value, etc).
// tolerance defines how small a singular value must be in order to be considered a "zero". 
// ****************************************************************************
long int matrix_nullity(double *A, long int rows, long int cols, double tolerance)
{
	double *U, *D, *V, max_val, min_val;
	long int i, rank_deficit;

	U=(double *)malloc(sizeof(double)*rows*cols);
	if (U==NULL) quit_error((char*)"Out of memory");
	D=(double *)malloc(sizeof(double)*cols*cols);
	if (D==NULL) quit_error((char*)"Out of memory");
	V=(double *)malloc(sizeof(double)*cols*cols);
	if (V==NULL) quit_error((char*)"Out of memory");
	
	// Calculate the Singular Value Decomposition
	svd(A, rows, cols, U, D, V);

	// Find the maximum and minimum magnitude singular values
	max_val=0.0;
	min_val=DBL_MAX;
	for (i=0; i<cols; i++)
	{
		if (fabs(D[i*cols+i])<fabs(min_val)) min_val=D[i*cols+i];
		if (fabs(D[i*cols+i])>fabs(max_val)) max_val=D[i*cols+i];
	}

	// Count the rank deficit of the matrix, and extract
	// relevant rows of the transpose(V) matrix. 
	rank_deficit=0;
	for (i=0; i<cols; i++)
	{
		if (fabs(D[i*cols+i])<fabs(tolerance)) rank_deficit++;
	}
	
	free(V);
	free(D);
	free(U);
	
	return rank_deficit;
}

