// tridiagonal_solver.cpp
#include <cstdio>
#include <cstdlib>
#include "memory_allocations.h"
#include "tridiagonal_solver.h"

// static inline functions...
static inline void inline_tridiagonal_solver_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd);
static inline void inline_tridiagonal_solver_periodic_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd, double *C);

// ****************************************************************************
// A tri-diagonal system of equations is:
//
// | b0 c0                              |   | x0     |   | d0     | 
// | a1 b1 c1                           |   | x1     |   | d1     | 
// |    a2 b2 c2                        |   | ...    |   |        |
// |     .......                        | x | ...    | = |        | 
// |       .......                      |   | ...    |   |        | 
// |        a(N-3) b(N-3) c(N-3)        |   | ...    |   |        | 
// |               a(N-2) b(N-2) c(N-2) |   | x(N-1) |   | d(N-1) |  
// |                      a(N-1) b(N-1) |   | x(N-1) |   | d(N-1) | 
//
// ****************************************************************************



// ****************************************************************************
// tridiagonal_solver_nondestructive()
// This does NOT over-write the input matrices...
// ****************************************************************************
static inline void inline_tridiagonal_solver_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd)
{
	int i;
	
	// --------------------
	// GAUSSIAN ELIMINATION
	// --------------------
	bb[0] = b[0];
	dd[0] = d[0];
	for(i=1 ; i<N ; i++)
	{
		bb[i] = b[i] - a[i]/bb[i-1] * c[i-1];
		dd[i] = d[i] - a[i]/bb[i-1] * dd[i-1];
	}
	
	// ---------------------
	// BACKWARD SUBSTITUTION
	// ---------------------
	x[N-1] = dd[N-1] / bb[N-1];
	for(i=N-2 ; i>=0 ; i--)
	{
		x[i] = (dd[i] - c[i]*x[i+1]) / bb[i];
	}
	
	return;
}
// ****************************************************************************
// tridiagonal_solver_nondestructive()
// ****************************************************************************
void tridiagonal_solver_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd)
{
	inline_tridiagonal_solver_nondestructive(N, a, b, c, d, x, bb, dd);
	return;
}

// ****************************************************************************
// tridiagonal_solver_nondestructive_static()
// ****************************************************************************
void tridiagonal_solver_nondestructive_static(int N, double *a, double *b, double *c, double *d, double *x)
{
	static int Nbefore = 10;
	static double *bb = malloc_1d_array_double(Nbefore);
	static double *dd = malloc_1d_array_double(Nbefore);
	if(N > Nbefore)
	{
		Nbefore = N;
		free(bb); bb = malloc_1d_array_double(N);
		free(dd); dd = malloc_1d_array_double(N);
	}
	
	inline_tridiagonal_solver_nondestructive(N, a, b, c, d, x, bb, dd);
	
	return;
}


// ****************************************************************************
// tridiagonal_solver_destructive()
// This does OVER-WRITES the input matrices...
// ****************************************************************************
void tridiagonal_solver_destructive(int N, double *a, double *b, double *c, double *d, double *x)
{
	int i;
	
	// --------------------
	// GAUSSIAN ELIMINATION
	// --------------------
	for(i=1 ; i<N ; i++)
	{
		b[i] = b[i] - a[i]*c[i-1]/b[i-1];
		d[i] = d[i] - a[i]*d[i-1]/b[i-1];
	}
	
	// ---------------------
	// BACKWARD SUBSTITUTION
	// ---------------------
	x[N-1] = d[N-1] / b[N-1];
	for(i=N-2 ; i>=0 ; i--)
	{
		x[i] = (d[i] - c[i]*x[i+1])/b[i];
	}
	
	return;
}

// ****************************************************************************
// A periodic tri-diagonal system of equations is:
//
// | b0 c0                           a0 |   | x0     |   | d0     | 
// | a1 b1 c1                           |   | x1     |   | d1     | 
// |    a2 b2 c2                        |   | ...    |   |        |
// |     .......                        | x | ...    | = |        | 
// |       .......                      |   | ...    |   |        | 
// |        a(N-3) b(N-3) c(N-3)        |   | ...    |   |        | 
// |               a(N-2) b(N-2) c(N-2) |   | x(N-1) |   | d(N-1) |  
// | c(N-1)               a(N-1) b(N-1) |   | x(N-1) |   | d(N-1) | 
//
// ****************************************************************************



// ****************************************************************************
// inline_tridiagonal_solver_periodic_nondestructive()
// This does NOT over-write the input matrices...
// ****************************************************************************
static inline void inline_tridiagonal_solver_periodic_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd, double *C)
{
	double R;
	int i, j;
	
	// --------------------
	// GAUSSIAN ELIMINATION
	// --------------------
	
	// initialisation...
	bb[0] = b[0];
	dd[0] = d[0];
    C[0] = a[0];
	//aa[N-1] = a[N-1];
	bb[N-1] = b[N-1];
	dd[N-1] = d[N-1];
	R = c[N-1];
	
    for(i=0 ; i<=N-4 ; i++)
	{
		// ELIMINATE a[i+1]:
		// row_(i+1) -> row_(i+1) + a(i+1)/b(i) * row(i)
		bb[i+1] = b[i+1] - a[i+1]/bb[i] * c[i];
		dd[i+1] = d[i+1] - a[i+1]/bb[i] * dd[i];
		C[i+1] = 0       - a[i+1]/bb[i] * C[i];
		// ELIMINATE R:
		// row_(N-1) -> row_(N-1) + R/b(i) * row(i)
        bb[N-1] = bb[N-1] - R/bb[i] * C[i];
		dd[N-1] = dd[N-1] - R/bb[i] * dd[i];
		R = 0             - R/bb[i] * c[i];
    }
	i = N-3;
	{
		// ELIMINATE a[i+1]:
		// row_(i+1) -> row_(i+1) + a(i+1)/b(i) * row(i)
		bb[i+1] = b[i+1] - a[i+1]/bb[i] * c[i];
		dd[i+1] = d[i+1] - a[i+1]/bb[i] * dd[i];
		C[i+1] = c[i+1] - a[i+1]/bb[i] * C[i];
		// ELIMINATE R:
		// row_(N-1) -> row_(N-1) + R/b(i) * row(i)
        bb[N-1] = bb[N-1] - R/bb[i] * C[i];
		dd[N-1] = dd[N-1] - R/bb[i] * dd[i];
		R = a[N-1] - R/bb[i] * c[i];
	}
    i = N-2; 
	{
		// ELIMINATE a[N-1]:
		// row_(N-1) -> row_(N-1) + a(N-1)/b(i) * row(i)
		bb[N-1] = bb[N-1] - R/bb[i] * C[i];
		dd[N-1] = dd[N-1] - R/bb[i] * dd[i];
	}
	
    // ---------------------
	// BACKWARD SUBSTITUTION
	// ---------------------
	
	i = N-1;
	{
		// MAKE B[i] = 1:
		// row_(i) -> row_(i) / b(i)
		x[i] = dd[i] / bb[i];
		// ELIMINATE c[i-1]
		// row_(i-1) -> row_(i-1) - c(i-1) * row(i)
		for(j=i-1 ; j>= 0 ; j--)
		{
			// ELIMINATE C[j]
			// row_(j) -> row_(j) - C(j) * row(i)
			x[j] = dd[j] - C[j] * x[i];
		}
	}
	
	for (i=N-2 ; i>=1 ; i--)
	{
		// MAKE B[i] = 1:
		// row_(i) -> row_(i) / b(i)
		x[i] = x[i] / bb[i];
		// ELIMINATE c[i-1]
		// row_(i-1) -> row_(i-1) - c(i-1) * row(i)
		x[i-1] = x[i-1] - c[i-1] * x[i];
	}
	
	i = 0;
	{
		// MAKE B[i] = 1:
		x[i] = x[i] / bb[i];
	}
	
	return;
}

// ****************************************************************************
// tridiagonal_solver_periodic_nondestructive()
// This does NOT over-write the input matrices...
// ****************************************************************************
void tridiagonal_solver_periodic_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd, double *C)
{
	inline_tridiagonal_solver_periodic_nondestructive(N,a,b,c,d,x,bb,dd,C);
	return;
}

// ****************************************************************************
// tridiagonal_solver_periodic_nondestructive_static()
// This does NOT over-write the input matrices...
// ****************************************************************************
void tridiagonal_solver_periodic_nondestructive_static(int N, double *a, double *b, double *c, double *d, double *x)
{
	static int Nbefore = 10;
	static double *bb = malloc_1d_array_double(Nbefore);
	static double *dd = malloc_1d_array_double(Nbefore);
	static double *C = malloc_1d_array_double(Nbefore);
	if(N > Nbefore)
	{
		Nbefore = N;
		free(bb); bb = malloc_1d_array_double(N);
		free(dd); dd = malloc_1d_array_double(N);
		free(C) ; C  = malloc_1d_array_double(N);
	}
	
	inline_tridiagonal_solver_periodic_nondestructive(N,a,b,c,d,x,bb,dd,C);
	
	return;
}


// ****************************************************************************
// tridiagonal_solver_periodic_destructive()
// This OVER-WRITES the input matrices...
// ****************************************************************************
void tridiagonal_solver_periodic_destructive(int N, double *a, double *b, double *c, double *d, double *x)
{
	double R;
	int i, j;
	
	// --------------------
	// GAUSSIAN ELIMINATION
	// --------------------
	
	// initialisation...
	R = c[N-1];
	
    for(i=0 ; i<=N-4 ; i++)
	{
		// ELIMINATE a[i+1]:
		// row_(i+1) -> row_(i+1) + a(i+1)/b(i) * row(i)
		b[i+1] = b[i+1] - a[i+1]/b[i] * c[i];
		d[i+1] = d[i+1] - a[i+1]/b[i] * d[i];
		//C[i+1] = 0      - a[i+1]/b[i] * C[i];
		a[i+1] = 0      - a[i+1]/b[i] * a[i];
		// ELIMINATE R:
		// row_(N-1) -> row_(N-1) + R/b(i) * row(i)
        b[N-1] = b[N-1] - R/b[i] * a[i];
		d[N-1] = d[N-1] - R/b[i] * d[i];
		R = 0           - R/b[i] * c[i];
    }
	i = N-3;
	{
		// ELIMINATE a[i+1]:
		// row_(i+1) -> row_(i+1) + a(i+1)/b(i) * row(i)
		b[i+1] = b[i+1] - a[i+1]/b[i] * c[i];
		d[i+1] = d[i+1] - a[i+1]/b[i] * d[i];
		c[i+1] = c[i+1] - a[i+1]/b[i] * a[i];
		// ELIMINATE R:
		// row_(N-1) -> row_(N-1) + R/b(i) * row(i)
        b[N-1] = b[N-1] - R/b[i] * a[i];
		d[N-1] = d[N-1] - R/b[i] * d[i];
		a[N-1] = a[N-1] - R/b[i] * c[i];
	}
    i = N-2;
	{
		// ELIMINATE a[N-1]:
		// row_(N-1) -> row_(N-1) + a(N-1)/b(i) * row(i)
		b[N-1] = b[N-1] - a[N-1]/b[i] * c[i];
		d[N-1] = d[N-1] - a[N-1]/b[i] * d[i];
	}

    // ---------------------
	// BACKWARD SUBSTITUTION
	// ---------------------
	
	i = N-1;
	{
		// MAKE B[i] = 1:
		// row_(i) -> row_(i) / b(i)
		x[i] = d[i] / b[i];
		// ELIMINATE c[i-1]
		// row_(i-1) -> row_(i-1) - c(i-1) * row(i)
		x[i-1] = d[i-1] - c[i-1] * x[i];	
		for(j=i-2 ; j>= 0 ; j--)
		{
			// ELIMINATE C[j]
			// row_(j) -> row_(j) - C(j) * row(i)
			x[j] = d[j] - a[j] * x[i];
		}
	}
	
	for (i=N-2 ; i>=1 ; i--)
	{
		// MAKE B[i] = 1:
		// row_(i) -> row_(i) / b(i)
		x[i] = x[i] / b[i];
		// ELIMINATE c[i-1]
		// row_(i-1) -> row_(i-1) - c(i-1) * row(i)
		x[i-1] = x[i-1] - c[i-1] * x[i];
	}
	
	i = 0;
	{
		// MAKE B[i] = 1:
		x[i] = x[i] / b[i];
	}
	
	return;
}

// ****************************************************************************
// tridiagonal_solver_periodic_test()
// ****************************************************************************
void tridiagonal_solver_periodic_test(void)
{
	const int N=8;
	double a[N] = {1, 2, 3, 4, 5, 3, 2, 1};
	double b[N] = {2, 3, 4, 2, 2, 1, 2, 2};
	double c[N] = {6, 5, 3, 5, 4, 3, 3, 4};
	double d[N] = {4, 3, 3, 3, 3, 2, 1, 2};
	double a1[N] = {1, 2, 3, 4, 5, 3, 2, 1};
	double b1[N] = {2, 3, 4, 2, 2, 1, 2, 2};
	double c1[N] = {6, 5, 3, 5, 4, 3, 3, 4};
	double d1[N] = {4, 3, 3, 3, 3, 2, 1, 2};
	double C[N];
	double x[3][N];
	double bb[N], dd[N];
	
	printf("-----------------------------------\n");
	printf("Running tridiagonal_solver_periodic_test()...\n");
	printf("a = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",a[i]);	printf("}\n");
	printf("b = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",b[i]);	printf("}\n");
	printf("c = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",c[i]);	printf("}\n");
	printf("d = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",d[i]);	printf("}\n");

	tridiagonal_solver_periodic_nondestructive(N,a,b,c,d,x[0],bb,dd,C);
	tridiagonal_solver_periodic_nondestructive_static(N,a,b,c,d,x[1]);
	tridiagonal_solver_periodic_destructive(N,a1,b1,c1,d1,x[2]);

	printf("x_periodic             = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",x[0][i]);	printf("}\n");
	printf("x_periodic_static      = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",x[1][i]);	printf("}\n");
	printf("x_periodic_destructive = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",x[2][i]);	printf("}\n");
	printf("x_periodic_matlab      = { 0.4095 0.5196 0.1244 0.3145 0.3747 0.1696 0.2355 0.0633 }\n");
	printf("-----------------------------------\n");
	
	return;
}


// ****************************************************************************
// tridiagonal_solver_test()
// ****************************************************************************
void tridiagonal_solver_test(void)
{
	const int N=8;
	double a[N] = {0, 2, 3, 4, 5, 3, 2, 1};
	double b[N] = {2, 3, 4, 2, 2, 1, 2, 2};
	double c[N] = {6, 5, 3, 5, 4, 3, 3, 0};
	double d[N] = {4, 3, 3, 3, 3, 2, 1, 2};
	double a1[N] = {0, 2, 3, 4, 5, 3, 2, 1};
	double b1[N] = {2, 3, 4, 2, 2, 1, 2, 2};
	double c1[N] = {6, 5, 3, 5, 4, 3, 3, 0};
	double d1[N] = {4, 3, 3, 3, 3, 2, 1, 2};
	double C[N];
	double x[6][N];
	double bb[N], dd[N];
	
	printf("-----------------------------------\n");
	printf("Running tridiagonal_solver_test()...\n");
	printf("a = { "); for(int i=0 ; i<N ; i++) printf("%f ",a[i]);	printf("}\n");
	printf("b = { "); for(int i=0 ; i<N ; i++) printf("%f ",b[i]);	printf("}\n");
	printf("c = { "); for(int i=0 ; i<N ; i++) printf("%f ",c[i]);	printf("}\n");
	printf("d = { "); for(int i=0 ; i<N ; i++) printf("%f ",d[i]);	printf("}\n");

	tridiagonal_solver_periodic_nondestructive(N,a,b,c,d,x[0],bb,dd,C);
	tridiagonal_solver_periodic_nondestructive_static(N,a,b,c,d,x[1]);
	tridiagonal_solver_periodic_destructive(N,a,b,c,d,x[2]);
	tridiagonal_solver_nondestructive(N,a1,b1,c1,d1,x[3],bb,dd);
	tridiagonal_solver_nondestructive_static(N,a1,b1,c1,d1,x[4]);
	tridiagonal_solver_destructive(N,a1,b1,c1,d1,x[5]);	
	
	printf("x_periodic               = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",x[0][i]); printf("}\n");
	printf("x_periodic_static        = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",x[1][i]); printf("}\n");
	printf("x_periodic_destructive   = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",x[2][i]); printf("}\n");
	printf("x_non_destructive        = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",x[3][i]);	printf("}\n");
	printf("x_non_destructive_static = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",x[4][i]);	printf("}\n");
	printf("x_destructive            = { "); for(int i=0 ; i<N ; i++) printf("%.4f ",x[5][i]);	printf("}\n");
	printf("x_matlab                 = { 2.3540 -0.1180 -0.2708 1.4791 0.2250 -1.2114 0.8455 0.5773 }\n");
	printf("-----------------------------------\n");
	
	return;	
}


