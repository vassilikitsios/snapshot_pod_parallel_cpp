// tridiagonal_solver.h

extern void tridiagonal_solver_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd);
extern void tridiagonal_solver_nondestructive_static(int N, double *a, double *b, double *c, double *d, double *x);
extern void tridiagonal_solver_destructive(int N, double *a, double *b, double *c, double *d, double *x);
extern void tridiagonal_solver_periodic_nondestructive(int N, double *a, double *b, double *c, double *d, double *x, double *bb, double *dd, double *C);
extern void tridiagonal_solver_periodic_nondestructive_static(int N, double *a, double *b, double *c, double *d, double *x);
extern void tridiagonal_solver_periodic_destructive(int N, double *a, double *b, double *c, double *d, double *x);

extern void tridiagonal_solver_test(void);
extern void tridiagonal_solver_periodic_test(void);


