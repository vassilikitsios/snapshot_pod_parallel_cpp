//==========
//misc.h
//==========

#include <cstring>

//****************************************************************************
// defines
//****************************************************************************
//#define PI 3.14159265358979
#define MAX_VAL(X,Y) ((X)>=(Y)?(X):(Y))
#define MIN_VAL(X,Y) ((X)<=(Y)?(X):(Y))
#define MAX(X,Y) (((X) > (Y) ? (X) : (Y)))
#define MIN(X,Y) (((X) < (Y) ? (X) : (Y)))

#define true	1
#define false	0

//#define PI 3.1415926535897932384626433832795
#define PI			M_PI
//#define E 2.7182818284590452353602874713527
#define E			M_E

#define BIG_NUMBER		HUGE_VAL
#define ZERO				1.0e-14
#define SMALL_NUMBER		-HUGE_VAL

#define SQ(X) ( (X) * (X) )
#define CU(X) ( (X) * (X) * (X) )
#define ABS(X) sqrt( SQ(X) )
#define SIGN(X) ABS(X)/X

// The maximum allowable line length whenever fgets or whatever is used to read
// from a file. 
#define MAX_LINE_LENGTH 1000

// The maximum length for a string which is expected to be fairly short
#define MAX_STRING_LENGTH 1000

// Convert a string to upper case
#define ucase(STR) {long int __ucase_i;for(__ucase_i=0; __ucase_i<strlen(STR); __ucase_i++) STR[__ucase_i]=toupper(STR[__ucase_i]);}


//****************************************************************************
// general
//****************************************************************************

// Terminates the program with an error message. The arguments are in the same format as printf().
extern void quit_error(char *fmt, ...);
void quit_error_mpi(int myrank, int nprocs, char *fmt, ...);

// Prints a warning to stderr. The arguments are in the same format as printf().
extern void print_warning(char *fmt, ...);
