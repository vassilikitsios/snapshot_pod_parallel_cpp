//==========================
//file_opening_closing.h
//==========================

//****************************************************************************
// files
//****************************************************************************
/** 
Opens a file in the given mode.
**/
extern FILE *open_file(char *filename, char *mode);
/** 
Opens a file in the given mode. 
**/
extern FILE *open_file(char *filename, char *mode, bool verbose);
/** 
Closes a file.
**/
extern void close_file(FILE *fp, char *filename);
/** 
Closes a file.
**/
extern void close_file(FILE *fp, char *filename, bool verbose);


