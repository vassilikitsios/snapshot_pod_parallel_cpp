//==========================
//binary_read_or_write.h
//==========================

using namespace std;

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>

enum ReadOrWrite {Read, Write};

//****************************************************************************

extern void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, int *val);
extern void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, double *val);
extern void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, char *val);
extern void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, bool *val);

