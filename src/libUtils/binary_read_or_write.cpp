// binary_read_or_write.cpp

#include "misc.h"
#include "binary_read_or_write.h"

//****************************************************************************
//binary_read_or_write()
//****************************************************************************
void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, int *val)
{
	const size_t size = sizeof(int);
	if(read_or_write == Write)
	{
		fwrite((const void *)val, size, count, fp);
	}
	else
	{
		fread((void *)val, size, count, fp);
	}
	return;
}

//****************************************************************************
//binary_read_or_write()
//****************************************************************************
void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, double *val)
{
	const size_t size = sizeof(double);
	if(read_or_write == Write)
	{
		fwrite((const void *)val, size, count, fp);
	}
	else
	{
		fread((void *)val, size, count, fp);
	}
	return;
}

//****************************************************************************
//binary_read_or_write()
//****************************************************************************
void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, char *val)
{
	const size_t size = sizeof(char);
	if(read_or_write == Write)
	{
		fwrite((const void *)val, size, count, fp);
	}
	else
	{
		fread((void *)val, size, count, fp);
	}
	return;
}

//****************************************************************************
//binary_read_or_write()
//****************************************************************************
void binary_read_or_write(FILE *fp, ReadOrWrite read_or_write, int count, bool *val)
{
	const size_t size = sizeof(bool);
	if(read_or_write == Write)
	{
		fwrite((const void *)val, size, count, fp);
	}
	else
	{
		fread((void *)val, size, count, fp);
	}
	return;
}


