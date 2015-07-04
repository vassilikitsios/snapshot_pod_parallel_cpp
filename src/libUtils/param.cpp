//param.cpp

using namespace std;

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "memory_allocations.h"
#include "misc.h"
#include "file_opening_closing.h"
#include "string.h"
#include "param.h"


//****************************************************************************
//Param()
//****************************************************************************
Param_t::Param_t(char *name, char *line)
{
	this->name = copy_string(name);
	this->line = copy_string(line);
	this->numberOfAccesses = 0;
	this->next = NULL;
}


//****************************************************************************
//addParam()
//****************************************************************************
void ParamListHead_t::addParam(char *param_line)
{
	Param_t *p;
	int i;
	char buf_name[MAX_PARAM_LINE_LENGTH+1];
	char buf_line[MAX_PARAM_LINE_LENGTH+1];
	
	//********************************
	// Find the name of the param line
	//********************************
	i = 0;
	// Skip leading white space
	while( *param_line == ' ' || *param_line == '\t' )
		param_line++;
	// Get the name of the command
	while( *param_line != '\0' && *param_line != ' ' && *param_line != '=' && i < MAX_PARAM_LINE_LENGTH )
	{
		buf_name[i++] = *param_line;
		param_line++;
	}
	buf_name[i] = '\0';
	
	// Is this a blank line?
	if(i==0 && *param_line == '\0')
		return;
	
	// Error checking
	if(i==0) // If i==0 for any other reason there must be a problem.
		quit_error((char*)"length of buf_name = 0 in addParam()");
	if(i >= MAX_PARAM_LINE_LENGTH)
		quit_error((char*)"length of buf_name >= MAX_PARAM_LINE_LENGTH in addParam()");

	//***********************
	// Now find its arguments
	//***********************
	i = 0;
	// Skip leading white space and equals signs
	while( *param_line == ' ' || *param_line == '\t' || *param_line == '=' )
		param_line++;
	// Get the name of the command
	while( *param_line != '\0' && i < MAX_PARAM_LINE_LENGTH )
	{
		buf_line[i++] = *param_line;
		param_line++;
	}
	buf_line[i] = '\0';
	// Error checking
	//if(i == 0)
	//		quit_error((char*)"length of buf_line == 0 in addParam()");
	if(i >= MAX_PARAM_LINE_LENGTH)
		quit_error((char*)"length of buf_line >= MAX_PARAM_LINE_LENGTH in addParam()");
		
	
	//***********************
	// Make the new parameter
	//***********************
	p = new Param_t(buf_name, buf_line);
		
	// Place the new entry at the head of the list
	p->next = param_list;
	p->prev = NULL;

	// Update the prev pointer for the next node in the list
	if(this->param_list != NULL)
		this->param_list->prev = p;
		
	this->param_list = p;
	
}


//****************************************************************************
//print_params()
//****************************************************************************
void ParamListHead_t::print_params()
{
	Param_t *param_list = this->param_list;
	
	// Go to the end of the list
	if(param_list != NULL)
	{
		while(param_list->next != NULL)
		{
			param_list = param_list->next;
		}
	}
	
	
	fprintf(stdout,"\n*************** Param List ***************\n");
	while(param_list != NULL)
	{
		fprintf(stdout,"\"%s\" \"%s\"\n", param_list->name, param_list->line);
		param_list = param_list->prev;
	}
	fprintf(stdout,"************* End Param List *************\n\n");
	return;
}

//****************************************************************************
//print_param_accesses()
//****************************************************************************
void ParamListHead_t::print_param_accesses()
{
	Param_t *param_list = this->param_list;
	
	// Go to the end of the list
	if(param_list != NULL)
	{
		while(param_list->next != NULL)
		{
			param_list = param_list->next;
		}
	}
	
	
	fprintf(stdout,"\n************** Param Access Count **************\n");
	while(param_list != NULL)
	{
		fprintf(stdout,"\"%s\" \"%s\" %d\n", param_list->name, param_list->line, param_list->numberOfAccesses);
		param_list = param_list->prev;
	}
	fprintf(stdout,"************ End Param Access Count ************\n\n");	
	return;
}

//****************************************************************************
//get_params_from_file()
//****************************************************************************
void ParamListHead_t::get_params_from_file(char *filename)
{
	get_params_from_file(filename,true);
	return;
}

//****************************************************************************
//get_params_from_file()
//****************************************************************************
void ParamListHead_t::get_params_from_file(char *filename, bool verbose)
{
	char buf[MAX_PARAM_LINE_LENGTH+1];
	int i; // index into the buffer
	int c;
	FILE *fp;
	int finished = 0;
	int comment_line;

	fp = open_file(filename, (char*) "r", verbose);
	i = 0;
	comment_line = 0;
	while( finished == 0 )
	{
		c = fgetc(fp);
				
		if(c == EOF)
			finished = 1;
		
		if( i==0 && c == '#' )
			comment_line = 1;
		
		if(c == '\n' || c == EOF) // end of the current line?
		{
			if(i != 0 && comment_line == 0) // line isn't blank?
			{
				buf[i] = '\0';
				this->addParam(buf);
			}
			i = 0;
			comment_line = 0;
		}
		else // keep buffering 
		{
			if(i >= MAX_PARAM_LINE_LENGTH)
			{
				quit_error((char*)"Param line is too long in get_params_from_file()");
			}
			else
			{
				buf[i++] = c;
			}
		}
	}
	
	close_file(fp, filename, verbose);
	
	return;
}

//****************************************************************************
//get_params_from_command_line()
//****************************************************************************
void ParamListHead_t::get_params_from_command_line(int argc, char *argv[])
{
	char buf[MAX_PARAM_LINE_LENGTH+1];
	int i, j, len;
	
	//fprintf(stdout,"Reading in any parameters from the command line...\n");

	if(argc == 1)
		return;
	
	// skip over all the parameters until there is a "--"
	i = 1;
	if( strlen(argv[i]) < 2 || ! ( argv[i][0] == '-' && argv[i][1] == '-' ) )
	{
		i++;
		if( i >= argc )
			return;
	}
	
	// must have reached an argument with a --
	while( i < argc )
	{
		j = 0;
		len = strlen(argv[i]) - 2;
		if( j + len > MAX_PARAM_LINE_LENGTH )
			quit_error((char*)"Command line parameter is too long");
		strcpy(buf+j, argv[i]+2);
		j = len;
		
		i++;
		while( i < argc && ( strlen(argv[i]) < 2 || !(argv[i][0] == '-' && argv[i][1] == '-') ) ) 
		{
			len = strlen(argv[i]);
			if( j + 1 + len > MAX_PARAM_LINE_LENGTH )
				quit_error((char*)"Command line parameter is too long");
			buf[j] = ' ';
			strcpy(buf+j+1, argv[i]);
			j = j + 1 + len;
			i++;
		}
		
		if(j > 0) // 
			this->addParam(buf);
		
	}
	
	return;
}



//****************************************************************************
//check_param()
//This checks whether a particular parameter exists.
//****************************************************************************
bool ParamListHead_t::check_param(char *name)
{
	Param_t *param_list = this->param_list;
	
	while( param_list != NULL )
	{
		if(strcmp(param_list->name,name)==0)
		{
			param_list->numberOfAccesses++;
			
			return true;
		}
		else
		{
			param_list = param_list->next;
		}
	}
	return false;
}


//****************************************************************************
//get_int_param()
//****************************************************************************
int ParamListHead_t::get_int_param(char *name)
{
	
	if(this->param_list == NULL)
	{
		quit_error((char*)"param_list() is empty");
	}
	
	return this->param_list->get_int_param(name);
}

//****************************************************************************
//get_int_param()
//****************************************************************************
int Param_t::get_int_param(char *name)
{
	Param_t *param_list = this;
	
	while( param_list != NULL )
	{
		if(strcmp(param_list->name,name)==0)
		{
			param_list->numberOfAccesses++;
			return atoi(param_list->line);
		}
		else
		{
			param_list = param_list->next;
		}
	}
	fprintf(stderr,"In get_int_param, can't find param \"%s\"\n", name);
	quit_error((char*)"Unable to find param \"%s\"", name);
	return 0;
}

//****************************************************************************
//get_long_integer_param()
//****************************************************************************
int ParamListHead_t::get_long_int_param(char *name)
{
	
	if(this->param_list == NULL)
	{
		quit_error((char*)"param_list() is empty");
	}
	
	return this->param_list->get_long_int_param(name);
}

//****************************************************************************
//get_long_integer_param()
//****************************************************************************
int Param_t::get_long_int_param(char *name)
{
	Param_t *param_list = this;
	
	while( param_list != NULL )
	{
		if(strcmp(param_list->name,name)==0)
		{
			param_list->numberOfAccesses++;
			return atol(param_list->line);
		}
		else
		{
			param_list = param_list->next;
		}
	}
	fprintf(stderr,"In get_int_param, can't find param \"%s\"\n", name);
	quit_error((char*)"Unable to find param \"%s\"", name);
	return 0;
}


//****************************************************************************
//get_int_param()
//****************************************************************************
int ParamListHead_t::get_int_param(char *name, int defaultValue)
{
	
	if(this->param_list == NULL)
	{
		quit_error((char*)"param_list() is empty");
	}
	
	return this->param_list->get_int_param(name, defaultValue);
}

//****************************************************************************
//get_int_param()
//****************************************************************************
int Param_t::get_int_param(char *name, int defaultValue)
{
	Param_t *p;
	
	p = this->get_param(name);
	
	if(p == NULL)
	{
		printf("No value supplied, so using default value for %s = %d\n", name, defaultValue);
		return defaultValue;
	}
	else
		return p->get_int_param(name);
}

//****************************************************************************
//get_double_param()
//****************************************************************************
double ParamListHead_t::get_double_param(char *name)
{
	if(this->param_list == NULL)
	{
		quit_error((char*)"param_list() is empty");
	}
	
	return this->param_list->get_double_param(name);
}

//****************************************************************************
//get_nth_param()
//****************************************************************************
Param_t * ParamListHead_t::get_nth_param(int n, char *name)
{
	int i = 0;
	Param_t *p;
	
	if(this->param_list == NULL)
	{
		quit_error((char*)"param_list() is empty");
	}
	
	if( n<1 )
	{
		quit_error((char*)"n = %d < 1", n);
	}
	
	p = this->param_list->get_param(name);
	i++;
	
	
	while( i < n && p != NULL )
	{
		p = p->get_next_param(name);
		i++;
	}

	return p;
}

//****************************************************************************
//get_nth_int_param()
//****************************************************************************
int ParamListHead_t::get_nth_int_param(int n, char *name)
{
	Param_t *p = this->get_nth_param(n,name);
	if(p==NULL)	quit_error((char*)"Colour not find param '%s'", name);
	return atoi(p->line);
}

//****************************************************************************
//get_nth_long_int_param()
//****************************************************************************
long int ParamListHead_t::get_nth_long_int_param(int n, char *name)
{
	Param_t *p = this->get_nth_param(n,name);
	if(p==NULL)	quit_error((char*)"Colour not find param '%s'", name);
	return atol(p->line);
}

//****************************************************************************
//get_nth_double_param()
//****************************************************************************
double ParamListHead_t::get_nth_double_param(int n, char *name)
{
	Param_t *p = this->get_nth_param(n,name);
	if(p==NULL)	quit_error((char*)"Colour not find param '%s'", name);
	return atof(p->line);
}

//****************************************************************************
//get_double_param()
//****************************************************************************
double Param_t::get_double_param(char *name)
{
	Param_t *param_list = this;
	
	while( param_list != NULL )
	{
		if(strcmp(param_list->name,name)==0)
		{
			param_list->numberOfAccesses++;
			return atof(param_list->line);
		}
		else
		{
			param_list = param_list->next;
		}
	}
	fprintf(stderr,"In get_double_param, can't find param \"%s\"\n", name);
	quit_error((char*)"Unable to find param \"%s\"", name);
	return 0;
}

//****************************************************************************
//get_int_param()
//****************************************************************************
double ParamListHead_t::get_double_param(char *name, double defaultValue)
{
	if(this->param_list == NULL)
	{
		quit_error((char*)"param_list() is empty");
	}
	
	return this->param_list->get_double_param(name, defaultValue);
}

//****************************************************************************
//get_int_param()
//****************************************************************************
double Param_t::get_double_param(char *name, double defaultValue)
{
	Param_t *p;
	
	p = this->get_param(name);
	
	if(p == NULL)
	{
		fprintf(stdout, "No value supplied, so using default value for %s = %f\n", name, defaultValue);
		return defaultValue;
	}
	else
		return p->get_double_param(name);
}

//****************************************************************************
//get_string_param()
//****************************************************************************
char * ParamListHead_t::get_string_param(char *name)
{
	if(this->param_list == NULL)
	{
		quit_error((char*)"param_list() is empty");
	}
	
	return this->param_list->get_string_param(name);
}

//****************************************************************************
//get_param_line()
//****************************************************************************
char * ParamListHead_t::get_param_line(char *name)
{
	if(this->param_list == NULL)
	{
		quit_error((char*)"param_list() is empty");
	}
	
	return this->param_list->get_param_line(name);
}

//****************************************************************************
//get_string_param()
//****************************************************************************
char * Param_t::get_string_param(char *name)
{
	Param_t *p = this;
	char *str;
	
	while( p != NULL )
	{
		if(strcmp(p->name,name)==0)
		{
			p->numberOfAccesses++;
			str = copy_string(p->line);
			return str;
		}
		else
		{
			p = p->next;
		}
	}
	fprintf(stderr,"In get_double_param, can't find param \"%s\"\n", name);
	quit_error((char*)"Unable to find param \"%s\"", name);

	return NULL;
}

// ****************************************************************************
// get_param_line()
// ****************************************************************************
char * Param_t::get_param_line(char *name)
{
	Param_t *p = this;
	
	while( p != NULL )
	{
		if(strcmp(p->name,name)==0)
		{
			p->numberOfAccesses++;
			return p->line;
		}
		else
		{
			p = p->next;
		}
	}
	fprintf(stderr,"In get_param_line, can't find param \"%s\"\n", name);
	quit_error((char*)"Unable to find param \"%s\"", name);
	
	return NULL;
}

//****************************************************************************
//get_param()
// This starts at the head of the list it is supplied with and searches for
// a param with the specified name. If it can't find such a param, then NULL
// is returned.
//****************************************************************************
Param_t * Param_t::get_param(char *name)
{
	Param_t *param_list = this;
	
	while( param_list != NULL )
	{
		if(strcmp(param_list->name,name)==0)
		{
			param_list->numberOfAccesses++;			
			return param_list;
		}
		else
		{
			param_list = param_list->next;
		}
	}

	return NULL;
}




//****************************************************************************
//get_param()
// This starts at the head of the list it is supplied with and searches for
// a param with the specified name. If it can't find such a param, then NULL
// is returned.
//****************************************************************************
Param_t * ParamListHead_t::get_param(char *name)
{
	if(this->param_list == NULL)
	{
		quit_error((char*)"param_list() is empty");
	}
	
	return this->param_list->get_param(name);
}

//****************************************************************************
//get_next_param()
// This skips the head of the list it is supplied with and searches the rest for
// a param with the specified name. If it can't find such a param, then NULL
// is returned.
//****************************************************************************
Param_t * Param_t::get_next_param(char *name)
{
	Param_t *param_list = this;
	
	if(param_list == NULL)
		return NULL;
	
	// Skip the head
	param_list = param_list->next;
	
	while( param_list != NULL )
	{
		if(strcmp(param_list->name,name)==0)
		{
			param_list->numberOfAccesses++;			
			return param_list;
		}
		else
		{
			param_list = param_list->next;
		}
	}
	
	return NULL;
}

//****************************************************************************
//get_next_param()
// This skips the head of the list it is supplied with and searches the rest for
// a param with the specified name. If it can't find such a param, then NULL
// is returned.
//****************************************************************************
Param_t * ParamListHead_t::get_next_param(char *name)
{
	Param_t *param_list = this->param_list;
	
	if(param_list == NULL)
		return NULL;

	// Skip the head
	param_list = param_list->next;

	while( param_list != NULL )
	{
		if(strcmp(param_list->name,name)==0)
		{
			param_list->numberOfAccesses++;			
			return param_list;
		}
		else
		{
			param_list = param_list->next;
		}
	}
	
	return NULL;
}

//****************************************************************************
//ParamList::ParamList()
//****************************************************************************
ParamListHead_t::ParamListHead_t()
{
	// intitialise
	this->param_list = NULL;
	this->num_params = 0;
}


