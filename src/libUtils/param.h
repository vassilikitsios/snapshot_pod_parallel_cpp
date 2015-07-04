//param.h

#define MAX_PARAM_LINE_LENGTH 256

//typedef class Param_t;
//typedef class ParamListHead_t;

//****************************************************************************
//Param
//****************************************************************************
class Param_t
{
private:
	
public:
	//-------------------------------------------------------------------------
	//variables
	//-------------------------------------------------------------------------
	char *name;
	char *line;
	int numberOfAccesses;
	Param_t *next;
	Param_t *prev;
	//-------------------------------------------------------------------------
	//functions
	//-------------------------------------------------------------------------
	//functions
	Param_t(char *name, char *line); // constructor
	//-------------------------------------------------------------------------
	int get_int_param(char *name);
	int get_int_param(char *name, int defaultValue);
	//-------------------------------------------------------------------------
	int get_long_int_param(char *name);
	//-------------------------------------------------------------------------
	double get_double_param(char *name);
	//-------------------------------------------------------------------------
	char * get_string_param(char *name);
	char * get_param_line(char *name);
	//-------------------------------------------------------------------------
	double get_double_param(char *name, double defaultValue);
	//-------------------------------------------------------------------------
	Param_t * get_param(char *name);
	Param_t * get_next_param(char *name);
};

//****************************************************************************
//Param
//****************************************************************************
class ParamListHead_t
{
private:
	
public:
	//-------------------------------------------------------------------------
	//variables
	//-------------------------------------------------------------------------
	Param_t *param_list;
	int num_params;
	//-------------------------------------------------------------------------
	//functions
	//-------------------------------------------------------------------------
	ParamListHead_t(); //contructor
	//-------------------------------------------------------------------------
	void addParam(char *param_line);
	//-------------------------------------------------------------------------
	void print_params();
	void print_param_accesses();
	//-------------------------------------------------------------------------
	void get_params_from_file(char *filename);
	void get_params_from_file(char *filename, bool verbose);
	void get_params_from_command_line(int argc, char *argv[]);
	//-------------------------------------------------------------------------
	bool check_param(char *name);
	//-------------------------------------------------------------------------
	int get_int_param(char *name);
	int get_int_param(char *name, int defaultValue);
	//-------------------------------------------------------------------------
	int get_long_int_param(char *name);
	//-------------------------------------------------------------------------
	double get_double_param(char *name);
	double get_double_param(char *name, double defaultValue);
	//-------------------------------------------------------------------------
	Param_t * get_param(char *name);
	Param_t * get_next_param(char *name);
	char * get_string_param(char *name);
	char * get_param_line(char *name);
	//-------------------------------------------------------------------------
	Param_t * get_nth_param(int n, char *name);
	int get_nth_int_param(int n, char *name);
	long int get_nth_long_int_param(int n, char *name);
	double get_nth_double_param(int n, char *name);
};


