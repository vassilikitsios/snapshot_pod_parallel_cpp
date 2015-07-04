// time_stats.h
// ****************************************************************************
#include <ctime>
// ****************************************************************************
//typedef class TimeStats;
// ****************************************************************************
// TimeStats
// ****************************************************************************
class TimeStats
{
private:
	time_t time_start;
	long int computationStartStep;
	long int step_start;
	long int step_finish;
	long int total_steps;
	inline void convert_seconds(double total_seconds, int &days, int &hours, int &minutes, int &seconds);
public:
	char time_spent_string[40];
	char time_remaining_string[40];
	double percentage_complete;
	double timesteps_per_minute;
		
	TimeStats(long int step_start, long int step_finish);
	void reset(long int step_start, long int step_finish);
	void produce_stats(long int step_current);
};
// ****************************************************************************




