// time_stats.cpp
// ****************************************************************************
#include <cstdio>
// ****************************************************************************
#include "time_stats.h"
// ****************************************************************************

// ****************************************************************************
// TimeStats::TimeStats()
// ****************************************************************************
TimeStats::TimeStats(long int step_start, long int step_finish)
{
	reset(step_start, step_finish);
	return;
}

// ****************************************************************************
// TimeStats::TimeStats()
// ****************************************************************************
void TimeStats::reset(long int step_start, long int step_finish)
{
	this->step_start = step_start;
	this->step_finish = step_finish;
	total_steps = step_finish - step_start;	
	time(&(time_start));
	
	return;
}

// ****************************************************************************
// TimeStats::convert_seconds()
// ****************************************************************************
inline void TimeStats::convert_seconds(double total_seconds, int &days, int &hours, int &minutes, int &seconds)
{
	days = (int) (total_seconds / (60.0*60.0*24.0));
	total_seconds -= (double)(days * 24 * 60 * 60);
	hours = (int) (total_seconds / (60.0*60.0));
	total_seconds -= (double)(hours * 60 * 60);
	minutes = (int) (total_seconds / (60.0));
	total_seconds -= (double)(minutes * 60);
	seconds = (int) total_seconds;
	
	return;
}

// ****************************************************************************
// TimeStats::produce_stats()
// ****************************************************************************
void TimeStats::produce_stats(long int step_current)
{
	double steps_done = (double) (step_current - step_start);
	double steps_remaining = (double) (step_finish - step_current);
	int days[2], hours[2], minutes[2], seconds[2];
	double time_spent;
	double time_remaining;
	time_t time_current;

	// *********************
	// grab the current time
	// *********************
	time(&(time_current));

	time_spent = difftime(time_current, time_start);
	
	percentage_complete = steps_done / total_steps * 100.0;
	
	convert_seconds(time_spent, days[0], hours[0], minutes[0], seconds[0]);
	
	
	if( steps_done > 0.0 )
		time_remaining = time_spent / steps_done * steps_remaining;
	else
		time_remaining = 0.0;
	
	convert_seconds(time_remaining, days[1], hours[1], minutes[1], seconds[1]);
	
	if(time_spent > 0.0)
		timesteps_per_minute = steps_done / time_spent * 60.0;		
	else
		timesteps_per_minute = 0.0;
	
	if(days[0] == 0 && days[1] == 0)
	{
		sprintf(time_spent_string,"%02dh%02dm%02ds", hours[0], minutes[0], seconds[0]);
		sprintf(time_remaining_string,"%02dh%02dm%02ds", hours[1], minutes[1], seconds[1]);
	}
	else
	{
		sprintf(time_spent_string,"%dd%02dh%02dm%02ds", days[0], hours[0], minutes[0], seconds[0]);
		sprintf(time_remaining_string,"%dd%02dh%02dm%02ds", days[1], hours[1], minutes[1], seconds[1]);
	}

	return;
}



