#ifndef TIME_FUNCTIONS_DEF
#define TIME_FUNCTIONS_DEF
#include <string>
#include <time.h>
using namespace::std;
string timeMonthDayTime();
string timeMonthDayTime(time_t &rawTime);
string get_local_time();
#endif
