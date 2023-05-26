#include <string>
#include <time.h>
#include <sstream>
#include "TimeFunctions.h"

using namespace::std;
std::string timeMonthDayTime() {
    time_t rawTime;
    char timeChar[100];
    time(&rawTime);
    strftime(timeChar,80,"%b %d %H:%M:%SS",localtime(&rawTime));
    std::string timeString=timeChar;
    timeString.erase(timeString.end()-1,timeString.end());
    return timeString;
};

std::string timeMonthDayTime(time_t &rawTime) {
    char timeChar[100];
    strftime(timeChar,80,"%b %d %H:%M:%SS",localtime(&rawTime));
    std::string timeString=timeChar;
    timeString.erase(timeString.end()-1,timeString.end());
    return timeString;
}

string get_local_time(){
    time_t a;
    time(&a);
    struct tm* b=localtime(&a);
    int cur_year=b->tm_year+1900;
    int cur_mon=b->tm_mon+1;
    int cur_day=b->tm_mday;
    ostringstream out_time;
    out_time<<cur_year<<"-"<<cur_mon<<"-"<<cur_day<<"  "<<b->tm_hour<<":"<<b->tm_min<<":"<<b->tm_sec<<endl;
    string out=out_time.str();
    out.erase(out.end()-1);
    return out;
}
