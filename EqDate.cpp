#include "EqDate.hpp"

eqdate::eqdate(void) {  // Void Constructor
  year=0;
  month=1;
  mday=1;
  day=0;
}

eqdate::eqdate(unsigned int y,unsigned int m,unsigned int md) { // Constructor
  year=y;
  month=m;
  mday=md;
  day=y*360+(m-1)*30+(mday-1);
}

eqdate::eqdate(const eqdate &date){   // Copy Vector Constructor
  year=date.year;
  month=date.month;
  mday=date.mday;
  day=date.day;
}

unsigned int eqdate::GetYear(void) const { // Getter
  return year;
}
unsigned int eqdate::GetMonth(void) const { // Getter
  return month; 
}
unsigned int eqdate::GetMday(void) const { // Getter
  return mday;
}
unsigned int eqdate::GetDay(void) const { // Getter
  return day;
}

void eqdate::SetDay(unsigned int d) { // Setter
  day = d;
  year=d/360;
  month=(d%360)/30+1;
  mday=(d%360)%30+1;
}

void eqdate::SetDate(unsigned int y,unsigned int m,unsigned int md) { // Setter
  year=y;
  month=m;
  mday=md;
  day=y*360+(m-1)*30+(mday-1);
}

eqdate::~eqdate(){}    // Destructor
