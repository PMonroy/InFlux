#ifndef EQDAYMONTHDATE
#define EQDAYMONTHDATE

using namespace std;

class eqdate{

  unsigned int year,month,mday,day;   
public:
  
  eqdate(void);                          // Zero date Constructor
  eqdate(unsigned int y,unsigned int m,unsigned int d); // Constructor
  eqdate(const eqdate &d);            // Copy Vector Constructor
  
  unsigned int GetYear(void ) const; //Getter
  unsigned int GetMonth(void ) const; //Getter
  unsigned int GetMday(void ) const; //Getter
  unsigned int GetDay(void ) const; //Getter


  void SetDay(unsigned int d); //setter
  void SetDate(unsigned int y,unsigned int m,unsigned int md); //setter days

  ~eqdate();	                            // Destructor
};

#endif
