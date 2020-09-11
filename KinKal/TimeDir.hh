#ifndef KinKal_TimeDir_hh
#define KinKal_TimeDir_hh
//
//  Define directions of time
//
#include <ostream>
namespace KinKal { 
  enum struct TimeDir{forwards=0,backwards, end}; // time direction; forwards = increasing time, backwards = decreasing time
  TimeDir& operator ++ (TimeDir& tdir);
  std::ostream& operator <<(std::ostream& ost, TimeDir const& tdir);
}
#endif
