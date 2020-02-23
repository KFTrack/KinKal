#ifndef KinKal_TDir_hh
#define KinKal_TDir_hh
//
//  Define processing directions according to the flow of time
//
#include <ostream>
namespace KinKal { 
  enum struct TDir{forwards=0,backwards, end}; // time direction; forwards = increasing time, backwards = decreasing time
  TDir& operator ++ (TDir& tdir);
  std::ostream& operator <<(std::ostream& ost, TDir const& tdir);
}
#endif
