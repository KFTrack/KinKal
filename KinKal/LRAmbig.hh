#ifndef KinKal_LRAmbig_hh
#define KinKal_LRAmbig_hh
//
// define left-right ambiguity, = sign of angular momentum of particle trajectory WRT wire
// Null means to ignore drift information and constrain to the wire
//
#include <ostream>
namespace KinKal { 
  enum struct LRAmbig { left=-1, null=0, right=1}; 
  std::ostream& operator <<(std::ostream& ost, LRAmbig const& tdir);
}
#endif

