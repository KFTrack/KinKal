#ifndef KinKal_IntersectFlag_hh
#define KinKal_IntersectFlag_hh
//
// describe flag bits used for reco intersection
//
//
// Original author David Brown
//
#include <string>
#include <map>
namespace KinKal {
  struct IntersectFlag {
    IntersectFlag() : onsurface_(false), inbounds_(false) {}
    bool onsurface_; // intersection is on the surface
    bool inbounds_;  // intersection is inside the surface boundaries
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::IntersectFlag const& iflag);
#endif