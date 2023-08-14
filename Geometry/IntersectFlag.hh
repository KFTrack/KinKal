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
    bool onsurface_ = false; // intersection is on the surface
    bool inbounds_ = false;  // intersection is inside the surface boundaries
    bool operator ==(IntersectFlag const& other) const { return other.onsurface_ == onsurface_ && other.inbounds_ == inbounds_; }
    bool operator !=(IntersectFlag const& other) const { return other.onsurface_ != onsurface_ || other.inbounds_ != inbounds_; }
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::IntersectFlag const& iflag);
#endif
