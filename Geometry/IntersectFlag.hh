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
    bool inrange_ = false;  // intersection is inside the time range
    bool hascov_ = false; // the projected spatial covariance (Intersection::pcov_) is filled. This is optional, and independent of good()
    bool operator ==(IntersectFlag const& other) const { return other.onsurface_ == onsurface_ && other.inbounds_ == inbounds_ && other.inrange_ == inrange_ && other.hascov_ == hascov_;}
    bool operator !=(IntersectFlag const& other) const { return other.onsurface_ != onsurface_ || other.inbounds_ != inbounds_ || other.inrange_ != inrange_ || other.hascov_ != hascov_;}
    bool good() const { return onsurface_ && inbounds_ && inrange_; }
    friend std::ostream& operator <<(std::ostream& ost, KinKal::IntersectFlag const& iflag);
  };
}
#endif
