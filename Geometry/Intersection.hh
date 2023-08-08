//
//  payload for intersection calculation
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKal_Intersection_hh
#define KinKal_Intersection_hh
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/IntersectFlag.hh"
#include "KinKal/Geometry/Ray.hh"
#include <memory>

namespace KinKal {
  struct Intersection {
    Intersection() : time_(0.0) {}
    IntersectFlag flag_; // intersection status
    VEC3 pos_; // intersection position
    VEC3 norm_; // surface normal at intersection
    VEC3 pdir_; // particle direction at intersection
    TimeRange range_; // time range searched for this intersection
    double time_; // time at intersection (from particle)
    // simple utility functions
    bool inRange(TimeRange const& trange) const { return flag_.onsurface_ && trange.inRange(time_); }
    Ray ray() const { return Ray(pdir_,pos_); }
  };
}
#endif
