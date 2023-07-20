//
//  Generic payload for the intersection of a Trajectory with a surface
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKal_InterData_hh
#define KinKal_InterData_hh
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/IntersectFlag.hh"
#include "KinKal/Geometry/Ray.hh"

namespace KinKal {
  struct InterData {
    InterData() : time_(0.0) {}
    InterData(TimeRange const& trange) : time_(0.0), trange_(trange) {}
    IntersectFlag flag_; // intersection status
    VEC3 pos_; // intersection position
    VEC3 norm_; // surface normal at intersection
    VEC3 pdir_; // particle direction at intersection
    double time_; // time at intersection (from particle)
    TimeRange trange_; // time range used to search for this intersection
    bool inRange() const { return trange_.inRange(time_);}
    Ray ray() const { return Ray(pdir_,pos_); }
  };
}
#endif
