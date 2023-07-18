//
//  Data payload for the intersection point of a PiecewiseTrajectory with a surface
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKal_InterData_hh
#define KinKal_InterData_hh
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/IntersectFlag.hh"

namespace KinKal {
  struct InterData {
    InterData() : time_(0.0) {}
    IntersectFlag flag_; // intersection status
    VEC3 pos_; // intersection position
    VEC3 norm_; // surface normal at intersection
    VEC3 pdir_; // particle direction at intersection
    double time_; // time at intersection (from particle)
  };
}
#endif
