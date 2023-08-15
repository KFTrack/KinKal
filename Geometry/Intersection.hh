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
#include <ostream>

namespace KinKal {
  struct Intersection : public IntersectFlag {
    VEC3 pos_; // intersection position
    VEC3 norm_; // surface normal at intersection
    VEC3 pdir_; // particle direction at intersection
    double time_ = 0.0; // time at intersection (from particle)
    bool gap_ = false; // intersection is in a piecewise-trajectory gap
    Ray ray() const { return Ray(pdir_,pos_); }
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Intersection const& inter);
#endif
