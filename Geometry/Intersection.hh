//
//  Complete payload for intersection calculation
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKal_Intersection_hh
#define KinKal_Intersection_hh
#include "KinKal/Geometry/InterData.hh"
#include "KinKal/Geometry/Surface.hh"

namespace KinKal {
  // intersection product
  template <class KTRAJ> struct Intersection : public InterData {
    Intersection(KTRAJ const& ktraj, Surface const& surf,TimeRange const& trange, double tol) : InterData(trange), ktraj_(ktraj), surf_(surf), tol_(tol) {}
    KTRAJ const& ktraj_; // trajectory of this intersection
    Surface const& surf_; // surf of this intersection
    double tol_; // tol used in this intersection
  };
}
#endif
