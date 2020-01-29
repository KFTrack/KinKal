#ifndef KinKal_POCAFinder_hh
#define KinKal_POCAFinder_hh
///
//  This functor class finds the (spacetime) points of closest approach between two TTrajs.
//  Concrete instances are specializations and must be implemented explicity for
//  each trajectory type pair.
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/Types.hh"
#include <iostream>

namespace KinKal {

  struct TPOCA {
    Vec4 p1_, p2_; // spacetime points at GEOMETRIC closest approach
    void delta(Vec4& ds) const { ds = p2_-p1_; }
    double dt() const { return p2_.T() - p1_.T(); }
    void dVec(Vec3& ds) const { Vec3 p1(p1_); Vec3 p2(p2_); ds = p2-p1; }
  };

  template<class T1, class T2> void FindTPOCA(T1 const& traj1, T2 const& traj2, TPOCA& tpoca) {
    // no default implementation so throw: this must be specialized to be useful
    std::cout << "NoFindPOCA implementation for types " << typeid(T1).name()
    << " and " << typeid(T2).name() << std::endl; // should throw FIXME!
  }

}

#endif
