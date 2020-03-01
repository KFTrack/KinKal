#ifndef KinKal_TDMInter_hh
#define KinKal_TDMInter_hh
//
//  Describe the intersection of a trajectory with a piece of material
//  Used in the kinematic Kalman fit
//
#include "KinKal/TTraj.hh"
#include "KinKal/DMat.hh"
namespace KinKal {
  class TDMInter {
    public:
      TMat(double tinter, TTraj const& ttraj, DMat const& dmat) : tinter_(tinter),
      ttraj_(ttraj), dmat_(dmat) {}
    private:
      double tinter_; // time of the intersection
      TTraj const& ttraj_; // trajectory which intersects the material
      DMat const& dmat_; // detector matter intersected
  };
}
#endif
