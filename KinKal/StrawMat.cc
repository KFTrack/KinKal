#include "KinKal/StrawMat.hh"
#include <cmath>
#include <algorithm>
namespace KinKal {
  float StrawMat::gasPath(float doca, float ddoca, float adot) const {
    float retval;
  // if the uncertainty is large, just take the average over all possible impact parameters
    if(ddoca > rad_) {
      retval = 0.5*M_PI*rad_;
    } else if (doca-ddoca > rad_) {
    // no intersection within errors!
      retval = 0.0;
    } else {
      // test points +- 1 sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rdmin is handled naturally
      float rdmax = std::min(rdmax_,(doca+ddoca)/rad_);
      float rdmin = (doca-ddoca)/rad_;
      retval = rad2_*(asin(rdmax) - asin(rdmin) +
	  rdmax*sqrt(1.0-rdmax*rdmax) - rdmin*sqrt(1.0-rdmin*rdmin) )/(rdmax-rdmin);
    }
    // correct for the angle
    return retval/sqrt(1.0-adot*adot);
  }

  float StrawMat::wallPath(float doca, float ddoca, float adot) const{
    float retval;
    // if the uncertainty is large, just take the average over all possible impact parameters
    if(ddoca > rad_) {
      retval = M_PI*thick_;
    } else if (doca-ddoca > rad_) {
    // no intersection within errors!
      retval = 0.0;
    } else {
      // test points +- 1 sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rdmin is handled naturally
      float rdmax = std::min(rdmax_,(doca+ddoca)/rad_);
      float rdmin = (doca-ddoca)/rad_;
      retval = 2.0*rad2_*thick_*(asin(rdmax) - asin(rdmin))/(rdmax-rdmin);
    }
    // correct for the angle
    return retval/sqrt(1.0-adot*adot);
  }

  void StrawMat::intersect(TPocaBase const& tpoca, std::vector<MIsect>& misects) const {
  // calculate the angle between traj0 and traj1
    Vec3 dir0, dir1;
    tpoca.ttraj0().direction(tpoca.t0(),dir0);
    tpoca.ttraj1().direction(tpoca.t1(),dir1);
    double adot = dir0.Dot(dir1);
    return intersect(tpoca.doca(),tpoca.dDoca(),adot,misects);
  }

  void StrawMat::intersect(float doca, float ddoca, float adot, std::vector<MIsect>& misects) const {
    misects.clear();
    misects.push_back(MIsect(wallmat_,wallPath(doca,ddoca,adot)));
    misects.push_back(MIsect(gasmat_,gasPath(doca,ddoca,adot)));
// for now, take 0 path on the wire: this should be some probability based on doca and ddoca FIXME!
    if(rwire_<0.0) misects.push_back(MIsect(wiremat_,0.0));
  }

}
