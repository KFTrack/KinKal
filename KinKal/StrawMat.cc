#include "KinKal/StrawMat.hh"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace KinKal {
  float StrawMat::gasPath(float doca, float ddoca, float adot) const {
    doca = fabs(doca);
    double afac = 1.0/sqrt(1.0-adot*adot); // angle factor, =1.0/sin(theta)
    float retval;
  // if the uncertainty is large, just take the average over all possible impact parameters
    if(ddoca > rad_) {
      retval = 0.5*M_PI*rad_;
    } else if (ddoca < ddmax_) {  // small error: don't integrate
      if(doca < rad_)
	retval = 2.0*sqrt(rad2_-doca*doca); 
      else
	retval = 0.0;
    } else if ( doca - rad_ > ddoca ) {
    // near the edge, set to minimal value
      retval = 0.0;
    } else { 
      // integrate +- 1 sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rdmin is handled naturally
      float rdmax = std::min(rdmax_,(doca+ddoca)/rad_);
      float rdmin = (doca-ddoca)/rad_;
      retval = rad_*(asin(rdmax) - asin(rdmin) +
	  rdmax*sqrt(1.0-rdmax*rdmax) - rdmin*sqrt(1.0-rdmin*rdmin) )/(rdmax-rdmin);
    }
    // correct for the angle
    retval *= afac;
    if(isnan(retval))throw std::runtime_error("Invalid pathlength");
    return retval;
  }

  float StrawMat::wallPath(float doca, float ddoca, float adot) const{
    doca = fabs(doca);
    double afac = 1.0/sqrt(1.0-adot*adot); // angle factor, =1.0/sin(theta)
    float retval(-1.0);
    // if the uncertainty is large, just take the average over all possible impact parameters
    if(ddoca > rad_ ) {
      retval = M_PI*thick_;
    } else if (ddoca < ddmax_) {  // small error: don't integrate
      if(rad_ > doca + 0.25*thick_)
	retval = 2.0*thick_*rad_/sqrt(rad2_-doca*doca);
      else
	retval = wpmax_;
    } else if ( doca - rad_ > ddoca ) {
    // way outside: set to a maximal value
      retval = wpmax_;
    } else {
      // integrate +- 1 sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rdmin is handled naturally
      float rdmax = std::min(rdmax_,(doca+ddoca)/rad_);
      float rdmin = std::max(-rdmax_,(doca-ddoca)/rad_);
      retval = 2.0*thick_*(asin(rdmax) - asin(rdmin))/(rdmax-rdmin);
    }
    // correct for the angle
    retval *= afac;
    if(isnan(retval))throw std::runtime_error("Invalid pathlength");
    return retval;
  }

  void StrawMat::findXings(TPocaBase const& tpoca, std::vector<MatXing>& mxings) const {
  // calculate the angle between traj0 and traj1
    return findXings(tpoca.doca(),tpoca.dDoca(),tpoca.dirDot(),mxings);
  }

  void StrawMat::findXings(float doca, float ddoca, float adot, std::vector<MatXing>& mxings) const {
    mxings.clear();
    mxings.push_back(MatXing(wallmat_,wallPath(doca,ddoca,adot)));
    mxings.push_back(MatXing(gasmat_,gasPath(doca,ddoca,adot)));
// for now, take 0 path on the wire: this should be some probability based on doca and ddoca FIXME!
    if(rwire_<0.0) mxings.push_back(MatXing(wiremat_,0.0));
  }

}
