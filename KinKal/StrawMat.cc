#include "KinKal/StrawMat.hh"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace KinKal {
  float StrawMat::gasPath(float doca, float ddoca, float adot) const {
    doca = fabs(doca);
    double afac = 1.0/sqrt(1.0-adot*adot); // angle factor, =1.0/sin(theta)
    if(!isfinite(afac))throw std::runtime_error("Invalid angle");

    float retval;
  // if the uncertainty is large, just take the average over all possible impact parameters
    if(ddoca > srad_) {
      retval = 0.5*M_PI*srad_;
    } else if (ddoca < ddmax_) {  // small error: don't integrate
      if(doca < srad_)
	retval = 2.0*sqrt(srad2_-doca*doca); 
      else
	retval = 0.0;
    } else if ( doca - srad_ > ddoca ) {
    // near the edge, set to minimal value
      retval = 0.0;
    } else { 
      // integrate +- 1 sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rdmin is handled naturally
      float rdmax = std::min(rdmax_,(doca+ddoca)/srad_);
      float rdmin = (doca-ddoca)/srad_;
      retval = srad_*(asin(rdmax) - asin(rdmin) +
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
    if(ddoca > srad_ ) {
      retval = M_PI*thick_;
    } else if (ddoca < ddmax_) {  // small error: don't integrate
      if(srad_ > doca + 0.25*thick_)
	retval = 2.0*thick_*srad_/sqrt(srad2_-doca*doca);
      else
	retval = wpmax_;
    } else if ( doca - srad_ > ddoca ) {
    // way outside: set to a maximal value
      retval = wpmax_;
    } else {
      // integrate +- 1 sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rdmin is handled naturally
      float rdmax = std::min(rdmax_,(doca+ddoca)/srad_);
      float rdmin = std::max(-rdmax_,(doca-ddoca)/srad_);
      retval = 2.0*thick_*(asin(rdmax) - asin(rdmin))/(rdmax-rdmin);
    }
    // correct for the angle
    retval *= afac;
    if(isnan(retval))throw std::runtime_error("Invalid pathlength");
    return retval;
  }

  void StrawMat::findXings(float doca, float ddoca, float adot, std::vector<MatXing>& mxings) const {
    mxings.clear();
    float wpath = wallPath(doca,ddoca,adot);
    if(wpath > 0.0) mxings.push_back(MatXing(wallmat_,wpath));
    float gpath = gasPath(doca,ddoca,adot);
    if(gpath > 0.0) mxings.push_back(MatXing(gasmat_,gpath));
// for now, ignore the wire: this should be based on the probability that the wire was hit given doca and ddoca FIXME!
    if(wrad_<0.0) mxings.push_back(MatXing(wiremat_,0.0));
  }

}
