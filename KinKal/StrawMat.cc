#include "KinKal/StrawMat.hh"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace KinKal {
  double StrawMat::gasPath(double doca, double ddoca, double adot) const {
    doca = std::min(fabs(doca),srad_);
    double afac = 1.0/sqrt(1.0-adot*adot); // angle factor, =1.0/sin(theta)
    if(!isfinite(afac))throw std::runtime_error("Invalid angle");

    double retval;
    if (ddoca < ddmax_) {  // small error: don't integrate
      double rad = std::min(doca,srad_-thick_);
      retval = 2.0*sqrt(srad2_-rad*rad); 
    } else { 
      // integrate +- 1 sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rdmin is handled naturally
      double rdmax = std::min(rdmax_,(doca+ddoca)/srad_);
      double rdmin = std::min(std::max(wrad_,doca-ddoca)/srad_,rdmax-thick_);
      retval = srad_*(asin(rdmax) - asin(rdmin) +
	  rdmax*sqrt(1.0-rdmax*rdmax) - rdmin*sqrt(1.0-rdmin*rdmin) )/(rdmax-rdmin);
      if(isnan(retval))throw std::runtime_error("Invalid pathlength");
    }
    // correct for the angle
    retval *= afac;
    return retval;
  }

  double StrawMat::wallPath(double doca, double ddoca, double adot) const{
    doca = std::min(fabs(doca),srad_);
    double afac = 1.0/sqrt(1.0-adot*adot); // angle factor, =1.0/sin(theta)
    double retval(-1.0);
    if (ddoca < ddmax_) {  // small error: don't integrate
      double rad = std::min(doca,srad_-thick_);
      retval = 2.0*thick_*srad_/sqrt(srad2_-rad*rad);
    } else {
      // integrate +- 1 sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rdmin is handled naturally
      double rdmax = std::min(rdmax_,(doca+ddoca)/srad_);
      double rdmin = std::min(std::max(wrad_,doca-ddoca)/srad_,rdmax-thick_);
      retval = 2.0*thick_*(asin(rdmax) - asin(rdmin))/(rdmax-rdmin);
      if(isnan(retval))throw std::runtime_error("Invalid pathlength");
    }
    // correct for the angle
    retval *= afac;
    return retval;
  }

  void StrawMat::findXings(double doca, double ddoca, double adot, std::vector<MatXing>& mxings) const {
    mxings.clear();
    double wpath = wallPath(doca,ddoca,adot);
    if(wpath > 0.0) mxings.push_back(MatXing(*wallmat_,wpath));
    double gpath = gasPath(doca,ddoca,adot);
    if(gpath > 0.0) mxings.push_back(MatXing(*gasmat_,gpath));
// for now, ignore the wire: this should be based on the probability that the wire was hit given doca and ddoca FIXME!
    if(wrad_<0.0) mxings.push_back(MatXing(*wiremat_,0.0));
  }

}
