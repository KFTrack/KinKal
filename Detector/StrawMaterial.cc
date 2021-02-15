#include "KinKal/Detector/StrawMaterial.hh"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace KinKal {
  double StrawMaterial::gasPath(ClosestApproachData const& cadata,StrawXingConfig const& caconfig) const {
    double doca = std::min(fabs(cadata.doca()),srad_);
    double sigdoca = sqrt(cadata.docaVar());
    double adot = cadata.dirDot();
    doca = std::min(fabs(doca),srad_);
    double afac = 1.0/sqrt(1.0-adot*adot); // angle factor, =1.0/sin(theta)
    if(!isfinite(afac))throw std::runtime_error("Invalid angle");

    double retval;
    if (sigdoca < caconfig.minsigdoca_) {  // small error: don't integrate
      double rad = std::min(doca,srad_-thick_);
      retval = 2.0*sqrt(srad2_-rad*rad); 
    } else { 
      // integrate +- N sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rmin is handled naturally
      double rmax = std::min(srad_-thick_,(doca+caconfig.nsig_*sigdoca))/srad_;
      double rmin = std::min(std::max(wrad_,doca-caconfig.nsig_*sigdoca)/srad_,rmax-thick_);
      retval = srad_*(asin(rmax) - asin(rmin) +
	  rmax*sqrt(1.0-rmax*rmax) - rmin*sqrt(1.0-rmin*rmin) )/(rmax-rmin);
      if(isnan(retval))throw std::runtime_error("Invalid pathlength");
    }
    // correct for the angle
    retval *= afac;
    return retval;
  }

  double StrawMaterial::wallPath(ClosestApproachData const& cadata,StrawXingConfig const& caconfig) const{
    double doca = std::min(fabs(cadata.doca()),srad_);
    double sigdoca = sqrt(cadata.docaVar());
    double adot = cadata.dirDot();
    doca = std::min(fabs(doca),srad_);
    double afac = 1.0/sqrt(1.0-adot*adot); // angle factor, =1.0/sin(theta)
    double retval(-1.0);
    if (sigdoca < caconfig.minsigdoca_) {  // small error: don't integrate
      double rad = std::min(doca,srad_-thick_);
      retval = 2.0*thick_*srad_/sqrt(srad2_-rad*rad);
    } else {
      // integrate +- N sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rmin is handled naturally
      double rmax = std::min(srad_-thick_,(doca+caconfig.nsig_*sigdoca))/srad_;
      double rmin = std::min(std::max(wrad_,doca-caconfig.nsig_*sigdoca)/srad_,rmax-thick_);
      retval = 2.0*thick_*(asin(rmax) - asin(rmin))/(rmax-rmin);
      if(isnan(retval))throw std::runtime_error("Invalid pathlength");
    }
    // correct for the angle
    retval *= afac;
    return retval;
  }

  double StrawMaterial::wirePath(ClosestApproachData const& cadata,StrawXingConfig const& caconfig) const{
  // Model the wire as a diffuse gas, spread given by the DOCA sigma.  This gives a probabilistic interpretation of the material
    return -1.0; //TODO
  }

  void StrawMaterial::findXings(ClosestApproachData const& cadata,StrawXingConfig const& caconfig, std::vector<MaterialXing>& mxings) const {
    mxings.clear();
    double wpath = wallPath(cadata,caconfig);
    if(wpath > 0.0) mxings.push_back(MaterialXing(*wallmat_,wpath));
    double gpath = gasPath(cadata,caconfig);
    if(gpath > 0.0) mxings.push_back(MaterialXing(*gasmat_,gpath));
    double wirepath = wirePath(cadata,caconfig);
    if(wirepath > 0.0) mxings.push_back(MaterialXing(*wiremat_,wirepath));
  }

}
