#include "KinKal/Detector/StrawMaterial.hh"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace KinKal {
  void StrawMaterial::pathLengths(ClosestApproachData const& cadata,StrawXingConfig const& caconfig,
      double& wallpath, double& gaspath, double& wirepath) const {
    wallpath = gaspath = wirepath = 0.0;
    double doca = std::min(fabs(cadata.doca()),srad_-thick_);
    double sigdoca = sqrt(cadata.docaVar());
    if (sigdoca < caconfig.minsigdoca_) {  // small error: don't integrate
      double ddoca = doca*doca;
      gaspath = 2.0*sqrt(srad2_-ddoca);
      wallpath = 2.0*sqrt(srad2_ + 2*thick_*srad_ - ddoca) - gaspath;
   } else if (sigdoca*caconfig.nsig_ > srad_ ) {
      // errors are large WRT the size of the straw, or DOCA is very far from the wire: just take the average over all impact parameters
      wallpath = M_PI*thick_;
      gaspath = 0.5*M_PI*srad_;
    } else {
      // integrate +- N sigma from DOCA.  Restrict rmax to physical values
      // Note that negative rmin is handled naturally
      double rmax = std::min(srad_-thick_,(doca+caconfig.nsig_*sigdoca))/srad_;
      double rmin = std::max(-srad_+thick_,doca-caconfig.nsig_*sigdoca)/srad_;
      wallpath = 2.0*thick_*(asin(rmax) - asin(rmin))/(rmax-rmin);
      gaspath = srad_*(asin(rmax) - asin(rmin) +
          rmax*sqrt(1.0-rmax*rmax) - rmin*sqrt(1.0-rmin*rmin) )/(rmax-rmin);
    }
    if(isnan(wallpath) || isnan(gaspath))throw std::runtime_error("Invalid pathlength");
    // Model the wire as a diffuse gas, density constrained by DOCA TODO
    // correct for the angle WRT the axis
    double afac = angleFactor(cadata.dirDot());
    wallpath *= afac;
    gaspath *= afac;
  }

  double StrawMaterial::transitLength(ClosestApproachData const& cadata) const {
    double doca = std::min(fabs(cadata.doca()),srad_-thick_);
    double tlen = 2.0*sqrt(srad2_-doca*doca);
    // correct for the angle WRT the axis
    double afac = angleFactor(cadata.dirDot());
    tlen *= afac;
    return tlen;
  }

  double StrawMaterial::angleFactor(double dirdot) const {
    // protect against nearly parallel angles
    static const double maxfac(10.0); // beyond this the straw irregularities dominate and the estimate is unreliable
    static const double minsin2(1.0/(maxfac*maxfac)); // minimum sin^2(theta) this implies
    double sin2 = std::max(1.0 -dirdot*dirdot,minsin2);
    return 1.0/sqrt(sin2);
  }

  void StrawMaterial::findXings(ClosestApproachData const& cadata,StrawXingConfig const& caconfig, std::vector<MaterialXing>& mxings) const {
    mxings.clear();
    double wallpath, gaspath, wirepath;
    pathLengths(cadata,caconfig,wallpath, gaspath, wirepath);
    if(wallpath > 0.0) mxings.push_back(MaterialXing(*wallmat_,wallpath));
    if(gaspath > 0.0) mxings.push_back(MaterialXing(*gasmat_,gaspath));
    if(wirepath > 0.0) mxings.push_back(MaterialXing(*wiremat_,wirepath));
  }

}
