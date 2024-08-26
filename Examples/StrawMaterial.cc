#include "KinKal/Examples/StrawMaterial.hh"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace KinKal {
  void StrawMaterial::pathLengths(ClosestApproachData const& cadata,StrawXingConfig const& caconfig,
      double& wallpath, double& gaspath, double& wirepath) const {
    wallpath = gaspath = wirepath = 0.0;
    double adoca = fabs(cadata.doca());
    double sigdoca = sqrt(cadata.docaVar());
    if(adoca < caconfig.maxdoca_){
      if ((!caconfig.average_) && sigdoca < caconfig.minsigdoca_ && adoca < caconfig.maxddoca_) {  // use exact calculation based on DOCA
        double doca = std::min(adoca, grad_); // truncate
        double ddoca = doca*doca;
        gaspath = 2.0*sqrt(grad2_- ddoca);
        wallpath = 2.0*thick_*srad_/sqrt(srad2_-ddoca);
      } else {
        // errors are large WRT the size of the straw, or DOCA is very far from the wire: just take the average over all impact parameters
        gaspath = M_PI_2*srad_;
        wallpath = M_PI*thick_;
      }
      if(isnan(wallpath) || isnan(gaspath))throw std::runtime_error("Invalid StrawMaterial pathlength");
      // Model the wire as a diffuse gas, density constrained by DOCA TODO
      // correct for the angle WRT the axis
      double afac = angleFactor(cadata.dirDot());
      wallpath *= afac;
      gaspath *= afac;
    }
  }

  double StrawMaterial::transitLength(ClosestApproachData const& cadata) const {
    double doca = std::min(fabs(cadata.doca()),grad_);
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
