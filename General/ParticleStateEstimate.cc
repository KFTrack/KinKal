#include "KinKal/General/ParticleStateEstimate.hh"
#include <stdexcept>
namespace KinKal {
  double ParticleStateEstimate::momentumVariance() const {
    auto momdir = momentum3().Unit();
    DVEC dMdm(0.0, 0.0, 0.0, momdir.X(), momdir.Y(), momdir.Z());
    return ROOT::Math::Similarity(dMdm,scovar_);
  }

  double ParticleStateEstimate::positionVariance(MomBasis::Direction dir) const {
    if(dir == MomBasis::momdir_)throw std::invalid_argument("Variance along momentum direction is undefined");
    auto momdir = momentum3();
    auto pdir = MomBasis::direction(dir,momdir);
    DVEC dPdp(pdir.X(),pdir.Y(),pdir.Z(),0.0, 0.0, 0.0);
    return ROOT::Math::Similarity(dPdp,scovar_);
  }

  PMAT ParticleStateEstimate::planeCovariance() const {
    auto momdir = momentum3();
    auto udir = MomBasis::direction(MomBasis::perpdir_,momdir);
    auto vdir = MomBasis::direction(MomBasis::phidir_,momdir);
    SVEC6 uvec(udir.X(),udir.Y(),udir.Z(),0.0, 0.0, 0.0);
    SVEC6 vvec(vdir.X(),vdir.Y(),vdir.Z(),0.0, 0.0, 0.0);
    PPMAT ppmat;
    ppmat.Place_in_row(uvec,0,0);
    ppmat.Place_in_row(vvec,1,0);
    return ROOT::Math::Similarity(ppmat,scovar_);
  }
}
