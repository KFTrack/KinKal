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
}
