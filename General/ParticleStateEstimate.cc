#include "KinKal/General/ParticleStateEstimate.hh"
namespace KinKal {
  double ParticleStateEstimate::momentumVariance() const {
    auto momdir = momentum3().Unit();
    DVEC dMdm(0.0, 0.0, 0.0, momdir.X(), momdir.Y(), momdir.Z());
    return ROOT::Math::Similarity(dMdm,scovar_);
  }
}
