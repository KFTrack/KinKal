#ifndef KinKal_ElementXing_hh
#define KinKal_ElementXing_hh
//
//  Describe the material effects of a particle crossing a detector element (piece of the detector)
//  Used in the kinematic Kalman fit
//
#include "KinKal/General/MomBasis.hh"
#include "KinKal/Detector/MaterialXing.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/Detector/Hit.hh"
#include "KinKal/Fit/MetaIterConfig.hh"
#include <vector>
#include <stdexcept>
#include <array>
#include <limits>
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class ElementXing {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      ElementXing() {}
      virtual ~ElementXing() {}
      virtual void updateReference(KTRAJPTR const& ktrajptr) = 0;
      virtual void updateState(MetaIterConfig const& config,bool first) =0;
      virtual double time() const=0; // time the particle crosses thie element
      virtual double transitTime() const=0; // time to cross this element
      virtual KTRAJ const& referenceTrajectory() const =0; // trajectory WRT which the xing is defined
      virtual void print(std::ostream& ost=std::cout,int detail=0) const =0;
      // crossings  without material are inactive
      bool active() const { return mxings_.size() > 0; }
      std::vector<MaterialXing>const&  matXings() const { return mxings_; }
      std::vector<MaterialXing>&  matXings() { return mxings_; }
      // calculate the cumulative material effect from these crossings
      void materialEffects(TimeDir tdir, std::array<double,3>& dmom, std::array<double,3>& momvar) const;
      // sum radiation fraction
      double radiationFraction() const;
   private:
      std::vector<MaterialXing> mxings_; // Effect of each physical material component of this detector element on this trajectory
  };

  template <class KTRAJ> void ElementXing<KTRAJ>::materialEffects(TimeDir tdir, std::array<double,3>& dmom, std::array<double,3>& momvar) const {
    // compute the derivative of momentum to energy, at the reference trajectory
    double mom = referenceTrajectory().momentum(time());
    double mass = referenceTrajectory().mass();
    double dmFdE = sqrt(mom*mom+mass*mass)/(mom*mom); // dimension of 1/E
    if(tdir == TimeDir::backwards)dmFdE *= -1.0;
    // loop over crossings for this detector piece
    for(auto const& mxing : mxings_){
      // compute FRACTIONAL momentum change and variance on that in the given direction
      momvar[MomBasis::momdir_] += mxing.dmat_.energyLossVar(mom,mxing.plen_,mass)*dmFdE*dmFdE;
      dmom [MomBasis::momdir_]+= mxing.dmat_.energyLoss(mom,mxing.plen_,mass)*dmFdE;
      double scatvar = mxing.dmat_.scatterAngleVar(mom,mxing.plen_,mass);
      // scattering is the same in each direction and has no net effect, it only adds noise
      momvar[MomBasis::perpdir_] += scatvar;
      momvar[MomBasis::phidir_] += scatvar;
    }
  }

  template <class KTRAJ> double ElementXing<KTRAJ>::radiationFraction() const {
    double retval(0.0);
    for(auto const& mxing : mxings_)
      retval += mxing.dmat_.radiationFraction(mxing.plen_/10.0); // Ugly conversion to cm FIXME!!
    return retval;
  }

}
#endif
