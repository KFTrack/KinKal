#ifndef KinKal_DetectorXing_hh
#define KinKal_DetectorXing_hh
//
//  Describe the material effects of a kinematic trajectory crossing a piece of physical detector (material)
//  Used in the kinematic Kalman fit
//
#include "KinKal/MomBasis.hh"
#include "KinKal/MaterialXing.hh"
#include "KinKal/ParticleTrajectory.hh"
#include "KinKal/TimeDir.hh"
#include <vector>
#include <stdexcept>
#include <array>
#include <limits>
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class DetectorXing {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      // construct from a time
      DetectorXing(double time=-std::numeric_limits<double>::max()) : xtime_(time) {}
      virtual ~DetectorXing() {}
      virtual void update(PKTRAJ const& pktraj,double precision) =0;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const =0;
      // accessors
      double crossingTime() const { return xtime_; }
      double& crossingTime() { return xtime_; }
      std::vector<MaterialXing>const&  matXings() const { return mxings_; }
      std::vector<MaterialXing>&  matXings() { return mxings_; }
      // calculate the cumulative material effect from these crossings
      void materialEffects(PKTRAJ const& pktraj, TimeDir tdir, std::array<double,3>& dmom, std::array<double,3>& momvar) const;
    private:
      double xtime_; // time on the reference trajectory when the xing occured
      std::vector<MaterialXing> mxings_; // material crossings for this detector piece on this trajectory
  };

  template <class KTRAJ> void DetectorXing<KTRAJ>::materialEffects(PKTRAJ const& pktraj, TimeDir tdir, std::array<double,3>& dmom, std::array<double,3>& momvar) const {
    // compute the derivative of momentum to energy
    double mom = pktraj.momentumMag(xtime_);
    double mass = pktraj.mass();
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
    // correct for time direction
  }
}
#endif
