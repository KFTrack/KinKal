#ifndef KinKal_DXing_hh
#define KinKal_DXing_hh
//
//  Describe the material effects of a kinematic trajectory crossing a piece of physical detector (material)
//  Used in the kinematic Kalman fit
//
#include "KinKal/KInter.hh"
#include "KinKal/MatXing.hh"
#include "KinKal/TPocaBase.hh"
#include "KinKal/TDir.hh"
#include <vector>
#include <stdexcept>
#include <array>
#include <limits>
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class DXing {
    public:
      typedef PKTraj<KTRAJ> PKTRAJ;
      // construct from a time
      DXing(float time=-std::numeric_limits<float>::max()) : xtime_(time) {}
      virtual ~DXing() {}
      virtual void update(PKTRAJ const& pktraj) =0;
      virtual void update(PKTRAJ const& pktraj, float xtime) =0; // update including an estimate of the xing time
      virtual void print(std::ostream& ost=std::cout,int detail=0) const =0;
      // accessors
      float crossingTime() const { return xtime_; }
      float& crossingTime() { return xtime_; }
      std::vector<MatXing>const&  matXings() const { return mxings_; }
      // calculate the cumulative material effect from these crossings
      void momEffects(PKTRAJ const& pktraj, TDir tdir, std::array<float,3>& dmom, std::array<float,3>& momvar) const;
    protected:
      float xtime_; // time on the reference trajectory when the xing occured
      std::vector<MatXing> mxings_; // material crossings for this piece of matter
  };

  template <class KTRAJ> void DXing<KTRAJ>::momEffects(PKTRAJ const& pktraj, TDir tdir, std::array<float,3>& dmom, std::array<float,3>& momvar) const {
    // compute the derivative of momentum to energy
    double mom = pktraj.momentum(xtime_);
    double mass = pktraj.mass();
    double dmFdE = sqrt(mom*mom+mass*mass)/(mom*mom); // dimension of 1/E
    if(tdir == TDir::backwards)dmFdE *= -1.0;
    // loop over crossings for this detector piece
    for(auto const& mxing : mxings_){
      // compute FRACTIONAL momentum change and variance on that in the given direction
      momvar[KInter::momdir] += mxing.dmat_.energyLossVar(mom,mxing.plen_,mass)*dmFdE*dmFdE;
      dmom [KInter::momdir]+= mxing.dmat_.energyLoss(mom,mxing.plen_,mass)*dmFdE;
      double scatvar = mxing.dmat_.scatterAngleVar(mom,mxing.plen_,mass);
      // scattering is the same in each direction and has no net effect, it only adds noise
      momvar[KInter::theta1] += scatvar;
      momvar[KInter::theta2] += scatvar;
    }
    // correct for time direction
  }
}
#endif
