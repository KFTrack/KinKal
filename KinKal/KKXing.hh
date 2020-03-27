#ifndef KinKal_KKXing_hh
#define KinKal_KKXing_hh
//
//  Describe the a kinematic trajectory crossing a piece of physical detector (material)
//  Used in the kinematic Kalman fit
//
#include "KinKal/KInter.hh"
#include "KinKal/MatXing.hh"
#include "KinKal/TDir.hh"
#include <vector>
#include <stdexcept>
#include <array>
namespace KinKal {
// struct to describe a particular trajectory crossing a particular detector piece, potentially mulitple materials
  template <class KTRAJ> class KKXing {
    public:
      // construct from a trajectory and a time:
      KKXing(KTRAJ const& ktraj,double xtime, std::vector<MatXing> mxings) : ktraj_(ktraj), xtime_(xtime), mxings_(mxings) {}
      // append an interaction
      void addMatXing(MatXing const& mxing) { mxings_.push_back(mxing); }
      // accessors
      KTRAJ const& kTraj() const { return ktraj_; }
      double crossingTime() const { return xtime_; }
      std::vector<MatXing>const&  matXings() { return mxings_; }
      // calculate the cumulative material effect from these crossings
      void momEffects(TDir tdir, std::array<double,3>& dmom, std::array<double,3>& momvar) const;
    private:
      KTRAJ const& ktraj_; // kinematic trajectory which crosses the material
      double xtime_; // time of the crossing
      std::vector<MatXing> mxings_; // material crossings for this piece of matter
  };


  template <class KTRAJ> void KKXing<KTRAJ>::momEffects(TDir tdir, std::array<double,3>& dmom, std::array<double,3>& momvar) const {
    // compute the derivative of momentum to energy
    double mom = ktraj_.momentum(xtime_);
    double mass = ktraj_.mass();
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
