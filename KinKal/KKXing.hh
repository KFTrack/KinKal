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
      void momEffect(TDir tdir, KInter::MDir mdir, double& dmom, double& momvar) const;
    private:
      KTRAJ const& ktraj_; // kinematic trajectory which crosses the material
      double xtime_; // time of the crossing
      std::vector<MatXing> mxings_; // material crossings for this piece of matter
  };


  template <class KTRAJ> void KKXing<KTRAJ>::momEffect(TDir tdir, KInter::MDir mdir, double& dmom, double& momvar) const {
    dmom = momvar = 0.0;
    // compute the derivative of momentum to energy
    double mom = ktraj_.momentum(xtime_);
    double mass = ktraj_.mass();
    double dmFdE = sqrt(mom*mom+mass*mass)/(mom*mom); // dimension of 1/E
    for(auto const& mxing : mxings_){
      switch (mdir) {
	case KInter::momdir:
	  // compute FRACTIONAL momentum change and variance on that in the given direction
	  momvar += mxing.dmat_.energyLossVar(mom,mxing.plen_,mass)*dmFdE*dmFdE;
	  dmom += mxing.dmat_.energyLoss(mom,mxing.plen_,mass)*dmFdE;
	case KInter::theta1 : case KInter::theta2:
	  // scattering is the same in each direction and has no net effect, it only adds noise
	  momvar += mxing.dmat_.scatterAngleVar(mom,mxing.plen_,mass);
	  break;
	default :
	  throw std::invalid_argument("Invalid direction");
      }
    }
    // correct for time direction
    if(tdir ==TDir::backwards && mdir ==KInter::momdir) dmom *= -1.0;
  }
}
#endif
