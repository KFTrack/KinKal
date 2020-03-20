#ifndef KinKal_KTMIsect_hh
#define KinKal_KTMIsect_hh
//
//  Describe the intersection of a kinematic trajectory with a piece of material
//  Used in the kinematic Kalman fit
//
#include "MatEnv/DetMaterial.hh"
#include "KinKal/KInter.hh"
#include "KinKal/MIsect.hh"
#include "KinKal/TDir.hh"
#include <vector>
#include <stdexcept>
namespace KinKal {
// struct to describe a particular interaction with a particular material
  template <class KTRAJ> class KTMIsect {
    public:
      // construct from a trajectory and a time:
      KTMIsect(KTRAJ const& ktraj,double tisect, std::vector<MIsect> misects) : ktraj_(ktraj), tisect_(tisect), misects_(misects) {}
      // append an interaction
      void addMIsect(MIsect const& misect) { misects_.push_back(misect); }
      // accessors
      KTRAJ const& kTraj() const { return ktraj_; }
      double intersectTime() const { return tisect_; }
      std::vector<MIsect>const&  matIntersections() { return misects_; }
      // calculate the cumulative material effect from these intersections
      void momEffect(TDir tdir, KInter::MDir mdir, double& dmom, double& momvar) const;
    private:
      KTRAJ const& ktraj_; // kinematic trajectory which intersects the material
      double tisect_; // time of the intersection
      std::vector<MIsect> misects_; // material intersections for this piece of matter
  };


  template <class KTRAJ> void KTMIsect<KTRAJ>::momEffect(TDir tdir, KInter::MDir mdir, double& dmom, double& momvar) const {
    dmom = momvar = 0.0;
    // compute the derivative of momentum to energy
    double mom = ktraj_.momentum(tisect_);
    double mass = ktraj_.mass();
    double dmFdE = sqrt(mom*mom+mass*mass)/(mom*mom); // dimension of 1/E
    for(auto const& misect : misects_){
      switch (mdir) {
	case KInter::momdir:
	  // compute FRACTIONAL momentum change and variance on that in the given direction
	  momvar += misect.dmat_.energyLossVar(mom,misect.plen_,mass)*dmFdE*dmFdE;
	  dmom += misect.dmat_.energyLoss(mom,misect.plen_,mass)*dmFdE;
	case KInter::theta1 : case KInter::theta2:
	  // scattering is the same in each direction and has no net effect, it only adds noise
	  momvar += misect.dmat_.scatterAngleVar(mom,misect.plen_,mass);
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
