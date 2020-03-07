#ifndef KinKal_KTMIsect_hh
#define KinKal_KTMIsect_hh
//
//  Describe the intersection of a kinematic trajectory with a piece of material
//  Used in the kinematic Kalman fit
//
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
    double dmdE = sqrt(mom*mom+mass*mass)/mom;
    for(auto const& misect : misects_){
      switch (mdir) {
	case KInter::momdir:
	  switch (tdir) {
	    case TDir::forwards:
	      dmom += misect.dmat_.energyLoss(mom,misect.plen_,mass)*dmdE;
	      momvar += misect.dmat_.energyLossVar(mom,misect.plen_,mass)*dmdE*dmdE;
	      break;
	    case TDir::backwards:
	      dmom += misect.dmat_.energyGain(mom,misect.plen_,mass)*dmdE;
	      momvar += misect.dmat_.energyLossVar(mom,misect.plen_,mass)*dmdE*dmdE;
	      break;
	    default :
	      throw std::invalid_argument("Invalid direction");
	  }
	  break;
	  // scattering is the same in both directions
	case KInter::theta1 : case KInter::theta2:
	  // time direction doesn't matter for scattering
	  momvar += misect.dmat_.scatterAngleVar(mom,misect.plen_,mass);
	  break;
	default :
	  throw std::invalid_argument("Invalid direction");
      }
    }
  }

}
#endif
