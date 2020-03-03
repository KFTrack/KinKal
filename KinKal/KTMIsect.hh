#ifndef KinKal_KTMIsect_hh
#define KinKal_KTMIsect_hh
//
//  Describe the intersection of a kinematic trajectory with a piece of material
//  Used in the kinematic Kalman fit
//
#include "KinKal/KInter.hh"
#include "KinKal/TDir.hh"
#include "MatEnv/DetMaterial.hh"
#include <vector>
#include <stdexcept>
namespace KinKal {
// struct to describe a particular interaction with a particular material
  struct MInter {
    DetMaterial const& dmat_; // material
    double plen_; // path length through this material
    MInter(DetMaterial const& dmat,double plen) : dmat_(dmat), plen_(plen) {}
  };

  template <class KTRAJ> class KTMIsect {
    public:
      // construct from a trajectory and a time:
      KTMIsect(KTRAJ const& ktraj,double tinter) : ktraj_(ktraj), tinter_(tinter) {}
      // append an interaction
      void addMInter(MInter const& minter) { minters_.push_back(minter); }
      // accessors
      KTRAJ const& kTraj() const { return ktraj_; }

      // calculate the cumulative material effect from these intersections
      void momEffect(TDir tdir, KInter::MDir mdir, double& dmom, double& momvar) {
	dmom = momvar = 0.0;
	// compute the derivative of momentum to energy
	double mom = ktraj_.momentum(tinter_);
	double mass = ktraj_.mass();
	double dmdE = sqrt(mom*mom+mass*mass)/mom;
	for(auto const& minter : minters_){
	  switch (mdir) {
	    case momdir:
	      switch (tdir) {
		case TDir::forwards:
		  dmom += minter.dmat_.energyLoss(mom,minter.plen_,mass)*dmdE;
		  momvar += minter.dmat_.energyLossVar(mom,minter.plen_,mass)*dmdE*dmdE;
		  break;
		case TDir::backwards:
		  dmom += minter.dmat_.energyGain(mom,minter.plen_,mass)*dmdE;
		  momvar += minter.dmat_.energyLossVar(mom,minter.plen_,mass)*dmdE*dmdE;
		  break;
		default :
		  throw std::invalid_argument("Invalid direction");
	      }
	      break;
	      // scattering is the same in both directions
	    case theta1 : case theta2:
	      // time direction doesn't matter for scattering
	      dmom += minter.dmat_.scatterAngle(mom,minter.plen_,mass);
	      momvar += minter.dmat_.scatterAngleVar(mom,minter.plen_,mass);
	      break;
	    default :
	      throw std::invalid_argument("Invalid direction");
	  }
	}
      }
    private:
      KTRAJ const& ktraj_; // kinematic trajectory which intersects the material
      double tinter_; // time of the intersection
      std::vector<MInter> minters_; // material intersections for this piece of matter
  };
}
#endif
