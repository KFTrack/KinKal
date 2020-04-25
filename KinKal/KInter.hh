#ifndef KinKal_KInter_hh
#define KinKal_KInter_hh
//
//  Define the interface for kinematic information from a particle as a function of time.
//  This class is provided as a convenience: KKTrk can be instantiated on any particle trajectory class satisfying
//  the TTraj and KInter interface below, in addition to the functions described in TTraj.hh and the following :
//
//      void momentum(double t,Mom4& mom) const; // momentum in MeV/c, mass in MeV/c^2 as a function of time
//      void momentum(Vec4 const& pos, Mom4& mom) const { return momentum(pos.T(),mom); }
//      double momentum(float time) const; // momentum and energy magnitude in MeV/
//      double momentumVar(float time) const; // variance on momentum value
//      double energy(float time) const; 
//      void rangeInTolerance(TRange& range, BField const& bfield, double tol);
//
//  Used as part of the kinematic Kalman fit
//
#include <string>
namespace KinKal {
  class BField;
  class KInter {
    public:
      // define basis vectors WRT the local momentum direction.  theta2 is also perpendicular to z
      enum MDir {momdir=0,theta1,theta2,ndir};
      static std::string directionName(MDir tdir) {
	switch (tdir) {
	  case momdir:
	    return std::string("momdir");
	  case theta1:
	    return std::string("theta1");
	  case theta2:
	    return std::string("theta2");
	  default:
	    return std::string("unknown");
	}
      }
   // direct accessors
      double mass() const { return mass_;} // mass 
      int charge() const { return charge_;} // charge in proton charge units
    protected:
      KInter(double mass, int charge) : mass_(mass), charge_(charge) {}
      // kinematic parameters
      double mass_;  // in units of MeV/c^2
      int charge_; // charge in units of proton charge
    private:
      KInter() = delete; // no default construction
  };
}
#endif
