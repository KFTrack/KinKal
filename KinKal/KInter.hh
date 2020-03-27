#ifndef KinKal_KInter_hh
#define KinKal_KInter_hh
//
//  Define the interface for kinematic information from a particle as a function of time
//  Used as part of the kinematic Kalman fit
//
#include <string>
namespace KinKal {
  class BField;
  class KInter {
    public:
      // define basis vectors WRT the local momentum direction.  theta2 is also perpendicular to z
      enum MDir {momdir=0,theta1,theta2};
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
      // unit vectors in the different local directions
      virtual void dirVector(MDir mdir,double time,Vec3& unit) const = 0;
      
   // direct accessors
      double mass() const { return mass_;} // mass 
      int charge() const { return charge_;} // charge in proton charge units

    // kinematic accessors
      virtual void momentum(double t,Mom4& mom) const =0; // momentum in MeV/c, mass in MeV/c^2 as a function of time
      void momentum(Vec4 const& pos, Mom4& mom) const { return momentum(pos.T(),mom); }
      virtual double momentum(double time) const =0; // momentum and energy magnitude in MeV/
      virtual double energy(double time) const =0; 
      // reduce the end of the given range so that the trajectory position stays within the given spatial tolerance (mm), given the BField
      virtual void rangeInTolerance(TRange& range, BField const& bfield, double tol) const = 0;
      virtual ~KInter(){}
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
