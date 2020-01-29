#ifndef KinKal_KTraj_hh
#define KinKal_KTraj_hh
//
//  Base class for a trajectory used to represent a particle's path through space 
//  Used as part of the kinematic Kalman fit
//  Templated on the dimension of the parameter vector
//  The geometric and kinematic interpretation of the parameters is defined in the subclasses
//
#include "BTrk/KinKal/TTraj.hh"
#include "BTrk/KinKal/Constants.hh"
namespace KinKal {
  template <size_t PDIM> class KTraj : public TTraj<PDIM> {
    public:
      typedef typename TTraj<PDIM>::PVec PVec;
      typedef typename TTraj<PDIM>::PMat PMat;
      typedef typename TTraj<PDIM>::PDer PDer;
      // define local basis vector indices; along and perp to the local momentum.  theta2 is also perpendicular to z
      enum trajdir {momdir=0,theta1,theta2};
      // unit vectors in the different local directions
      virtual void dirVector(trajdir dir,double time,Vec3& unit) const = 0;
      
   // direct accessors
      double mass() const { return mass_;} // mass in MeV/c^2
      int charge() const { return charge_;} // charge in proton charge units

    // kinematic accessors
      virtual void momentum(double t,Mom4& mom) const =0; // momentum in MeV/c, mass in MeV/c^2 as a function of time
      void momentum(Vec4 const& pos, Mom4& mom) const { return momentum(pos.T(),mom); }

    // Parameter derivatives for a change in momentum along local basis directions
      virtual void momDeriv(trajdir dir, double time, PDer& der) const = 0;

    protected:
      KTraj(double mass, int charge) : mass_(mass), charge_(charge) {}
      // construct from particle properties  ters 
      KTraj(PVec const& pars, PMat const& pcov,double mass,int charge) :
	TTraj<PDIM>(pars, pcov), mass_(mass), charge_(charge) {}
      // kinematic parameters
      double mass_;  // in units of MeV
      int charge_; // charge in units of proton charge
    private:
      KTraj() : mass_(-1.0), charge_(0) {} //  invalid
  };
}
#endif
