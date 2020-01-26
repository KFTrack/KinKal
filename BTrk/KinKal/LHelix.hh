#ifndef KinKal_LHelix_hh
#define KinKal_LHelix_hh
// class desribing the looping helix basis for the kinematic Kalman fit
// It provides geometric, kinematic, and algebraic representation of
// a particule executing a multi-loop helix in a constant magnetic field.
//
// Original Author David Brown (LBNL) 1/2020
//

#include "KinKal/Types.hh"
#include "KinKal/KTraj.hh"
#include "KinKal/Context.hh"
namespace KinKal {

  class LHelix : pulbic KTraj<6> {
    public:
      // This struct must provide the following to be used to instantiate the 
      // classes implementing the Kalman fit
      // define the indices and names of the parameters
      enum paramindex {rad_=0,lamb_=1,cx_=2,cy_=3,phi0_=4,t0_=5,_npars=6};
      static size_t nParams() { return _npars; }
      static const& std::vector<std::string> paramNames { return paramNames_; }
      static const& std::vector<std::string> paramTitles { return paramTitles_; }
      static const& std::string paramName(paramIndex index) { return paramNames_[static_cast<size_t>(index)];}
      static const& std::string paramTitle(paramIndex index) { return paramTitle_[static_cast<size_t>(index)];}

      // construct from momentum, position, and particle properties
      LHelix(FourV pos, FourV mom, double charge, Context& context);

      // particle position and momentum as a function of time
      override void position(FourV& pos) const; // time is input 
      override void mometum(double t,FourV& mom) const;

      // accessors
      double pbar() const { return  sqrt(pars_[rad_]*pars_[rad_] + pars_[lamb_]*pars_[lamb_] ); } // momentum in mm
      double ebar() const { return  sqrt(pars_[rad_]*pars_[rad_] + pars_[lamb_]*pars_[lamb_] + mbar_*mbar_); } // energy in mm
      // angular rotation frequency
      double omega() const { return copysign(c,mbar_)/ebar(); } // rotational velocity, sign set by magnetic force 
      double beta() const { return pbar()/ebar(); }
      double phi(double t) const { return omega()*(t - pars_[t0_]) + pars_[phi0_]; }
      double time(double zpos) const { return pars_[t0_] + zpos/(omega()*pars_[lamb_]); }

      // flip the helix in time; this also reverses the charge
      void reverse() {
	mbar_ *= -1.0;
	charge_ *= -1;
	pars_[t0_] *= -1.0;
      }
      //
    private :
      double mbar_;  // reduced mass in units of mm (computed from the mass);
  };
}
#endif
