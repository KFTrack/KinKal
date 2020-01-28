#ifndef KinKal_LHelix_hh
#define KinKal_LHelix_hh
// class desribing the looping helix basis for the kinematic Kalman fit
// It provides geometric, kinematic, and algebraic representation of
// a particule executing a multi-loop helix in a constant magnetic field.
//
// Original Author David Brown (LBNL) 1/2020
//

#include "BTrk/KinKal/Types.hh"
#include "BTrk/KinKal/KTraj.hh"
#include "BTrk/KinKal/Context.hh"
#include <vector>
#include <string>

namespace KinKal {

  class LHelix : public KTraj<6> {
    public:
      // This struct must provide the following to be used to instantiate the 
      // classes implementing the Kalman fit
      // define the indices and names of the parameters
      enum paramIndex {rad_=0,lam_=1,cx_=2,cy_=3,phi0_=4,t0_=5,npars_=6};
      static size_t nParams() { return npars_; }
      static std::vector<std::string> const& paramNames(); 
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(paramIndex index);
      static std::string const& paramTitle(paramIndex index);

      // construct from momentum, position, and particle properties
      LHelix(Vec4 const& pos, Mom4 const& mom, int charge, Context const& context);

      // particle position and momentum as a function of time
      void position(Vec4& pos) const override; // time is input 
      void position(double t,Vec3& pos) const override; // time is input 
      void momentum(double t,Mom4& mom) const override;
      void velocity(double time, Vec3& vel) const override;
      void direction(double tval,Vec3& dir) const override;

      // local basis
      void dirVector(trajdir dir,double time,Vec3& unit) const override;

      // momentum change derivatives
      void momDeriv(trajdir dir, double time, PDer& der) const override;

      // local accessors
      double pbar() const { return  sqrt(pars_[rad_]*pars_[rad_] + pars_[lam_]*pars_[lam_] ); } // momentum in mm
      double ebar() const { return  sqrt(pars_[rad_]*pars_[rad_] + pars_[lam_]*pars_[lam_] + mbar_*mbar_); } // energy in mm
      double mbar() const { return mbar_; } // mass in mm
      double Q() const { return mbar_/mass_; } // reduced charge,
      double omega() const { return copysign(c_,mbar_)/ebar(); } // rotational velocity, sign set by magnetic force 
      double beta() const { return pbar()/ebar(); }
      double dphi(double t) const { return omega()*(t - pars_[t0_]); }
      double phi(double t) const { return dphi(t) + pars_[phi0_]; }
      double time(double zpos) const { return pars_[t0_] + zpos/(omega()*pars_[lam_]); }

      // flip the helix in time; this also reverses the charge
      void reverse() {
	mbar_ *= -1.0;
	charge_ *= -1;
	pars_[t0_] *= -1.0;
      }
      //
    private :
      double mbar_;  // reduced mass in units of mm (computed from the mass);

      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
  };
}
#endif
