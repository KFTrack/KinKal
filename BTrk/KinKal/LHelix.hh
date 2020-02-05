#ifndef KinKal_LHelix_hh
#define KinKal_LHelix_hh
// class desribing the looping helix basis for the kinematic Kalman fit
// It provides geometric, kinematic, and algebraic representation of
// a particule executing a multi-loop helix in a constant magnetic field.
//
// Original Author David Brown (LBNL) 1/2020
//

#include "BTrk/KinKal/Types.hh"
#include "BTrk/KinKal/TTraj.hh"
#include "BTrk/KinKal/KTraj.hh"
#include "BTrk/KinKal/TData.hh"
#include "BTrk/KinKal/Context.hh"
#include "BTrk/KinKal/Constants.hh"
#include <vector>
#include <string>

namespace KinKal {

  class LHelix : public TTraj, public KTraj {
    public:
      // This class must provide the following to be used to instantiate the 
      // classes implementing the Kalman fit
      // define the indices and names of the parameters
      enum paramIndex {rad_=0,lam_=1,cx_=2,cy_=3,phi0_=4,t0_=5,npars_=6};
      constexpr static size_t NParams() { return npars_; }
      typedef TData<npars_> TDATA; // Data payload for this class
      static std::vector<std::string> const& paramNames(); 
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(paramIndex index);
      static std::string const& paramTitle(paramIndex index);

      // construct from momentum, position, and particle properties.
      // This also requires the BField, through Context
      LHelix(Vec4 const& pos, Mom4 const& mom, int charge, Context const& context);

      // particle position and momentum as a function of time
      void position(Vec4& pos) const override; // time is input 
      void position(double t,Vec3& pos) const override; // time is input 
      void momentum(double t,Mom4& mom) const override;
      void velocity(double time, Vec3& vel) const override;
      void direction(double tval,Vec3& dir) const override;

      // local momentum direction basis
      void dirVector(trajdir dir,double time,Vec3& unit) const override;

      // momentum change derivatives; this is required to instantiate a KalTrk using this KTraj
      typedef ROOT::Math::SMatrix<double,npars_,1> PDer; // derivative of parameters
      void momDeriv(trajdir dir, double time, PDer& der) const;

     // named parameter accessors
      double param(size_t index) const { return pars_.vec()[index]; }
      TDATA const& params() const { return pars_; }
      double rad() const { return param(rad_); }
      double lam() const { return param(lam_); }
      double cx() const { return param(cx_); }
      double cy() const { return param(cy_); }
      double phi0() const { return param(phi0_); }
      double t0() const { return param(t0_); }
      
      // simple functions 
      double pbar() const { return  sqrt(rad()*rad() + lam()*lam() ); } // momentum in mm
      double ebar() const { return  sqrt(rad()*rad() + lam()*lam() + mbar_*mbar_); } // energy in mm
      double mbar() const { return mbar_; } // mass in mm
      double Q() const { return mbar_/mass_; } // reduced charge,
      double omega() const { return copysign(c_,mbar_)/ebar(); } // rotational velocity, sign set by magnetic force 
      double beta() const { return pbar()/ebar(); }
      double dphi(double t) const { return omega()*(t - t0()); }
      double phi(double t) const { return dphi(t) + phi0(); }
      double ztime(double zpos) const { return t0() + zpos/(omega()*lam()); }
      double zphi(double zpos) const { return zpos/lam() + phi0(); }

      // flip the helix in time and charge; it remains unchanged geometrically
      void invertCT() {
	mbar_ *= -1.0;
	charge_ *= -1;
	pars_.vec()[t0_] *= -1.0;
      }
      //
    private :
      TDATA pars_; // parameters
      double mbar_;  // reduced mass in units of mm, computed from the mass and nominal field
      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
  };
}
#endif
