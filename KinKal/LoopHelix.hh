#ifndef KinKal_LoopHelix_hh
#define KinKal_LoopHelix_hh
//
// class desribing the looping helix basis for the kinematic Kalman fit
// It provides geometric, kinematic, and algebraic representation of
// a particle executing a multi-loop helix in a constant magnetic field.
// Original Author David Brown (LBNL) 1/2020
//

#include "KinKal/Vectors.hh"
#include "KinKal/TimeRange.hh"
#include "KinKal/Parameters.hh"
#include "KinKal/MomBasis.hh"
#include "KinKal/ParticleState.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "Math/Rotation3D.h"
#include <vector>
#include <string>
#include <ostream>

namespace KinKal {

  class LoopHelix {
    public:
      // This class must provide the following to be used to instantiate the 
      // classes implementing the Kalman fit
      // define the indices and names of the parameters
      enum ParamIndex {rad_=0,lam_=1,cx_=2,cy_=3,phi0_=4,t0_=5,npars_=6};
      constexpr static ParamIndex t0Index() { return t0_; }
      
      static std::vector<std::string> const& paramNames();
      static std::vector<std::string> const& paramUnits();
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(ParamIndex index);
      static std::string const& paramUnit(ParamIndex index);
      static std::string const& paramTitle(ParamIndex index);
      static std::string const& trajName();

      // interface needed for KKTrk instantiation
      // construct from momentum, position, and particle properties.
      // This also requires the nominal BField, which can be a vector (3d) or a scalar (B along z)
      LoopHelix(VEC4 const& pos, MOM4 const& mom, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange());
      LoopHelix(VEC4 const& pos, MOM4 const& mom, int charge, double bnom, TimeRange const& range=TimeRange()); // do I really need this?
      // construct from the particle state at a given time, plus mass and charge
      LoopHelix(ParticleState const& pstate, double time, double mass, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange());
      // same, including covariance information
      LoopHelix(ParticleStateMeasurement const& pstate, double time, double mass, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange());
      // copy payload and adjust parameters to correspond to a different BField at a particular time
      LoopHelix(LoopHelix const& other, VEC3 const& bnom, double trot);
      // copy payload and override the parameters; Is this used?
      LoopHelix(Parameters const& pdata, LoopHelix const& other);
      VEC4 pos4(double time) const;
      void position(VEC4& pos) const; // time of pos is input 
      VEC3 position(double time) const;
      VEC3 velocity(double time) const;
      double speed(double time) const  {  return CLHEP::c_light*beta(); }
      void print(std::ostream& ost, int detail) const;
      TimeRange const& range() const { return trange_; }
      TimeRange& range() { return trange_; }
      void setRange(TimeRange const& trange) { trange_ = trange; }
      // allow resetting the BField.  Note this is time-dependent
      void setBNom(double time, VEC3 const& bnom);
      bool inRange(double time) const { return trange_.inRange(time); }
      MOM4 momentum(double time) const;
      double momentumMag(double time) const  { return  fabs(mass_*betaGamma()); }
      double momentumVar(double time) const;
      double energy(double time) const  { return  fabs(mass_*ebar()/mbar_); }
      VEC3 direction(double time, MomBasis::Direction mdir= MomBasis::momdir_) const;
      double mass() const { return mass_;} // mass 
      int charge() const { return charge_;} // charge in proton charge units
      double paramVal(size_t index) const { return pars_.parameters()[index]; }
      Parameters const& params() const { return pars_; }
      Parameters& params() { return pars_; }
      // named parameter accessors
      double rad() const { return paramVal(rad_); }
      double lam() const { return paramVal(lam_); }
      double cx() const { return paramVal(cx_); }
      double cy() const { return paramVal(cy_); }
      double phi0() const { return paramVal(phi0_); }
      double t0() const { return paramVal(t0_); }
      // express fit results as a state vector (global coordinates)
      ParticleState state(double time) const;
      ParticleStateMeasurement measurementState(double time) const;
      // simple functions; these can be cached if they cause performance problems
      double sign() const { return copysign(1.0,mbar_); } // combined bending sign including Bz and charge
      double pbar2() const { return  rad()*rad() + lam()*lam(); } 
      double pbar() const { return  sqrt(pbar2()); } // momentum in mm
      double ebar2() const { return  pbar2() + mbar_*mbar_; }
      double ebar() const { return  sqrt(ebar2()); } // energy in mm
      double mbar() const { return mbar_; } // mass in mm; includes charge information!
      double Q() const { return mass_/mbar_; } // reduced charge
      double omega() const { return CLHEP::c_light*sign()/ ebar(); } // rotational velocity, sign set by magnetic force
      double beta() const { return pbar()/ebar(); } // relativistic beta
      double gamma() const { return fabs(ebar()/mbar_); } // relativistic gamma
      double betaGamma() const { return fabs(pbar()/mbar_); } // relativistic betagamma
      double dphi(double t) const { return omega()*(t - t0()); }
      double phi(double t) const { return dphi(t) + phi0(); }
      double ztime(double zpos) const { return t0() + zpos/(omega()*lam()); }
      double zphi(double zpos) const { return zpos/lam() + phi0(); }
      VEC3 const& bnom(double time=0.0) const { return bnom_; }
      double bnomR() const { return bnom_.R(); }
      // flip the helix in time and charge; it remains unchanged geometrically
      void invertCT() {
	mbar_ *= -1.0;
	charge_ *= -1;
	pars_.parameters()[t0_] *= -1.0;
      }
      // functions related to euclidean space to parameter space derivatives
      DPDV dPardX(double time) const; // return the derivative of the parameters WRT the (global) position vector
      DPDV dPardM(double time) const; // return the derivative of the parameters WRT the (global) momentum vector
      DVDP dXdPar(double time) const; // return the derivative of the (global) position vector WRT the parameters
      DVDP dMdPar(double time) const; // return the derivative of the (global) momentum vector WRT parameters
      DSDP dPardState(double time) const; // derivative of parameters WRT global state
      DPDS dStatedPar(double time) const; // derivative of global state WRT parameters
      DVEC momDeriv(double time, MomBasis::Direction mdir) const; // projection of M derivatives onto direction basis
      // package the above for full (global) state
      // Parameter derivatives given a change in BField
      DVEC dPardB(double time) const; // parameter derivative WRT change in BField magnitude
      DVEC dPardB(double time, VEC3 const& BPrime) const; // parameter change given a new BField vector
    private :
// local coordinate system functions, used internally
      VEC3 localDirection(double time, MomBasis::Direction mdir= MomBasis::momdir_) const;
      VEC3 localMomentum(double time) const;
      VEC3 localPosition(double time) const;
      DPDV dPardXLoc(double time) const; // return the derivative of the parameters WRT the local (unrotated) position vector
      DPDV dPardMLoc(double time) const; // return the derivative of the parameters WRT the local (unrotated) momentum vector
      DSDP dPardStateLoc(double time) const; // derivative of parameters WRT local state

      TimeRange trange_;
      Parameters pars_; // parameters
      double mass_;  // in units of MeV/c^2
      int charge_; // charge in units of proton charge
      double mbar_;  // reduced mass in units of mm, computed from the mass and nominal field
      VEC3 bnom_; // nominal BField, in global coordinate system
      ROOT::Math::Rotation3D l2g_, g2l_; // rotations between local and global coordinates 
      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
      static std::vector<std::string> paramUnits_;
      static std::string trajName_;
      // non-const accessors
      double& param(size_t index) { return pars_.parameters()[index]; }
 };
  std::ostream& operator <<(std::ostream& ost, LoopHelix const& lhel);
}
#endif
