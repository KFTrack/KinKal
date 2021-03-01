#ifndef KinKal_CentralHelix_hh
#define KinKal_CentralHelix_hh
//
// class desribing the looping helix basis for the kinematic Kalman fit
// It provides geometric, kinematic, and algebraic representation of
// a particule executing a multi-loop helix in a constant magnetic field.
// Original Author Roberto Soleti (LBNL) 1/2020
//

#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/General/MomBasis.hh"
#include "KinKal/General/PhysicalConstants.h"
#include "Math/Rotation3D.h"
#include <vector>
#include <string>
#include <ostream>

namespace KinKal {

  class CentralHelix {
    public:
      // This class must provide the following to be used to instantiate the
      // classes implementing the Kalman fit
      // define the indices and names of the parameters
      enum ParamIndex {d0_=0,phi0_=1,omega_=2,z0_=3,tanDip_=4,t0_=5,npars_=6};
      constexpr static ParamIndex t0Index() { return t0_; }

      static std::vector<std::string> const &paramNames();
      static std::vector<std::string> const &paramUnits();
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(ParamIndex index);
      static std::string const& paramUnit(ParamIndex index);
      static std::string const& paramTitle(ParamIndex index);
      static std::string const& trajName();

      // interface needed for KKTrk instantiation
      // construct from momentum, position, and particle properties.
      // This also requires the nominal BFieldMap, which can be a vector (3d) or a scalar (B along z)
      CentralHelix(VEC4 const& pos, MOM4 const& mom, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange());
      CentralHelix(VEC4 const& pos, MOM4 const& mom, int charge, double bnom, TimeRange const& range=TimeRange());
      // copy payload and adjust for a different BFieldMap and range
      CentralHelix(CentralHelix const& other, VEC3 const& bnom, double trot);
      // copy and override parameters
      CentralHelix(Parameters const &pdata, CentralHelix const& other);
      // construct from the particle state at a given time, plus mass and charge
      CentralHelix(ParticleState const& pstate, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange()); // TODO
      // same, including covariance information
      CentralHelix(ParticleStateEstimate const& pstate, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange()); //TODO
      // particle position and momentum as a function of time
      VEC4 position4(double time) const;
      VEC3 position3(double time) const;
      MOM4 momentum4(double time) const;
      VEC3 momentum3(double time) const;
      VEC3 velocity(double time) const;
      VEC3 direction(double time, MomBasis::Direction mdir= MomBasis::momdir_) const;
      // scalar momentum and energy in MeV/c units
      double momentum(double time=0) const  { return fabs(mass_ * pbar() / mbar_); }
      double momentumVar(double time=0) const  { return -1.0; }//TODO
      double energy(double time=0) const  { return fabs(mass_ * ebar() / mbar_); }
      // speed in mm/ns
      double speed(double time=0) const  { return CLHEP::c_light * beta(); }
      // local momentum direction basis
      void print(std::ostream& ost, int detail) const;
      TimeRange const& range() const { return trange_; }
      TimeRange& range() { return trange_; }
      void setRange(TimeRange const& trange) { trange_ = trange; }
      void setBNom(double time, VEC3 const& bnom) {} // TODO
      bool inRange(double time) const { return trange_.inRange(time); }

      // momentum change derivatives; this is required to instantiate a KalTrk using this KTraj
      DVEC momDeriv(double time, MomBasis::Direction mdir) const;
      double mass() const { return mass_;} // mass
      int charge() const { return charge_;} // charge in proton charge units

      // named parameter accessors
      double paramVal(size_t index) const { return pars_.parameters()[index]; }
      Parameters const &params() const { return pars_; }
      Parameters &params() { return pars_; }
      double d0() const { return paramVal(d0_); }
      double phi0() const { return paramVal(phi0_); }
      double omega() const { return paramVal(omega_); } // rotational velocity, sign set by magnetic force
      double z0() const { return paramVal(z0_); }
      double tanDip() const { return paramVal(tanDip_); }
      double t0() const { return paramVal(t0_); }
      // express fit results as a state vector (global coordinates)
      ParticleState state(double time) const {  return ParticleState(position4(time),momentum4(time)); }
      ParticleStateEstimate stateEstimate(double time) const { return ParticleStateEstimate(); } // TODO

      // simple functions
      double sign() const { return copysign(1.0,mbar_); } // combined bending sign including Bz and charge
      double pbar() const { return 1./ (omega() * cosDip() ); } // momentum in mm
      double ebar() const { return sqrt(pbar()*pbar() + mbar_ * mbar_); } // energy in mm
      double cosDip() const { return 1./sqrt(1.+ tanDip() * tanDip() ); }
      double sinDip() const { return tanDip()*cosDip(); }
      double mbar() const { return mbar_; } // mass in mm; includes charge information!
      double vt() const { return vt_; }
      double Q() const { return mass_/mbar_; } // reduced charge
      double beta() const { return fabs(pbar()/ebar()); } // relativistic beta
      double gamma() const { return fabs(ebar()/mbar_); } // relativistic gamma
      double betaGamma() const { return fabs(pbar()/mbar_); } // relativistic betagamma
      double Omega() const { return omega()*vt(); } // true angular velocity
      double dphi(double t) const { return Omega()*(t - t0()); } // rotation WRT 0 at a given time
      double phi(double t) const { return dphi(t) + phi0(); } // absolute azimuth at a given time
      double ztime(double zpos) const { return t0() + zpos / (vt()*tanDip()); } // time the particle reaches given z value
      VEC3 const &bnom(double time=0.0) const { return bnom_; }
      double bnomR() const { return bnom_.R(); }
      DPDV dPardX(double time) const; // TODO
      DPDV dPardM(double time) const; // TODO
      DVDP dXdPar(double time) const; // TODO
      DVDP dMdPar(double time) const; // TODO
      PSMAT dPardState(double time) const { return PSMAT(); } // TODO
      PSMAT dStatedPar(double time) const { return PSMAT(); } // TODO
      // package the above for full (global) state
      // Parameter derivatives given a change in BFieldMap
      DVEC dPardB(double time) const { return DVEC(); } // TODO
      DVEC dPardB(double time, VEC3 const& BPrime) const { return DVEC(); } //TODO

      // flip the helix in time and charge; it remains unchanged geometrically
      void invertCT()
      {
        mbar_ *= -1.0;
        charge_ *= -1;
        pars_.parameters()[t0_] *= -1.0;
      }
      //
    private :
      DPDV dPardMLoc(double time) const; // return the derivative of the parameters WRT the local (unrotated) momentum vector
      DPDV dPardXLoc(double time) const;
      TimeRange trange_;
      Parameters pars_; // parameters
      double mass_;  // in units of MeV/c^2
      int charge_; // charge in units of proton charge
      double mbar_;  // reduced mass in units of mm, computed from the mass and nominal field
      VEC3 bnom_;    // nominal BField vector, from the map
      ROOT::Math::Rotation3D l2g_, g2l_; // rotations between local and global coordinates
      const static std::vector<std::string> paramTitles_;
      const static std::vector<std::string> paramNames_;
      const static std::vector<std::string> paramUnits_;
      const static std::string trajName_;
      double vt_; // transverse velocity
      // non-const accessors
      double &param(size_t index) { return pars_.parameters()[index]; }
  };
  std::ostream& operator <<(std::ostream& ost, CentralHelix const& hhel);
}
#endif
