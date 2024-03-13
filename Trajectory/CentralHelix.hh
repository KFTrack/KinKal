#ifndef KinKal_CentralHelix_hh
#define KinKal_CentralHelix_hh
//
// class desribing a helix in terms of its impact parameter to the z axis, direction, and curvature.
// It provides geometric, kinematic, and algebraic representation
// a particule moving along an arc of a helix in a constant magnetic field.
// Original Author Roberto Soleti (LBNL) 1/2020
//

#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/ParticleStateEstimate.hh"
#include "KinKal/General/MomBasis.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/PhysicalConstants.h"
#include "KinKal/Geometry/Ray.hh"
#include "KinKal/Geometry/Plane.hh"
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
      constexpr static ParamIndex phi0Index() { return phi0_; }
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
      // construct from explicit parametric and kinematic info
      CentralHelix(Parameters const &pdata, double mass, int abscharge, double bnom, TimeRange const& range);
      // copy payload and adjust for a different BFieldMap and range
      CentralHelix(CentralHelix const& other, VEC3 const& bnom, double trot);
      // copy and override parameters
      CentralHelix(Parameters const &pdata, CentralHelix const& other);
      // construct from the particle state.  Requires the BField
      explicit CentralHelix(ParticleState const& pstate, VEC3 const& bnom, TimeRange const& range=TimeRange());
      // same, including covariance information
      explicit CentralHelix(ParticleStateEstimate const& pstate, VEC3 const& bnom, TimeRange const& range=TimeRange());
      // particle position and momentum as a function of time
      VEC4 position4(double time) const;
      VEC3 position3(double time) const;
      MOM4 momentum4(double time) const;
      VEC3 momentum3(double time) const;
      VEC3 velocity(double time) const;
      double speed(double time=0) const  { return CLHEP::c_light * beta(); }
      double acceleration() const { return CLHEP::c_light*CLHEP::c_light/(omega()*ebar2()); }
      VEC3 acceleration(double time) const;
      VEC3 direction(double time, MomBasis::Direction mdir= MomBasis::momdir_) const;
      // scalar momentum and energy in MeV/c units
      double momentum(double time=0) const  { return fabs(mass_ * pbar() / mbar()); }
      double momentumVariance(double time=0) const;
      double positionVariance(double time,MomBasis::Direction dir) const;
      PMAT planeCovariance(double time,Plane const& plane) const;
      double energy(double time=0) const  { return fabs(mass_ * ebar() / mbar()); }
      // local momentum direction basis
      void print(std::ostream& ost, int detail) const;
      TimeRange const& range() const { return trange_; }
      TimeRange& range() { return trange_; }
      void setRange(TimeRange const& trange) { trange_ = trange; }
      void setBNom(double time, VEC3 const& bnom);
      // change the BField.  This also resets the transforms
      void resetBNom(VEC3 const& bnom);
      bool inRange(double time) const { return trange_.inRange(time); }

      // momentum change derivatives; this is required to instantiate a KalTrk using this KTraj
      DVEC momDeriv(double time, MomBasis::Direction mdir) const;
      double mass() const { return mass_;} // mass
      int charge() const { return (-1*omega() < 0) ? -1*abscharge_ : abscharge_;}
      // named parameter accessors
      double paramVal(size_t index) const { return pars_.parameters()[index]; }
      double paramVar(size_t index) const { return pars_.covariance()(index,index); }
      Parameters const &params() const { return pars_; }
      Parameters &params() { return pars_; }
      double d0() const { return paramVal(d0_); }
      double phi0() const { return paramVal(phi0_); }
      double omega() const { return paramVal(omega_); } // rotational velocity, sign set by magnetic force
      double z0() const { return paramVal(z0_); }
      double tanDip() const { return paramVal(tanDip_); }
      double t0() const { return paramVal(t0_); }
      // express fit results as a state vector (global coordinates)
      ParticleState state(double time) const {  return ParticleState(position4(time),momentum4(time),charge()); }
      ParticleStateEstimate stateEstimate(double time) const;

      // simple functions
      double sign() const { return copysign(1.0,omega()); } // combined bending sign including Bz and charge
      double parameterSign() const { return copysign(1.0,omega()); }
      // helicity is defined as the sign of the projection of the angular momentum vector onto the linear momentum vector
      double helicity() const { return copysign(1.0,tanDip()); } // needs to be checked TODO
      double pbar() const { return 1./ (omega() * cosDip() ); } // momentum in mm
      double ebar2() const { return pbar()*pbar() + mbar() * mbar(); }
      double ebar() const { return sqrt(ebar2()); } // energy in mm
      double mbar() const { return fabs(mass_/Q()); } // mass in mm
      double Q() const { return -BFieldMap::cbar()*charge()*bnom_.R(); } // reduced charge
      double cosDip() const { return 1./sqrt(1.+ tanDip() * tanDip() ); }
      double sinDip() const { return tanDip()*cosDip(); }
      double omegaZ() const { return omega()/(CLHEP::c_light*beta()*tanDip()); } // dPhi/dz
      double beta() const { return fabs(pbar()/ebar()); } // relativistic beta
      double gamma() const { return fabs(ebar()/mbar()); } // relativistic gamma
      double betaGamma() const { return fabs(pbar()/mbar()); } // relativistic betagamma
      double Omega() const { return Q()*CLHEP::c_light/energy(); } // true angular velocity
      double dphi(double t) const { return Omega()*(t - t0()); } // rotation WRT 0 at a given time
      double phi(double t) const { return dphi(t) + phi0(); } // absolute azimuth at a given time
      VEC3 const& bnom(double time=0.0) const { return bnom_; }
      double bnomR() const { return bnom_.R(); }
      DPDV dPardX(double time) const;
      DPDV dPardM(double time) const;
      DVDP dXdPar(double time) const;
      DVDP dMdPar(double time) const;
      PSMAT dPardState(double time) const;
      PSMAT dStatedPar(double time) const;
      // package the above for full (global) state
      // Parameter derivatives given a change in BFieldMap
      DVEC dPardB(double time) const;
      DVEC dPardB(double time, VEC3 const& BPrime) const;
      PSMAT dPardPardB(double time,VEC3 const& db) const; // Parameter covariance rotation for a change in BField

      // flip the helix in time and charge; it remains unchanged geometrically
      void invertCT() {
        pars_.parameters()[omega_] *= -1.0;
        pars_.parameters()[tanDip_] *= -1.0;
        pars_.parameters()[d0_] *= -1.0;
        pars_.parameters()[phi0_] += M_PI;
        pars_.parameters()[t0_] *= -1.0;
      }
      // helix interface
      VEC3 center(double time) const; // helix center in global coordinates
      Ray axis(double time) const; // helix axis in global coordinates
      double axisSpeed() const; // speed along the axis direction (always positive)
      double bendRadius() const { return fabs(1.0/omega()); }
    private :
      VEC3 localDirection(double time, MomBasis::Direction mdir= MomBasis::momdir_) const;
      VEC3 localMomentum(double time) const;
      VEC3 localPosition(double time) const;
      DPDV dPardMLoc(double time) const; // return the derivative of the parameters WRT the local (unrotated) momentum vector
      DPDV dPardXLoc(double time) const;
      PSMAT dPardStateLoc(double time) const; // derivative of parameters WRT local state
      double rc() const { return -1.0/omega() - d0(); }
      VEC3 center() const { return VEC3(rc()*sin(phi0()), -rc()*cos(phi0()), 0.0); } // local circle center
      void setTransforms(); // define global to local and local to global given BNom
      TimeRange trange_;
      Parameters pars_; // parameters
      double mass_;  // in units of MeV/c^2
      int abscharge_; // absolute value of charge in units of proton charge; we need to take the sign of the charge from omega
      VEC3 bnom_;    // nominal BField vector, from the map
      ROOT::Math::Rotation3D l2g_, g2l_; // rotations between local and global coordinates
      const static std::vector<std::string> paramTitles_;
      const static std::vector<std::string> paramNames_;
      const static std::vector<std::string> paramUnits_;
      const static std::string trajName_;
      // DO NOT CACHE ANYTHING that depends on parameters.  it will break the parameter-based constructors
      // non-const accessors
      double &param(size_t index) { return pars_.parameters()[index]; }
  };
  std::ostream& operator <<(std::ostream& ost, CentralHelix const& hhel);
}
#endif
