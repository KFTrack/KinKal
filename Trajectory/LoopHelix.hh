#ifndef KinKal_LoopHelix_hh
#define KinKal_LoopHelix_hh
//
// class desribing the looping helix basis for the kinematic Kalman fit
// It provides geometric, kinematic, and algebraic representation of
// a particle executing a multi-loop helix in a constant magnetic field.
// Original Author David Brown (LBNL) 1/2020
//

#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/MomBasis.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/ParticleStateEstimate.hh"
#include "KinKal/General/PhysicalConstants.h"
#include "KinKal/Geometry/Plane.hh"
#include "KinKal/Geometry/Ray.hh"
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
      constexpr static ParamIndex phi0Index() { return phi0_; }
      constexpr static ParamIndex t0Index() { return t0_; }

      static std::vector<std::string> const& paramNames();
      static std::vector<std::string> const& paramUnits();
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(ParamIndex index);
      static std::string const& paramUnit(ParamIndex index);
      static std::string const& paramTitle(ParamIndex index);
      static std::string const& trajName();
      // default constructor, needed for persistance
      LoopHelix();
      // interface needed for KKTrk instantiation
      // construct from momentum, position, and particle properties.
      // This also requires the nominal BField, which can be a vector (3d) or a scalar (B along z)
      LoopHelix(VEC4 const& pos, MOM4 const& mom, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange());
      LoopHelix(VEC4 const& pos, MOM4 const& mom, int charge, double bnom, TimeRange const& range=TimeRange()); // do I really need this?
      // construct from the particle state at a given time, plus mass and charge. Parameter covariance matrix is undefined
      explicit LoopHelix(ParticleState const& pstate, VEC3 const& bnom, TimeRange const& range=TimeRange());
      // same, including covariance information
      explicit LoopHelix(ParticleStateEstimate const& pstate, VEC3 const& bnom, TimeRange const& range=TimeRange());
      // copy payload and adjust parameters to correspond to a different BField at a particular time
      LoopHelix(LoopHelix const& other, VEC3 const& bnom, double tref);
      // create from parameters and kinematics separately
      LoopHelix( Parameters const& pars, double mass, int charge, VEC3 const& bnom, TimeRange const& trange=TimeRange() );
      // copy payload and override the parameters; Is this used?
      LoopHelix(Parameters const& pdata, LoopHelix const& other);
      // synchronize phi0, which has a 2pi wrapping
      void syncPhi0(LoopHelix const& other);
      VEC4 position4(double time) const;
      VEC3 position3(double time) const;
      VEC3 velocity(double time) const;
      double speed(double time=0.0) const  {  return CLHEP::c_light*beta(); }
      double transverseSpeed() const {  return fabs(CLHEP::c_light*lam()/ebar()); } // speed perpendicular to the axis
      double acceleration() const { return rad()*CLHEP::c_light*CLHEP::c_light/ebar2(); }
      VEC3 acceleration(double time) const;
      void print(std::ostream& ost, int detail) const;
      TimeRange const& range() const { return trange_; }
      TimeRange& range() { return trange_; }
      void setRange(TimeRange const& trange) { trange_ = trange; }
      // allow resetting the BField.  Note this is time-dependent
      void setBNom(double time, VEC3 const& bnom);
      // change the BField.  This also resets the transforms
      void resetBNom(VEC3 const& bnom);
      bool inRange(double time) const { return trange_.inRange(time); }
      VEC3 momentum3(double time) const;
      MOM4 momentum4(double time) const;
      double momentum(double time=0) const  { return  mass_*betaGamma(); }
      double momentumVariance(double time=0) const;
      double positionVariance(double time,MomBasis::Direction dir) const;
      PMAT planeCovariance(double time,Plane const& plane) const;
      double energy(double time=0) const  { return  fabs(ebar()*Q()); }
      VEC3 direction(double time, MomBasis::Direction mdir= MomBasis::momdir_) const;
      double mass() const { return mass_;} // mass
      int charge() const { return charge_;} // charge in proton charge units
      double paramVal(size_t index) const { return pars_.parameters()[index]; }
      double paramVar(size_t index) const { return pars_.covariance()(index,index); }
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
      ParticleState state(double time) const { return ParticleState(position4(time),momentum4(time),charge()); }
      ParticleStateEstimate stateEstimate(double time) const;
      // simple functions
      double sign() const { return copysign(1.0,charge()); } // charge sign
      double parameterSign() const { return copysign(1.0,rad()); }
      // helicity is defined as the sign of the projection of the angular momentum vector onto the linear momentum vector
      double helicity() const { return copysign(1.0,lam()); }
      double pbar2() const { return  rad()*rad() + lam()*lam(); }
      double pbar() const { return  sqrt(pbar2()); } // momentum in mm
      double ebar2() const { double mb = mbar(); return  pbar2() + mb*mb; }
      double ebar() const { return  sqrt(ebar2()); } // energy in mm
      double mbar() const { return fabs(mass_/Q()); } // mass in mm
      double Q() const { return -BFieldMap::cbar()*charge()*bnom_.R(); } // reduced charge
      double omega() const { return -CLHEP::c_light*sign()/ ebar(); } // rotational velocity, sign set by magnetic force
      double omegaZ() const { return 1.0/lam(); } // dPhi/dz
      double beta() const { return pbar()/ebar(); } // relativistic beta
      double gamma() const { return ebar()/mbar(); } // relativistic gamma
      double betaGamma() const { return pbar()/mbar(); } // relativistic betagamma
      double dphi(double t) const { return omega()*(t - t0()); }
      double phi(double t) const { return dphi(t) + phi0(); }
      double zphi(double zpos) const { return zpos/lam() + phi0(); }
      VEC3 const& bnom(double time=0.0) const { return bnom_; }
      double bnomR() const { return bnom_.R(); }
      // flip the helix in time and charge; it remains unchanged geometrically
      void invertCT() {
        charge_ *= -1;
        pars_.parameters()[t0_] *= -1.0;
      }
      // functions related to euclidean space to parameter space derivatives
      DPDV dPardX(double time) const; // return the derivative of the parameters WRT the (global) position vector
      DPDV dPardM(double time) const; // return the derivative of the parameters WRT the (global) momentum vector
      DVDP dXdPar(double time) const; // return the derivative of the (global) position vector WRT the parameters
      DVDP dMdPar(double time) const; // return the derivative of the (global) momentum vector WRT parameters
      PSMAT dPardState(double time) const; // derivative of parameters WRT global state
      PSMAT dStatedPar(double time) const; // derivative of global state WRT parameters
      DVEC momDeriv(double time, MomBasis::Direction mdir) const; // projection of M derivatives onto direction basis
      // package the above for full (global) state
      DVEC dPardB(double time, VEC3 const& db) const; // parameter change given a change in BField vector; this includes the magnitude and direction changes
      PSMAT dPardPardB(double time,VEC3 const& db) const; // Parameter covariance rotation for a change in BField
      // helix interface
      VEC3 center(double time) const; // helix center in global coordinates
      Ray axis(double time) const; // helix axis in global coordinates
      // linear approximation
      Ray linearize(double time) const { return axis(time); }
      // convenience accessors
      VEC2 center() const { return VEC2(cx(),cy()); }
      double minAxisDist() const { return fabs(center().R()-fabs(rad())); } // minimum distance to the axis
      double maxAxisDist() const { return center().R()+fabs(rad()); } // maximum distance to the axis
      double axisSpeed() const; // speed along the axis direction (always positive)
      double bendRadius() const { return fabs(rad());}
      double sagitta(double range) const; // compute maximum sagitta over a time range
    private :
      // local coordinate system functions, used internally
      VEC3 localDirection(double time, MomBasis::Direction mdir= MomBasis::momdir_) const;
      VEC3 localMomentum(double time) const;
      VEC3 localPosition(double time) const;
      DPDV dPardXLoc(double time) const; // return the derivative of the parameters WRT the local (unrotated) position vector
      DPDV dPardMLoc(double time) const; // return the derivative of the parameters WRT the local (unrotated) momentum vector
      PSMAT dPardStateLoc(double time) const; // derivative of parameters WRT local state
      void setTransforms(); // define global to local and local to global given BNom

      TimeRange trange_;
      Parameters pars_; // parameters
      double mass_;  // in units of MeV/c^2
      int charge_; // charge in units of proton charge
      VEC3 bnom_; // nominal BField, in global coordinate system
      ROOT::Math::Rotation3D l2g_, g2l_; // rotations between local and global coordinates
      const static std::vector<std::string> paramTitles_;
      const static std::vector<std::string> paramNames_;
      const static std::vector<std::string> paramUnits_;
      const static std::string trajName_;
      // non-const accessors
      double& param(size_t index) { return pars_.parameters()[index]; }
  };
  std::ostream& operator <<(std::ostream& ost, LoopHelix const& lhel);
}
#endif
