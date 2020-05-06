#ifndef KinKal_IPHelix_hh
#define KinKal_IPHelix_hh
//
// class desribing the looping helix basis for the kinematic Kalman fit
// It provides geometric, kinematic, and algebraic representation of
// a particule executing a multi-loop helix in a constant magnetic field.
// Original Author Roberto Soleti (LBNL) 1/2020
//

#include "KinKal/Vectors.hh"
#include "KinKal/TRange.hh"
#include "KinKal/PData.hh"
#include "KinKal/LocalBasis.hh"
#include "KinKal/BField.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "Math/Rotation3D.h"
#include <vector>
#include <string>
#include <ostream>

namespace KinKal {

  class IPHelix {
    public:
      // This class must provide the following to be used to instantiate the 
      // classes implementing the Kalman fit
      // define the indices and names of the parameters
      enum ParamIndex {d0_=0,phi0_=1,omega_=2,z0_=3,tanDip_=4,t0_=5,npars_=6};
      constexpr static ParamIndex t0Index() { return t0_; }
      constexpr static size_t NParams() { return npars_; }
      typedef PData<npars_> PDATA; // Data payload for this class
      typedef typename PDATA::DVEC DVEC; // derivative of parameters type
      static std::vector<std::string> const &paramNames();
      static std::vector<std::string> const &paramUnits();
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(ParamIndex index);
      static std::string const& paramUnit(ParamIndex index);
      static std::string const& paramTitle(ParamIndex index);
      static std::string const& trajName();

      // construct from momentum, position, and particle properties.
      // This also requires the BField
      IPHelix(Vec4 const &pos, Mom4 const &mom, int charge, Vec3 const &bnom, TRange const &range = TRange());
      // construct from parameters
      IPHelix(PDATA const &pdata, double mass, int charge, Vec3 const &bnom, TRange const &range = TRange());
      IPHelix(PDATA::DVEC const &pvec, PDATA::DMAT const &pcov, double mass, int charge, Vec3 const &bnom, TRange const &range = TRange());
      // particle position and momentum as a function of time
      void position(Vec4& pos) const ; // time is input
      Vec3 position(double time) const ; // time is input
      Mom4 momentum(double time) const ;
      Vec3 velocity(double time) const ;
      Vec3 direction(double time, LocalBasis::LocDir mdir= LocalBasis::momdir) const;
      // scalar momentum and energy in MeV/c units
      double momentumMag(double time) const  { return mass_ * pbar() / mbar_; }
      double momentumVar(double time) const  { return -1.0; }//FIXME! 
      double energy(double time) const  { return mass_ * ebar() / mbar_; }
      // speed in mm/ns
      double speed(double time) const  { return CLHEP::c_light * beta(); }
      void rangeInTolerance(TRange &range, BField const &bfield, double tol) const ;
      // local momentum direction basis
      void print(std::ostream& ost, int detail) const  {} // FIXME!
      TRange const& range() const { return trange_; }
      TRange& range() { return trange_; }
      void setRange(TRange const& trange) { trange_ = trange; }
      bool inRange(double time) const { return trange_.inRange(time); }

      // momentum change derivatives; this is required to instantiate a KalTrk using this KTraj
      void momDeriv(double time, LocalBasis::LocDir mdir, DVEC &der,Vec3& unit) const;
      double mass() const { return mass_;} // mass 
      int charge() const { return charge_;} // charge in proton charge units

      // named parameter accessors
      double paramVal(size_t index) const { return pars_.parameters()[index]; }
      PDATA const &params() const { return pars_; }
      PDATA &params() { return pars_; }
      double d0() const { return paramVal(d0_); }
      double phi0() const { return paramVal(phi0_); }
      double omega() const { return paramVal(omega_); }
      double z0() const { return paramVal(z0_); }
      double tanDip() const { return paramVal(tanDip_); }
      double t0() const { return paramVal(t0_); }

      // simple functions
      double pbar() const { return 1 / omega() * sqrt( 1 + tanDip() * tanDip() ); } // momentum in mm
      double ebar() const { return sqrt(pbar()*pbar() + mbar_ * mbar_); } // energy in mm
      double cosDip() const { return 1./sqrt(1.+ tanDip() * tanDip() ); }
      double sinDip() const { return tanDip()*cosDip(); }
      double mbar() const { return mbar_; } // mass in mm; includes charge information!
      double vt() const { return vt_; }
      double vz() const { return vz_; }
      double Q() const { return mass_/mbar_; } // reduced charge
      double beta() const { return pbar()/ebar(); } // relativistic beta
      double gamma() const { return fabs(ebar()/mbar_); } // relativistic gamma
      double dphi(double t) const { return omega()*vt()*(t - t0()); }
      double phi(double t) const { return dphi(t) + phi0(); }
      double deltaPhi(double &phi, double refphi=0.) const;
      double angle(const double &f) const;
      double translen(const double &f) const { return cosDip() * f; }
      double arc(const double &f) const { return translen(f) * omega(); }
      double ztime(double zpos) const { return t0() + zpos / vz(); }
      Vec3 const &bnom(double time=0.0) const { return bnom_; }
      double bnomR() const { return bnom_.R(); }
      // flip the helix in time and charge; it remains unchanged geometrically
      void invertCT()
      {
        mbar_ *= -1.0;
        charge_ *= -1;
        pars_.parameters()[t0_] *= -1.0;
      }
      //
    private :
      TRange trange_;
      PDATA pars_; // parameters
      double mass_;  // in units of MeV/c^2
      int charge_; // charge in units of proton charge
      double mbar_;  // reduced mass in units of mm, computed from the mass and nominal field
      Vec3 bnom_;    // nominal BField
      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
      static std::vector<std::string> paramUnits_;
      static std::string trajName_;
      double vt_; // transverse velocity
      double vz_; // z velocity
      // non-const accessors
      double &param(size_t index) { return pars_.parameters()[index]; }
  };
  std::ostream& operator <<(std::ostream& ost, IPHelix const& hhel);
}
#endif
