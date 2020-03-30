#ifndef KinKal_IPHelix_hh
#define KinKal_IPHelix_hh
//
// class desribing the looping helix basis for the kinematic Kalman fit
// It provides geometric, kinematic, and algebraic representation of
// a particule executing a multi-loop helix in a constant magnetic field.
// Original Author David Brown (LBNL) 1/2020
//

#include "KinKal/Vectors.hh"
#include "KinKal/TTraj.hh"
#include "KinKal/KInter.hh"
#include "KinKal/PData.hh"
#include "KinKal/BField.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include <vector>
#include <string>
#include <ostream>

namespace KinKal {

  class IPHelix : public TTraj, public KInter {
    public:
      // This class must provide the following to be used to instantiate the 
      // classes implementing the Kalman fit
      // define the indices and names of the parameters
      enum paramIndex {d0_=0,phi0_=1,omega_=2,z0_=3,tanDip_=4,t0_=5,npars_=6};
      constexpr static size_t NParams() { return npars_; }
      typedef PData<npars_> PDATA; // Data payload for this class
      static std::vector<std::string> const &paramNames();
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(paramIndex index);
      static std::string const& paramTitle(paramIndex index);

      // construct from momentum, position, and particle properties.
      // This also requires the BField
      IPHelix(Vec4 const &pos, Mom4 const &mom, int charge, Vec3 const &bnom, TRange const &range = TRange());
      // construct from parameters
      IPHelix(PDATA const &pdata, double mass, int charge, Vec3 const &bnom, TRange const &range = TRange());
      IPHelix(PDATA::DVEC const &pvec, PDATA::DMAT const &pcov, double mass, int charge, Vec3 const &bnom, TRange const &range = TRange());
      // particle position and momentum as a function of time
      void position(Vec4& pos) const override; // time is input
      void position(double t,Vec3& pos) const override; // time is input
      void momentum(double t,Mom4& mom) const override;
      void velocity(double time, Vec3& vel) const override;
      void direction(double tval,Vec3& dir) const override;
      // scalar momentum and energy in MeV/c units
      double momentum(double time) const override { return mass_ * pbar() / mbar_; }
      double energy(double time) const override { return mass_ * ebar() / mbar_; }
      // speed in mm/ns
      double speed(double time) const override { return CLHEP::c_light * beta(); }
      void rangeInTolerance(TRange &range, BField const &bfield, double tol) const override;
      // local momentum direction basis
      virtual void dirVector(MDir dir, double time, Vec3 &unit) const override;

      // momentum change derivatives; this is required to instantiate a KalTrk using this KTraj
      typedef ROOT::Math::SMatrix<double, npars_, 1> PDER; // derivative of parameters type in a particular direction
      void momDeriv(MDir dir, double time, PDER &der) const;

      // named parameter accessors
      double param(size_t index) const { return pars_.parameters()[index]; }
      PDATA const &params() const { return pars_; }

      double d0() const { return param(d0_); }
      double phi0() const { return param(phi0_); }
      double omega() const { return param(omega_); }
      double z0() const { return param(z0_); }
      double tanDip() const { return param(tanDip_); }
      double t0() const { return param(t0_); }

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
      Vec3 const &bnom() const { return bnom_; }
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
      PDATA pars_; // parameters
      double mbar_;  // reduced mass in units of mm, computed from the mass and nominal field
      Vec3 bnom_;    // nominal BField
      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
      double vt_; // transverse velocity
      double vz_; // z velocity
      // non-const accessors
      double &param(size_t index) { return pars_.parameters()[index]; }
  };
  std::ostream& operator <<(std::ostream& ost, IPHelix const& hhel);
}
#endif
