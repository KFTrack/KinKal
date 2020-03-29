#ifndef KinKal_LHelix_hh
#define KinKal_LHelix_hh
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
#include "CLHEP/Units/PhysicalConstants.h"
#include "Math/Rotation3D.h"
#include <vector>
#include <string>
#include <ostream>

namespace KinKal {

  class LHelix : public TTraj, public KInter {
    public:
      // This class must provide the following to be used to instantiate the 
      // classes implementing the Kalman fit
      // define the indices and names of the parameters
      enum ParamIndex {rad_=0,lam_=1,cx_=2,cy_=3,phi0_=4,t0_=5,npars_=6};
      constexpr static ParamIndex t0Index() { return t0_; }
      constexpr static size_t NParams() { return npars_; }
      typedef PData<npars_> PDATA; // Data payload for this class
      typedef ROOT::Math::SVector<double,npars_> PDER; // derivative of parameters type 
      static std::vector<std::string> const& paramNames(); 
      static std::vector<std::string> const& paramUnits(); 
      static std::vector<std::string> const& paramTitles();
      static std::string const& paramName(ParamIndex index);
      static std::string const& paramUnit(ParamIndex index);
      static std::string const& paramTitle(ParamIndex index);

      // construct from momentum, position, and particle properties.
      // This also requires the nominal BField, which can be a vector (3d) or a scalar (B along z)
      LHelix(Vec4 const& pos, Mom4 const& mom, int charge, Vec3 const& bnom, TRange const& range=TRange());
      LHelix(Vec4 const& pos, Mom4 const& mom, int charge, double bnom, TRange const& range=TRange());
      // construct from parameters
      LHelix(PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range=TRange());
      LHelix(PDATA const& pdata, double mass, int charge, double bnom, TRange const& range=TRange());
      virtual ~LHelix() {} 
      // particle position and momentum as a function of time
      void position(Vec4& pos) const override; // time is input 
      void position(double t,Vec3& pos) const override; // time is input 
      void velocity(double time, Vec3& vel) const override;
      void direction(double tval,Vec3& dir) const override;
      void momentum(double t,Mom4& mom) const override;
      // scalar momentum and energy in MeV/c units
      double momentum(double time) const override { return  mass_*pbar()/mbar_; }
      double energy(double time) const override { return  mass_*ebar()/mbar_; }
      // speed in mm/ns
      double speed(double time) const override {  return CLHEP::c_light*beta(); }
      void rangeInTolerance(TRange& range, BField const& bfield, double tol) const override;

      // local momentum direction basis
      virtual void dirVector(MDir dir,double time,Vec3& unit) const override;

      // momentum change derivatives; this is required to instantiate a KalTrk using this KInter
      void momDeriv(MDir mdir, double time, PDER& der) const;
      // position change derivatives; this is required to instantiate a KalTrk using this KInter
      // we only care about this along the momentum direction
      void posDeriv(double time, PDER& der) const;

     // named parameter accessors
      double param(size_t index) const { return pars_.parameters()[index]; }
      PDATA const& params() const { return pars_; }
      PDATA& params() { return pars_; }
      double rad() const { return param(rad_); }
      double lam() const { return param(lam_); }
      double cx() const { return param(cx_); }
      double cy() const { return param(cy_); }
      double phi0() const { return param(phi0_); }
      double t0() const { return param(t0_); }
      
      // simple functions; these can be cached if they cause performance problems
      double pbar2() const { return  rad()*rad() + lam()*lam(); } 
      double pbar() const { return  sqrt(pbar2()); } // momentum in mm
      double ebar2() const { return  rad()*rad() + lam()*lam() + mbar_*mbar_; }
      double ebar() const { return  sqrt(ebar2()); } // energy in mm
      double mbar() const { return mbar_; } // mass in mm; includes charge information!
      double Q() const { return mbar_/mass_; } // reduced charge
      double omega() const { return copysign(CLHEP::c_light,mbar_)/ebar(); } // rotational velocity, sign set by magnetic force 
      double beta() const { return pbar()/ebar(); } // relativistic beta
      double gamma() const { return fabs(ebar()/mbar_); } // relativistic gamma
      double dphi(double t) const { return omega()*(t - t0()); }
      double phi(double t) const { return dphi(t) + phi0(); }
      double ztime(double zpos) const { return t0() + zpos/(omega()*lam()); }
      double zphi(double zpos) const { return zpos/lam() + phi0(); }
      int charge() const { return charge_; }
      Vec3 const& bnom() const { return bnom_; }
      double bnomR() const { return bnom_.R(); }
      // flip the helix in time and charge; it remains unchanged geometrically
      void invertCT() {
	mbar_ *= -1.0;
	charge_ *= -1;
	pars_.parameters()[t0_] *= -1.0;
      }
      //
    private :
      PDATA pars_; // parameters
      double mbar_;  // reduced mass in units of mm, computed from the mass and nominal field
      Vec3 bnom_; // nominal BField
      bool needsrot_; // logical flag if Bnom is parallel to global Z or not
      ROOT::Math::Rotation3D brot_; // rotation from the internal coordinate system (along B) to the global
      static std::vector<std::string> paramTitles_;
      static std::vector<std::string> paramNames_;
      static std::vector<std::string> paramUnits_;
      // non-const accessors
      double& param(size_t index) { return pars_.parameters()[index]; }
 };
  std::ostream& operator <<(std::ostream& ost, LHelix const& lhel);
}
#endif
