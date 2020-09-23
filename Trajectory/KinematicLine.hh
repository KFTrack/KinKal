#ifndef KinKal_KinematicLine_hh
#define KinKal_KinematicLine_hh
/*
  KinematicLine is the Linear Tracjectory specialization of the KTRAJ
  Original Author: S Middleton 2020
*/
#include "KinKal/General/PhysicalConstants.h"
#include "KinKal/Detector/BFieldMap.hh"
#include "KinKal/General/MomBasis.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/Vectors.hh"
#include "Math/Rotation3D.h"
#include <stdexcept>
#include <vector>
namespace KinKal {

class KinematicLine {
  public:
    enum ParamIndex {
      d0_ = 0,
      phi0_ = 1,
      z0_ = 2,
      cost_ = 3,
      t0_ = 4,
      mom_ = 5,
      npars_ = 6
    };
    constexpr static size_t NParams() { return npars_; }

    static std::vector<std::string> const &paramNames();
    static std::vector<std::string> const &paramUnits();
    static std::vector<std::string> const &paramTitles();
    static std::string const &paramName(ParamIndex index);
    static std::string const &paramTitle(ParamIndex index);
    static std::string const &paramUnit(ParamIndex index);

    constexpr static ParamIndex t0Index() { return t0_; }
    static std::string const &trajName();

    // This also requires the nominal BField, which can be a vector (3d) or a
    // scalar (B along z)
    KinematicLine(VEC4 const& pos, MOM4 const& mom, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange());
    KinematicLine(VEC4 const& pos, MOM4 const& mom, int charge, double bnom, TimeRange const& range=TimeRange());
    // copy payload and adjust for a different BField and range 
    KinematicLine(KinematicLine const& other, VEC3 const& bnom, double trot);

    // copy and override parameters
    KinematicLine(Parameters const &pdata, KinematicLine const& other); 
    KinematicLine(ParticleState const& pstate, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange());
    // same, including covariance information
    KinematicLine(ParticleStateMeasurement const& pstate, int charge, VEC3 const& bnom, TimeRange const& range=TimeRange());

    virtual ~KinematicLine() {}

    // particle momentum as a function of time
    MOM4 momentum4(double time) const;
    VEC3 momentum3(double time) const;

    // scalar momentum and energy in MeV/c units --> Needed for KKTrk:
    double momentum(double time) const { return mom(); }

    double momentumVar(double time) const { return -1.0; }
    double energy() const { double momval = mom(); return sqrt(momval*momval+mass_*mass_); }
    double energy(double time) const { return energy(); }

    // named parameter accessors
    double paramVal(size_t index) const { return pars_.parameters()[index]; }
    Parameters const &params() const { return pars_; }
    Parameters &params() { return pars_; }
    double d0() const { return paramVal(d0_); }
    double phi0() const { return paramVal(phi0_); }
    double z0() const { return paramVal(z0_); }
    double cost() const { return paramVal(cost_); }
    double t0() const { return paramVal(t0_); }
    double mom() const { return paramVal(mom_); }

    // express fit results as a state vector (global coordinates)
    ParticleState state(double time) const;
    ParticleStateMeasurement measurementState(double time) const;

    double translen(const double &f) const { return sinTheta() * f; }
    // simple functions
    double cosTheta() const { return cost(); }
    double sinTheta() const { return sqrt(1.0 - cost() * cost()); }
    double cosPhi0() const { return cos(phi0()); }
    double sinPhi0() const { return sin(phi0()); }
    double theta() const { return acos(cost()); }
    double tanTheta() const { return sqrt(1.0 - cost() * cost()) / cost(); }
    VEC3 pos0() const { return VEC3(-d0()*sinPhi0(),d0()*cosPhi0(),z0()); }
    double flightLength(double t)const { return (t-t0())*speed(); }
    VEC3 dir() const { double st = sinTheta(); return VEC3(st*cosPhi0(),st*sinPhi0(),cosTheta()); }

    TimeRange const &range() const { return trange_; }
    TimeRange &range() {return trange_; }
    virtual void setRange(TimeRange const &trange) { trange_ = trange; }
    void setBNom(double time, VEC3 const& bnom) { bnom_ = bnom; }
    bool inRange(double time) const { return trange_.inRange(time); }

    double speed() const {  return ( mom()/ energy()) * CLHEP::c_light; }
    double speed(double t) const { return speed(); }

    VEC3 position3(double time) const;
    VEC4 position4(double time) const;

    VEC3 velocity(double time) const { return dir() * speed(); }
    void print(std::ostream &ost, int detail) const;
    void rangeInTolerance(TimeRange &range, BFieldMap const &bfield, double tol) const {}; 

    // local momentum direction basis
    VEC3 direction(double time, MomBasis::Direction mdir = MomBasis::momdir_) const;
    // momentum change derivatives; this is required to instantiate a KalTrk
    DVEC momDeriv(double time, MomBasis::Direction mdir) const;

    // some possibly useful equations:
    double mass() const { return mass_; }
    double ztime(double zpos) const {
      return (t0() + zpos / ((speed() * dir()).z()));
    } // time to travel Z

    int charge() const { return charge_; }
    double beta() const { return (speed() / CLHEP::c_light); }
    double gamma() const { return energy()/mass_; }
    double betaGamma() const { return beta() * gamma(); }
    VEC3 const &bnom(double time = 0.0) const { return bnom_; }

    DPDV dPardX(double time) const;
    DPDV dPardM(double time) const;
    DVDP dXdPar(double time) const;
    DVDP dMdPar(double time) const;
    PSMAT dPardState(double time) const;
    PSMAT dStatedPar(double time) const;
    // package the above for full (global) state
    // Parameter derivatives given a change in BField.  These return null for KinematicLine
    DVEC dPardB(double time) const { return DVEC(); }
    DVEC dPardB(double time, VEC3 const& BPrime) const { return DVEC(); }

    void invertCT() {
      charge_ *= -1;
      pars_.parameters()[t0_] *= -1.0;
      // need to invert direction vector too FIXME!
    }

  private:
    const static std::string trajName_;
  // non-parametric variables
    VEC3 bnom_; // nominal BField: not used by this parameterization 
    double mass_;   // mass in MeV/c2
    int charge_; // charge in proton charge unites
    TimeRange trange_;  // valid range
    Parameters pars_;      // parameters
    const static std::vector<std::string> paramTitles_;
    const static std::vector<std::string> paramNames_;
    const static std::vector<std::string> paramUnits_;

    // nonconst accessors
    double &param(size_t index) { return pars_.parameters()[index]; }
  };
  std::ostream &operator<<(std::ostream &ost, KinematicLine const &line);

} // namespace KinKal
#endif
