#ifndef KinKal_KTLine_hh
#define KinKal_KTLine_hh
/*
  KTLine is the Linear Tracjectory specialization of the KTRAJ
  Original Author: S Middleton 2020
*/
#include "CLHEP/Units/PhysicalConstants.h"
#include "KinKal/BField.hh"
#include "KinKal/LocalBasis.hh"
#include "KinKal/PData.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TRange.hh"
#include "KinKal/Vectors.hh"
#include "Math/Rotation3D.h"
#include <stdexcept>
#include <vector>
namespace KinKal {

class KTLine {
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
    typedef PData<npars_> PDATA;       // Data payload for this class
    typedef typename PDATA::DVEC DVEC; // derivative of parameters type

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
    KTLine(Vec4 const& pos, Mom4 const& mom, int charge, Vec3 const& bnom, TRange const& range=TRange());
    KTLine(Vec4 const& pos, Mom4 const& mom, int charge, double bnom, TRange const& range=TRange());
    // copy and override parameters
    KTLine(PDATA const &pdata, KTLine const& other); 

    virtual ~KTLine() {}

    // particle momentum as a function of time
    Mom4 momentum(double time) const;
    void momentum(double t, Mom4 &mom) const;

    // scalar momentum and energy in MeV/c units --> Needed for KKTrk:
    double momentumMag(double time) const { return mom(); }

    double momentumVar(double time) const { return -1.0; }
    double energy() const { double momval = mom(); return sqrt(momval*momval+mass_*mass_); }
    double energy(double time) const { return energy(); }

    // named parameter accessors
    double paramVal(size_t index) const { return pars_.parameters()[index]; }
    PDATA const &params() const { return pars_; }
    PDATA &params() { return pars_; }
    double d0() const { return paramVal(d0_); }
    double phi0() const { return paramVal(phi0_); }
    double z0() const { return paramVal(z0_); }
    double cost() const { return paramVal(cost_); }
    double t0() const { return paramVal(t0_); }
    double mom() const { return paramVal(mom_); }
    double translen(const double &f) const { return sinTheta() * f; }
    // simple functions
    double cosTheta() const { return cost(); }
    double sinTheta() const { return sqrt(1.0 - cost() * cost()); }
    double cosPhi0() const { return cos(phi0()); }
    double sinPhi0() const { return sin(phi0()); }
    double theta() const { return acos(cost()); }
    double tanTheta() const { return sqrt(1.0 - cost() * cost()) / cost(); }
    Vec3 pos0() const { return Vec3(-d0()*sinPhi0(),d0()*cosPhi0(),z0()); }
    double flightLength(double t)const { return (t-t0())*speed(); }
    Vec3 dir() const { double st = sinTheta(); return Vec3(st*cosPhi0(),st*sinPhi0(),cosTheta()); }

    TRange const &range() const { return trange_; }
    TRange &range() {return trange_; }
    virtual void setRange(TRange const &trange) { trange_ = trange; }
    bool inRange(double time) const { return trange_.inRange(time); }

    double speed() const {  return ( mom()/ energy()) * CLHEP::c_light; }
    double speed(double t) const { return speed(); }

    void position(Vec4 &pos) const;
    Vec3 position(double time) const;

    Vec3 velocity(double time) const { return dir() * speed(); }
    void print(std::ostream &ost, int detail) const;
    void rangeInTolerance(TRange &range, BField const &bfield, double tol) const {}; 

    // local momentum direction basis
    Vec3 direction(double time, LocalBasis::LocDir mdir = LocalBasis::momdir) const;
    Vec4 pos4(double time) const;
    // momentum change derivatives; this is required to instantiate a KalTrk
    KTLine::DVEC momDeriv(double time, LocalBasis::LocDir mdir) const;

    // some possibly useful equations:
    double mass() const { return mass_; }
    double ztime(double zpos) const {
      return (t0() + zpos / ((speed() * dir()).z()));
    } // time to travel Z

    int charge() const { return charge_; }
    double beta() const { return (speed() / CLHEP::c_light); }
    double gamma() const { return energy()/mass_; }
    double betaGamma() const { return beta() * gamma(); }
    Vec3 const &bnom(double time = 0.0) const { return bnom_; }
    void invertCT() {
      charge_ *= -1;
      pars_.parameters()[t0_] *= -1.0;
      // need to invert direction vector too FIXME!
    }

  private:
    static std::string trajName_;
  // non-parametric variables
    Vec3 bnom_;
    double mass_;                 // mass in MeV/c2
    int charge_;
    TRange trange_;
    PDATA pars_;      // parameters
    static std::vector<std::string> paramTitles_;
    static std::vector<std::string> paramNames_;
    static std::vector<std::string> paramUnits_;

    // nonconst accessors
    double &param(size_t index) { return pars_.parameters()[index]; }
  };
  std::ostream &operator<<(std::ostream &ost, KTLine const &line);

} // namespace KinKal
#endif
