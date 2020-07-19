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
      npars_ = 5
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
    double momentumMag(double time) const {
      return gamma() * mass_ * beta();
    } 

    Mom4 mom() const { return mom_; }
    double momentumVar(double time) const { return -1.0; }
    double energy(double time) { return mom_.E(); }

    // named parameter accessors
    double paramVal(size_t index) const { return pars_.parameters()[index]; }
    PDATA const &params() const { return pars_; }
    PDATA &params() { return pars_; }
    double d0() const { return paramVal(d0_); }
    double phi0() const { return paramVal(phi0_); }
    double z0() const { return paramVal(z0_); }
    double cost() const { return paramVal(cost_); }
    double t0() const { return paramVal(t0_); }
    double translen(const double &f) const { return sinTheta() * f; }
    // simple functions
    double cosTheta() const { return cost(); }
    double sinTheta() const { return sqrt(1.0 - cost() * cost()); }
    double cosPhi0() const { return cos(phi0()); }
    double sinPhi0() const { return sin(phi0()); }
    double theta() const { return acos(cost()); }
    double tanTheta() const { return sqrt(1.0 - cost() * cost()) / cost(); }

    Vec3 const &dir() const { return dir_; }

    TRange const &range() const { return trange_; }
    TRange &range() {return trange_; }
    virtual void setRange(TRange const &trange) { trange_ = trange; }
    bool inRange(double time) const { return trange_.inRange(time); }

    // access position and direction
    Vec3 const &pos0() const { return pos0_; }
    double speed() const { return speed_; }
    double speed(double t) const { return speed_; }

    void position(Vec4 &pos) const;
    Vec3 position(double time) const;

    Vec3 velocity(double time) const { return dir_ * speed(); }
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
    double gamma() const {
      return (1 / sqrt(1 - ((speed() / CLHEP::c_light) *
                            (speed() / CLHEP::c_light))));
    }
    double betaGamma() const { return beta() * gamma(); }

    double energyBG(double time) const {
      return (sqrt(mass_ * mass_ + betaGamma() * betaGamma() * mass_ * mass_));
    } // in MeV
    Vec3 const &bnom(double time = 0.0) const { return bnom_; }

    void invertCT() {
      charge_ *= -1;
      pars_.parameters()[t0_] *= -1.0;
    }

  private:
    Vec3 bnom_;     // should be 0,0,0
    bool needsrot_; // logical flag if Bnom is parallel to global Z or not
    ROOT::Math::Rotation3D brot_; // rotation from the internal coordinate system
                // (along B) to the global
    Mom4 pos40_, mom_;            // 4 momentum vector - px,py,pz,m
    double mass_;                 // mass in MeV/c2
    int charge_;
    ROOT::Math::Rotation3D l2g_,
    g2l_; // rotations between local and global coordinates
    static std::string trajName_;

    TRange trange_;
    PDATA pars_;      // parameters
    double speed_;    // signed linear velocity, translates time to distance along
    Vec3 pos0_, dir_; // caches
    bool forcerange_; // if set, strictly enforce the range
    double vt_; // transverse velocity
    double vz_; // z velocity
    double amsign_;
    static std::vector<std::string> paramTitles_;
    static std::vector<std::string> paramNames_;
    static std::vector<std::string> paramUnits_;

    // nonconst accessors
    double &param(size_t index) { return pars_.parameters()[index]; }
  };
  std::ostream &operator<<(std::ostream &ost, KTLine const &line);

} // namespace KinKal
#endif
