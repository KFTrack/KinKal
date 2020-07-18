/*
  KTLine is the Linear Trajectory Specialization of KTRAJ - the kinematic trajectory.
  Original Author: S Middleton 2020
*/

#include "KinKal/KTLine.hh"
#include "KinKal/BField.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  vector<string> KTLine::paramTitles_ = {
      "Transverse DOCA to Z Axis (d_{0})", "Azimuth of POCA (#phi_{0})",
      "Z at POCA (z_{0})", "Cos #theta", "Time at POCA (t_{0})"};

  vector<string> KTLine::paramNames_ = {"d_{0}", "#phi_{0}", "z_{0}",
                                       "cos(#theta)", "t_{0}"};

  vector<string> KTLine::paramUnits_ = {"mm", "radians", "mm", "", "ns"};

  std::vector<std::string> const &KTLine::paramUnits() { return paramUnits_; }
  std::vector<std::string> const &KTLine::paramNames() { return paramNames_; }
  std::vector<std::string> const &KTLine::paramTitles() { return paramTitles_; }

  std::string const &KTLine::paramName(ParamIndex index) {
    return paramNames_[static_cast<size_t>(index)];
  }
  std::string const &KTLine::paramTitle(ParamIndex index) {
    return paramTitles_[static_cast<size_t>(index)];
  }
  std::string const &KTLine::paramUnit(ParamIndex index) {
    return paramUnits_[static_cast<size_t>(index)];
  }

  string KTLine::trajName_("KTLine");
  string const &KTLine::trajName() { return trajName_; }

  KTLine::KTLine(Vec4 const &pos0, Vec3 const &svel, TRange const &range, bool forcerange)
      : KTLine(pos0.Vect(), svel, pos0.T(), range, forcerange) {
    std::cout << " KTLine Constructor 1 " << trange_ << std::endl;
  }

  KTLine::KTLine(Vec3 const &pos0, Vec3 const &svel, double tmeas, TRange const &range, bool forcerange)
      : trange_(range), pars_(), speed_(sqrt(svel.Mag2())), pos0_(pos0),
        dir_(svel.Unit()), forcerange_(forcerange) {
   
    static const Vec3 zdir(0.0, 0.0, 1.0);
    double zddot = zdir.Dot(dir_);
    param(cost_) = zddot;
    param(d0_) = pos0_.Rho();
    param(phi0_) = atan2(pos0_.Y(), pos0_.X());
    param(z0_) = pos0_.Z();
    param(t0_) = tmeas;
    cout << "In KTLine. Params set to: " << pars_.parameters() << endl;
  }

  /*
  KTLine can take in Momentum externally as a 4-vector or calculate it based. You
  can initialize the line with an origin (pos0) or the trajectory parameters
  (pdata)
  */

  KTLine::KTLine(Vec4 const &pos0, Mom4 const &mom0, int charge, double bnom,TRange const &range)
      : KTLine(pos0, mom0, charge, Vec3(0.0, 0.0, bnom), range) {}

  KTLine::KTLine(Vec4 const &pos40, Mom4 const &mom0, int charge, Vec3 const &bnom, TRange const &range)
      : KTLine(pos40.Vect(), (mom0.Vect() / mom0.E()) * CLHEP::c_light, pos40.T(),range) {
        bnom_ = bnom;
        pos40_ = pos40;
        mom_ = mom0;
        charge_ = charge;
        mass_ = mom0.M();
  }

  KTLine::KTLine(PDATA const &pdata, double mass, int charge, double bnom,TRange const &range)
      : KTLine(pdata, mass, charge, Vec3(0.0, 0.0, bnom), range) {
  }

  KTLine::KTLine(PDATA const &pdata, double mass, int charge, Vec3 const &bnom, TRange const &range)
      : KTLine(pdata.parameters(), pdata.covariance(), mass, charge, bnom, range) {}

  KTLine::KTLine(PDATA::DVEC const &pvec, PDATA::DMAT const &pcov, double mass, int charge, Vec3 const &bnom, TRange const &trange) : KTLine(pvec, pcov) {
    bnom_ = bnom;
    mass_ = mass;
    charge_ = charge;
  }

  KTLine::KTLine(PDATA const &pdata, KTLine const &other) : KTLine(other) {
    pars_ = pdata;
  }

  void KTLine::position(Vec4 &pos) const {
    Vec3 pos3 = position(pos.T());
    pos.SetXYZT(pos3.X(), pos3.Y(), pos3.Z(), pos.T());
  }

  Vec3 KTLine::position(double time) const {
    if (forcerange_){
      range().forceRange(time);
    }
    return (pos0() + ((time - t0()) * speed()) * dir());
  }

  Vec4 KTLine::pos4(double time) const {
    Vec3 temp = position(time);
    return Vec4(temp.X(), temp.Y(), temp.Z(), time);
  }

  void KTLine::momentum(double tval, Mom4 &mom) const {
    Vec3 dir = direction(tval);
    mom.SetPx(momentumMag(tval) * dir.x());
    mom.SetPy(momentumMag(tval) * dir.y());
    mom.SetPz(momentumMag(tval) * dir.z());
    mom.SetM(mass_);
  }

  Mom4 KTLine::momentum(double tval) const {
    Mom4 mom;
    Vec3 dir = direction(tval);
    mom.SetPx(momentumMag(tval) * dir.x());
    mom.SetPy(momentumMag(tval) * dir.y());
    mom.SetPz(momentumMag(tval) * dir.z());
    mom.SetM(mass_);
    return mom_;
  }

  /*
  The effects for changes in 2 perpendicular directions (theta1 = theta and
  theta2 = phi()*sin(theta) can sometimes be added, as scattering in these
  are uncorrelated. These axes are track specific. as cosmics are not always
  coming along the same track direction it is necessary to have difference
  parameterization than that used for the helix case.
  alt dir = a test with the "BTrk parameterization" - just changes signs due to
  swithc in cos<->sin
  */
  Vec3 KTLine::direction(double t, LocalBasis::LocDir mdir) const {

    switch (mdir) {
    case LocalBasis::perpdir: // purely polar change theta 1 = theta
      return l2g_(
          Vec3(cosTheta() * sinPhi0(), cosTheta() * cosPhi0(), -1 * sinTheta()));
      // alt dir
      // l2g_(Vec3(cosTheta()*cosPhi0(),cosTheta()*sinPhi0(),-1*sinTheta()));
      break;
    case LocalBasis::phidir: // purely transverse theta2 = -phi()*sin(theta)
      return l2g_(Vec3(-cosPhi0(), sinPhi0(), 0.0));
      // alt dir = l2g_(Vec3(-sinPhi0(),cosPhi0(),0.0));//;
      break;
    case LocalBasis::momdir: // along momentum: check.
      return l2g_(Vec3(dir().x(), dir().y(), dir().z()));
      break;
    default:
      throw std::invalid_argument("Invalid direction");
    }
  }

  // derivatives of momentum projected along the given basis WRT the 5 parameters
  KTLine::DVEC KTLine::momDeriv(double time, LocalBasis::LocDir mdir) const {
    // compute some useful quantities
    double vz = CLHEP::c_light * mom().z() / mom().E();
    double l = speed() * time;
    KTLine::DVEC pder;
    //cout << "Mom deriv start params " << pder << endl;
    // cases
    switch (mdir) {
    case LocalBasis::perpdir:
      // polar bending: change in Theta
      pder[cost_] = -sinTheta();
      pder[d0_] = 0;
      pder[phi0_] = 0;
      pder[z0_] = l * cosTheta(); // alt dir =-l*cosTheta();
      pder[t0_] = pder[z0_] / vz;
      //cout << "Mom deriv perpdir params " << pder << endl;
      break;
    case LocalBasis::phidir:
      // change in dP/dtheta1 = dP/dphi0*(-1/sintheta)
      pder[cost_] = 0;
      pder[d0_] = l;                 // alt dir = -l;
      pder[phi0_] = -1 / sinTheta(); // alt dir = -1/sinTheta();
      pder[z0_] = d0() / (sinTheta() *
                          tanTheta()); // alt dir = -d0()/(sinTheta()*tanTheta());
      pder[t0_] = pder[z0_] / vz;
      //cout << "Mom deriv phidir params " << pder << endl;
      break;
    case LocalBasis::momdir:
      pder[cost_] = 0;
      pder[d0_] = 0;
      pder[phi0_] = 0;
      pder[z0_] = 0;
      pder[t0_] = pder[z0_] / vz;
      //cout << "Mom deriv momdir params " << pder << endl;
      break;

    default:
      throw std::invalid_argument("Invalid direction");
    }
    return pder;
  }

  void KTLine::print(ostream &ost, int detail) const {
    auto perr = params().diagonal();
    ost << " KTLine " << range() << " parameters: ";
    for (size_t ipar = 0; ipar < KTLine::npars_; ipar++) {
      ost << KTLine::paramName(static_cast<ParamIndex>(ipar)) << " "
          << paramVal(ipar) << " +- " << perr(ipar);
      if (ipar < KTLine::npars_ - 1)
        ost << " ";
    }
    ost << " with rotation around Bnom " << bnom_ << endl;
  }

  ostream &operator<<(ostream &ost, KTLine const &lhel) {
    lhel.print(ost, 0);
    return ost;
  }

} // namespace KinKal
