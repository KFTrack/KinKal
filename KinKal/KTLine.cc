/*
  KTLine is the Linear Trajectory Specialization of KTRAJ - the kinematic trajectory.
  Original Author: S Middleton 2020
*/

#include "KinKal/KTLine.hh"
#include "KinKal/BField.hh"
#include "KinKal/POCAUtil.hh"
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

 KTLine::KTLine( Vec4 const& pos0, Mom4 const& mom0, int charge, double bnom, TRange const& range) : KTLine(pos0,mom0,charge,Vec3(0.0,0.0,bnom),range) {}

  KTLine::KTLine(Vec4 const &pos0, Mom4 const &mom0, int charge, Vec3 const &bnom,
  TRange const &trange) :  bnom_(bnom), mass_(mom0.M()), charge_(charge), trange_(trange)
  {
    pos40_ = pos0;
    mom_ = mom0;
    speed_ = (sqrt(((mom0.Vect() / mom0.E()) * CLHEP::c_light).Mag2()));
    dir_ = ((mom0.Vect() / mom0.E()) * CLHEP::c_light).Unit();  
    Vec4 pos(pos0);
    Mom4 mom(mom0);
    g2l_ = Rotation3D(AxisAngle(Vec3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    if(fabs(g2l_(bnom_).Theta()) > 1.0e-6)throw invalid_argument("Rotation Error");
    pos = g2l_(pos);
    mom = g2l_(mom);
    // create inverse rotation; this moves back into the original coordinate system
    l2g_ = g2l_.Inverse();

    double pt = sqrt(mom.perp2());
    vt_ = CLHEP::c_light * pt / mom.E();
    vz_ = CLHEP::c_light * mom.z() / mom.E();
    static const Vec3 zdir(0.0, 0.0, 1.0);
    static const Vec3 zpos(0.0, 0.0, 0.0);

    const Vec3 pos3(pos0.x(), pos0.y(), pos0.z());
    double zddot = zdir.Dot(dir_);

    POCAUtil *poca = new POCAUtil(pos3, dir_, zpos, zdir);
    Vec3 const& p = poca->point1(); 
    amsign_ = copysign(1.0, p.X());
    param(d0_) = -amsign_*poca->dca();
    param(phi0_) = atan2(-amsign_*p.X(), amsign_*p.Y());
    param(z0_) = p.Z();
    param(cost_) = zddot;
    param(t0_) = pos.T() - (p.Z() - param(z0_)) / (cosTheta() * CLHEP::c_light * beta());
    cout << "In KTLine. Params set to: " << pars_.parameters() << endl;
  }

  /*
  KTLine can take in Momentum externally as a 4-vector or calculate it based. You
  can initialize the line with an origin (pos0) or the trajectory parameters
  (pdata)
  */

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
      break;
    case LocalBasis::phidir: // purely transverse theta2 = -phi()*sin(theta)
      return l2g_(Vec3(-cosPhi0(), sinPhi0(), 0.0));
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
      pder[z0_] = -l * cosTheta(); 
      pder[t0_] = pder[z0_] / vz + 1/tanTheta() * (time - t0()) * sinTheta()*sinTheta()*tanTheta();
      break;
    case LocalBasis::phidir:
      // change in dP/dtheta1 = dP/dphi0*(-1/sintheta)
      pder[cost_] = 0;
      pder[d0_] = l;               
      pder[phi0_] = 1 / sinTheta();
      pder[z0_] = d0() / (sinTheta()*tanTheta()); 
      pder[t0_] = pder[z0_] / vz;
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
