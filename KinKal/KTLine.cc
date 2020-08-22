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
      "Z at POCA (z_{0})", "Cos #theta", "Time at POCA (t_{0})", "Momentum"};

  vector<string> KTLine::paramNames_ = {"d_{0}", "#phi_{0}", "z_{0}",
                                       "cos(#theta)", "t_{0}", "mom"};

  vector<string> KTLine::paramUnits_ = {"mm", "radians", "mm", "", "ns","MeV/c"};

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
    double mommag = mom0.R();
    double speed = ( mommag/ mom0.E()) * CLHEP::c_light;
    Vec3 dir = mom0.Vect().Unit();
    
    static const Vec3 zdir(0.0, 0.0, 1.0);
    static const Vec3 zpos(0.0, 0.0, 0.0);

// calculate POCA to the Z axis.  This is the reference for the parameters
    POCAUtil poca(pos0.Vect(), dir, zpos, zdir);
    double flen = dir.Dot(pos0.Vect()-poca.point1()); // flight length from reference to POCA
    Vec3 pca = poca.point1()-poca.point2(); // vector from Z axis to POCA
// doca is signed by the angular momentum around the Z axis
    double amsign = copysign(1.0, -(zdir.Cross(pca)).Dot(dir));
    param(d0_) = amsign*poca.dca(); // dca2d and dca are the same for POCA to the Z axis 
    param(phi0_) = dir.Phi(); // same as at POCA
    param(z0_) = poca.point1().Z();
    param(cost_) = dir.Z();
    param(t0_) = pos0.T() - flen/speed;
    param(mom_) = mommag;
    cout << "In KTLine. poca=: " << poca.point1() << " " << poca.point2() << endl;
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
    return (pos0() + flightLength(time) * dir());
  }

  Vec4 KTLine::pos4(double time) const {
    Vec3 temp = position(time);
    return Vec4(temp.X(), temp.Y(), temp.Z(), time);
  }

  void KTLine::momentum(double tval, Mom4 &mom4) const {
    Vec3 dir = direction(tval);
    double momval = mom();
    mom4.SetPx(momval * dir.x());
    mom4.SetPy(momval * dir.y());
    mom4.SetPz(momval * dir.z());
    mom4.SetM(mass_);
  }

  Mom4 KTLine::momentum(double tval) const {
    Mom4 momvec;
    momentum(tval,momvec);
    return momvec;
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
      return Vec3(cosTheta() * sinPhi0(), cosTheta() * cosPhi0(), -1 * sinTheta());
      break;
    case LocalBasis::phidir: // purely transverse theta2 = -phi()*sin(theta)
      return Vec3(-cosPhi0(), sinPhi0(), 0.0);
      break;
    case LocalBasis::momdir: // along momentum
      return dir();
      break;
    default:
      throw std::invalid_argument("Invalid direction");
    }
  }

  // derivatives of momentum projected along the given basis WRT the 5 parameters
  KTLine::DVEC KTLine::momDeriv(double time, LocalBasis::LocDir mdir) const {
    double momval = mom();
    double gv = gamma();
    double vz = speed()*dir().Z(); 
    double deltat = time-t0();
    double tlen = translen(CLHEP::c_light * beta() *deltat );
    KTLine::DVEC pder;
    double sinT = sinTheta();
    double tanT = tanTheta();
    double cosT = cosTheta();
    // cases
    switch (mdir) {
    case LocalBasis::perpdir:
      // polar bending: change in Theta
      pder[cost_] = -sinT; 
      pder[d0_] = 0;
      pder[phi0_] = 0;
      pder[z0_] = tlen*cosT/sinT; //TODO - I think here is the issue!
      pder[t0_] = pder[z0_] / vz + 1/tanTheta() * (time - t0()) * sinT*sinT*tanT;
      pder[mom_] = 0.0;
      break;
    case LocalBasis::phidir:
      // change in dP/dtheta1 = dP/dphi0*(-1/sintheta)
      pder[cost_] = 0;
      pder[d0_] = -tlen/sinT;            
      pder[phi0_] = 1.0/ sinT;
      pder[z0_] = 0;
      pder[t0_] = pder[z0_] / vz;
      pder[mom_] = 0.0;
      break;
    case LocalBasis::momdir:
      pder[cost_] = 0;
      pder[d0_] = 0;
      pder[phi0_] = 0;
      pder[z0_] = 0;
      pder[t0_] = deltat/(gv*gv);
      pder[mom_] = 1.0*momval;
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
