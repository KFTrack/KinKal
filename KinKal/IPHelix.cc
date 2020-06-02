#include "KinKal/IPHelix.hh"
#include "KinKal/BField.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  vector<string> IPHelix::paramTitles_ = {
    "Distance of closest approach d_{0}",
    "Angle in the xy plane at closest approach #phi_{0}",
    "xy plane curvature of the track #omega",
    "Distance from the closest approach to the origin z_{0}",
    "Tangent of the track dip angle in the #rho - z projection tan#lambda",
    "Time at Z=0 Plane"};
  vector<string> IPHelix::paramNames_ = {
  "d0","phi0","omega","z0","tanDip","Time0"};
  vector<string> IPHelix::paramUnits_ = {
      "mm", "rad", "rad", "mm", "", "ns"};
  string IPHelix::trajName_("IPHelix");  
  std::vector<std::string> const& IPHelix::paramNames() { return paramNames_; }
  std::vector<std::string> const& IPHelix::paramUnits() { return paramUnits_; }
  std::vector<std::string> const& IPHelix::paramTitles() { return paramTitles_; }
  std::string const& IPHelix::paramName(ParamIndex index) { return paramNames_[static_cast<size_t>(index)];}
  std::string const& IPHelix::paramUnit(ParamIndex index) { return paramUnits_[static_cast<size_t>(index)]; }
  std::string const& IPHelix::paramTitle(ParamIndex index) { return paramTitles_[static_cast<size_t>(index)];}
  string const& IPHelix::trajName() { return trajName_; }

  IPHelix::IPHelix( Vec4 const& pos0, Mom4 const& mom0, int charge, double bnom, TRange const& range) : IPHelix(pos0,mom0,charge,Vec3(0.0,0.0,bnom),range) {}
  IPHelix::IPHelix(Vec4 const &pos0, Mom4 const &mom0, int charge, Vec3 const &bnom,
                 TRange const &trange) : trange_(trange), mass_(mom0.M()), charge_(charge), bnom_(bnom)
  {
    // Transform into the system where Z is along the Bfield.  This is a pure rotation about the origin
    Vec4 pos(pos0);
    Mom4 mom(mom0);
    g2l_ = Rotation3D(AxisAngle(Vec3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    if(fabs(g2l_(bnom_).Theta()) > 1.0e-6)throw invalid_argument("Rotation Error");
    pos = g2l_(pos);
    mom = g2l_(mom);
    // create inverse rotation; this moves back into the original coordinate system
    l2g_ = g2l_.Inverse();
    double momToRad = 1.0/(BField::cbar()*charge_*bnom_.R());
    mbar_ = -mass_ * momToRad;

    double pt = sqrt(mom.perp2());
    double radius = fabs(pt*momToRad);

    double lambda = -mom.z()*momToRad;
    double amsign = copysign(1.0, mbar_);

    Vec3 center = Vec3(pos.x() + mom.y()*momToRad, pos.y() - mom.x()*momToRad, 0.0);
    double rcent = sqrt(center.perp2());
    double fcent = center.phi();
    double centerx = rcent*cos(fcent);
    double centery = rcent*sin(fcent);

    param(omega_) = amsign/radius;
    param(tanDip_) = amsign*lambda/radius;
    param(d0_) = amsign*(rcent - radius);
    param(phi0_) = atan2(-amsign * centerx, amsign * centery);

    Vec3 pos3 = Vec3(pos.x(), pos.y(), pos.z());
    double fz0 = (pos3 - center).phi() - pos3.z() / lambda;
    deltaPhi(fz0);
    double refphi = fz0+amsign*M_PI_2;
    double phival = phi0();
    double dphi = deltaPhi(phival, refphi);
    param(z0_) = dphi * tanDip() / omega();

    param(t0_) = pos.T() - (pos.Z() - param(z0_)) / (sinDip() * CLHEP::c_light * beta());
    vt_ = CLHEP::c_light * pt / mom.E();
    vz_ = CLHEP::c_light * mom.z() / mom.E();
    // test position and momentum function
    Vec4 testpos(pos0);
    // std::cout << "Testpos " << testpos << std::endl;
    position(testpos);
    Mom4 testmom = momentum(testpos.T());
    auto dp = testpos.Vect() - pos0.Vect();
    auto dm = testmom.Vect() - mom0.Vect();
    if(dp.R() > 1.0e-5 || dm.R() > 1.0e-5)throw invalid_argument("Rotation Error");
  }

  double IPHelix::deltaPhi(double &phi, double refphi) const
  {
    double dphi = phi - refphi;
    static const double twopi = 2 * M_PI;
    while (dphi > M_PI)
    {
      dphi -= twopi;
      phi -= twopi;
    }
    while (dphi <= -M_PI)
    {
      dphi += twopi;
      phi += twopi;
    }
    return dphi;
  }

  IPHelix::IPHelix(PDATA const &pdata, IPHelix const& other) : IPHelix(other) {
    pars_ = pdata;
  }

  void IPHelix::position(Vec4 &pos) const
 {
    Vec3 pos3 = position(pos.T());
    pos.SetPx(pos3.X());
    pos.SetPy(pos3.Y());
    pos.SetPz(pos3.Z());
  }

  Vec4 IPHelix::pos4(double time) const
 {
    Vec3 pos3 = position(time);
    return Vec4( pos3.X(), pos3.Y(), pos3.Z(), time);
  }

  Vec3 IPHelix::position(double time) const
  {
    double cDip = cosDip();
    double phi00 = phi0();
    double l = CLHEP::c_light * beta() * (time - t0()) * cDip;
    double ang = phi00 + l * omega();
    double cang = cos(ang);
    double sang = sin(ang);
    double sphi0 = sin(phi00);
    double cphi0 = cos(phi00);

    return l2g_(Vec3((sang - sphi0) / omega() - d0() * sphi0, -(cang - cphi0) / omega() + d0() * cphi0, z0() + l * tanDip()));
  }

  Mom4 IPHelix::momentum(double time) const
  {

    Vec3 dir = direction(time);
    double bgm = betaGamma()*mass_;
    return Mom4(bgm*dir.X(), bgm*dir.Y(), bgm*dir.Z(), mass_);
  }

  void IPHelix::rangeInTolerance(TRange &brange, BField const &bfield, double tol) const
  {
    // precompute some factors
    double fact = 0.5 * sqrt(1./omega() * tol * bnom().R()) / CLHEP::c_light;
    // Limit to this traj's range
    brange.high() = std::min(brange.high(), range().high());
    // compute the BField difference in the middle of the range
    Vec3 midpos = position(brange.mid());
    Vec3 bvec  = bfield.fieldVect(midpos);
    auto db = bvec - bnom();
    double dt = fact / sqrt(db.R());
    // truncate the range if necessary
    if (dt < brange.range())
      brange.high() = brange.low() + dt;
  }

  double IPHelix::angle(const double &f) const
  {
    return phi0() + arc(f);
  }

  Vec3 IPHelix::velocity(double time) const
  {
    Mom4 mom = momentum(time);
    return mom.Vect() * (CLHEP::c_light * fabs(Q() / ebar()));
  }

  Vec3 IPHelix::direction(double time,LocalBasis::LocDir mdir) const
  {
    double cosval = cosDip();
    double sinval = sinDip();
    double phival = phi(time);

    switch ( mdir ) {
      case LocalBasis::perpdir:
        return l2g_(Vec3(-sinval * cos(phival), -sinval * sin(phival), cosval));
      case LocalBasis::phidir:
        return l2g_(Vec3(-sin(phival), cos(phival), 0.0));
      case LocalBasis::momdir:
        return l2g_(Vec3(Q() / omega() * cos(phival),
                         Q() / omega() * sin(phival),
                         Q() / omega() * tanDip()).Unit());
      default:
        throw std::invalid_argument("Invalid direction");
    }
  }

  IPHelix::DVEC IPHelix::momDeriv(double time, LocalBasis::LocDir mdir) const
  {
    // compute some useful quantities
    double tanval = tanDip();
    double cosval = cosDip();
    double omval = omega();
    double l = translen(CLHEP::c_light * beta() * (time - t0()));
    double d0val = d0();
    DVEC pder;
    // cases
    switch ( mdir ) {
      case LocalBasis::perpdir:
        // polar bending: only momentum and position are unchanged
        pder[d0_] = tanval*(1-cos(omval*l))/omval;
        pder[phi0_] = -tanval * sin(omval * l) / (1 + omval * d0val);
        pder[omega_] = omval * tanval;
        pder[z0_] = - l - tanval * tanval * sin(omval * l) / (omval * (1 + omval * d0val));
        pder[tanDip_] = 1 / (cosval * cosval);
        pder[t0_] = pder[z0_] / vz() + pder[tanDip_] * (time - t0()) * cosval * cosval / tanval;
        break;
      case LocalBasis::phidir:
        // Azimuthal bending: R, Lambda, t0 are unchanged
        pder[d0_] = -sin(omval * l) / (omval * cosval);
        pder[phi0_] = cos(omval * l) / (cosval * (1 + omval * d0val));
        pder[omega_] = 0;
        pder[z0_] = -tanval / (omval * cosval) * (1 - cos(omval * l) / (1 + omval * d0val));
        pder[tanDip_] = 0;
        pder[t0_] = pder[z0_] / vz();
        break;
      case LocalBasis::momdir:
        // fractional momentum change: position and direction are unchanged
        pder[d0_] = -(1 - cos(omval * l)) / omval;
        pder[phi0_] = sin(omval * l) / (1 + omval * d0val);
        pder[omega_] = -omval;
        pder[z0_] = -tanval * (l - sin(omval * l) / (omval * (1 + omval * d0val)));
        pder[tanDip_] = 0;
        pder[t0_] = pder[z0_] / vz();
        break;
      default:
        throw std::invalid_argument("Invalid direction");
    }
    return pder;
  }

  std::ostream& operator <<(std::ostream& ost, IPHelix const& hhel) {
    ost << " IPHelix parameters: ";
    for(size_t ipar=0;ipar < IPHelix::npars_;ipar++){
      ost << IPHelix::paramName(static_cast<IPHelix::ParamIndex>(ipar) ) << " : " << hhel.paramVal(ipar);
      if(ipar < IPHelix::npars_-1) ost << " , ";
    }
    return ost;
  }

} // KinKal namespace
