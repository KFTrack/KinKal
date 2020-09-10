#include "KinKal/CentralHelix.hh"
#include "KinKal/BFieldUtils.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  vector<string> CentralHelix::paramTitles_ = {
    "Distance of closest approach d_{0}",
    "Angle in the xy plane at closest approach #phi_{0}",
    "xy plane curvature of the track #omega",
    "Distance from the closest approach to the origin z_{0}",
    "Tangent of the track dip angle in the #rho - z projection tan#lambda",
    "Time at Z=0 Plane"};
  vector<string> CentralHelix::paramNames_ = {
  "d_{0}","#phi_{0}","#omega","z_{0}","tan#lambda","t_{0}"};
  vector<string> CentralHelix::paramUnits_ = {
      "mm", "rad", "rad", "mm", "", "ns"};
  string CentralHelix::trajName_("CentralHelix");  
  std::vector<std::string> const& CentralHelix::paramNames() { return paramNames_; }
  std::vector<std::string> const& CentralHelix::paramUnits() { return paramUnits_; }
  std::vector<std::string> const& CentralHelix::paramTitles() { return paramTitles_; }
  std::string const& CentralHelix::paramName(ParamIndex index) { return paramNames_[static_cast<size_t>(index)];}
  std::string const& CentralHelix::paramUnit(ParamIndex index) { return paramUnits_[static_cast<size_t>(index)]; }
  std::string const& CentralHelix::paramTitle(ParamIndex index) { return paramTitles_[static_cast<size_t>(index)];}
  string const& CentralHelix::trajName() { return trajName_; }

  CentralHelix::CentralHelix( VEC4 const& pos0, MOM4 const& mom0, int charge, double bnom, TimeRange const& range) : CentralHelix(pos0,mom0,charge,VEC3(0.0,0.0,bnom),range) {}
  CentralHelix::CentralHelix(VEC4 const &pos0, MOM4 const &mom0, int charge, VEC3 const &bnom,
                 TimeRange const &trange) : trange_(trange), mass_(mom0.M()), charge_(charge), bnom_(bnom)
  {
    // Transform into the system where Z is along the Bfield.  This is a pure rotation about the origin
    VEC4 pos(pos0);
    MOM4 mom(mom0);
    g2l_ = Rotation3D(AxisAngle(VEC3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    if(fabs(g2l_(bnom_).Theta()) > 1.0e-6)throw invalid_argument("Rotation Error");
    pos = g2l_(pos);
    mom = g2l_(mom);
    // create inverse rotation; this moves back into the original coordinate system
    l2g_ = g2l_.Inverse();
    double momToRad = 1.0/(BFieldUtils::cbar()*charge_*bnom_.R());
    mbar_ = -mass_ * momToRad;

    double pt = sqrt(mom.perp2());
    double radius = fabs(pt*momToRad);

    double lambda = -mom.z()*momToRad;
    double amsign = copysign(1.0, mbar_);

    VEC3 center = VEC3(pos.x() + mom.y()*momToRad, pos.y() - mom.x()*momToRad, 0.0);
    double rcent = sqrt(center.perp2());
    double fcent = center.phi();
    double centerx = rcent*cos(fcent);
    double centery = rcent*sin(fcent);

    param(omega_) = amsign/radius;
    param(tanDip_) = amsign*lambda/radius;
    param(d0_) = amsign*(rcent - radius);
    param(phi0_) = atan2(-amsign * centerx, amsign * centery);

    VEC3 pos3 = VEC3(pos.x(), pos.y(), pos.z());
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
    VEC4 testpos(pos0);
    // std::cout << "Testpos " << testpos << std::endl;
    position(testpos);
    MOM4 testmom = momentum(testpos.T());
    auto dp = testpos.Vect() - pos0.Vect();
    auto dm = testmom.Vect() - mom0.Vect();
    if(dp.R() > 1.0e-5 || dm.R() > 1.0e-5)throw invalid_argument("Rotation Error");
  }

  double CentralHelix::deltaPhi(double &phi, double refphi) const
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

  CentralHelix::CentralHelix(CentralHelix const& other, VEC3 const& bnom, double trot) : CentralHelix(other) {
    mbar_ *= bnom_.R()/bnom.R();
    bnom_ = bnom;
    pars_.parameters() += other.dPardB(trot,bnom);
    g2l_ = Rotation3D(AxisAngle(VEC3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    l2g_ = g2l_.Inverse();
  }

  CentralHelix::CentralHelix(Parameters const &pdata, CentralHelix const& other) : CentralHelix(other) {
    pars_ = pdata;
  }

  CentralHelix::CentralHelix(ParticleState const& pstate, int charge, VEC3 const& bnom, TimeRange const& range) :
    CentralHelix(pstate.position4(),pstate.momentum4(),charge,bnom,range) 
  {}

  CentralHelix::CentralHelix(ParticleStateMeasurement const& pstate, int charge, VEC3 const& bnom, TimeRange const& range) :
  CentralHelix(pstate.stateVector(),charge,bnom,range) {
  // derive the parameter space covariance from the global state space covariance
    DPDS dpds = dPardState(pstate.stateVector().time());
    pars_.covariance() = ROOT::Math::Similarity(dpds,pstate.stateCovariance());
  }

  void CentralHelix::position(VEC4 &pos) const
 {
    VEC3 pos3 = position(pos.T());
    pos.SetPx(pos3.X());
    pos.SetPy(pos3.Y());
    pos.SetPz(pos3.Z());
  }

  VEC4 CentralHelix::pos4(double time) const
 {
    VEC3 pos3 = position(time);
    return VEC4( pos3.X(), pos3.Y(), pos3.Z(), time);
  }

  VEC3 CentralHelix::position(double time) const
  {
    double cDip = cosDip();
    double phi00 = phi0();
    double l = CLHEP::c_light * beta() * (time - t0()) * cDip;
    double ang = phi00 + l * omega();
    double cang = cos(ang);
    double sang = sin(ang);
    double sphi0 = sin(phi00);
    double cphi0 = cos(phi00);

    return l2g_(VEC3((sang - sphi0) / omega() - d0() * sphi0, -(cang - cphi0) / omega() + d0() * cphi0, z0() + l * tanDip()));
  }

  MOM4 CentralHelix::momentum(double time) const
  {

    VEC3 dir = direction(time);
    double bgm = betaGamma()*mass_;
    return MOM4(bgm*dir.X(), bgm*dir.Y(), bgm*dir.Z(), mass_);
  }

  double CentralHelix::angle(const double &f) const
  {
    return phi0() + arc(f);
  }

  VEC3 CentralHelix::velocity(double time) const
  {
    MOM4 mom = momentum(time);
    return mom.Vect() * (CLHEP::c_light * fabs(Q() / ebar()));
  }

  VEC3 CentralHelix::direction(double time,MomBasis::Direction mdir) const
  {
    double cosval = cosDip();
    double sinval = sinDip();
    double phival = phi(time);

    switch ( mdir ) {
      case MomBasis::perpdir_:
        return l2g_(VEC3(-sinval * cos(phival), -sinval * sin(phival), cosval));
      case MomBasis::phidir_:
        return l2g_(VEC3(-sin(phival), cos(phival), 0.0));
      case MomBasis::momdir_:
        return l2g_(VEC3(Q() / omega() * cos(phival),
                         Q() / omega() * sin(phival),
                         Q() / omega() * tanDip()).Unit());
      default:
        throw std::invalid_argument("Invalid direction");
    }
  }

  DVEC CentralHelix::momDeriv(double time, MomBasis::Direction mdir) const
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
      case MomBasis::perpdir_:
        // polar bending: only momentum and position are unchanged
        pder[d0_] = tanval*(1-cos(omval*l))/omval;
        pder[phi0_] = -tanval * sin(omval * l) / (1 + omval * d0val);
        pder[omega_] = omval * tanval;
        pder[z0_] = - l - tanval * tanval * sin(omval * l) / (omval * (1 + omval * d0val));
        pder[tanDip_] = 1 / (cosval * cosval);
        pder[t0_] = pder[z0_] / vz() + pder[tanDip_] * (time - t0()) * cosval * cosval / tanval;
        break;
      case MomBasis::phidir_:
        // Azimuthal bending: R, Lambda, t0 are unchanged
        pder[d0_] = -sin(omval * l) / (omval * cosval);
        pder[phi0_] = cos(omval * l) / (cosval * (1 + omval * d0val));
        pder[omega_] = 0;
        pder[z0_] = -tanval / (omval * cosval) * (1 - cos(omval * l) / (1 + omval * d0val));
        pder[tanDip_] = 0;
        pder[t0_] = pder[z0_] / vz();
        break;
      case MomBasis::momdir_:
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

  std::ostream& operator <<(std::ostream& ost, CentralHelix const& hhel) {
    ost << " CentralHelix parameters: ";
    for(size_t ipar=0;ipar < CentralHelix::npars_;ipar++){
      ost << CentralHelix::paramName(static_cast<CentralHelix::ParamIndex>(ipar) ) << " : " << hhel.paramVal(ipar);
      if(ipar < CentralHelix::npars_-1) ost << " , ";
    }
    return ost;
  }

} // KinKal namespace
