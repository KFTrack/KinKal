#include "KinKal/IPHelix.hh"
#include "KinKal/BField.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  vector<string> IPHelix::paramTitles_ = {
    "Distance of closest approach",
    "Angle in the xy plane at closest approach",
    "xy plane curvature of the track",
    "Distance from the closest approach to the origin",
    "Tangent of the track dip angle in the rho-z projection",
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

  IPHelix::IPHelix(Vec4 const &pos, Mom4 const &mom, int charge, Vec3 const &bnom,
                 TRange const &trange) : trange_(trange), mass_(mom.M()), charge_(charge), bnom_(bnom)
  {
    double momToRad = 1000.0 / (charge_ * bnom_.R() * CLHEP::c_light);
    mbar_ = -mass_ * momToRad;

    double pt = sqrt(mom.perp2());
    double radius = fabs(pt*momToRad);

    double lambda = -mom.z()*momToRad;
    double amsign = copysign(1.0, -charge * bnom_.R());

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

  IPHelix::IPHelix(PDATA const &pdata, double mass, int charge, Vec3 const &bnom, TRange const &range) : IPHelix(pdata.parameters(),pdata.covariance(),mass,charge,bnom,range) {}


  IPHelix::IPHelix(PDATA::DVEC const &pvec, PDATA::DMAT const &pcov, double mass, int charge, Vec3 const &bnom,
                 TRange const &trange) :  trange_(trange), pars_(pvec, pcov), mass_(mass), charge_(charge), bnom_(bnom)
  {
    double momToRad = 1000.0 / (charge_ * bnom_.R() * CLHEP::c_light);
    mbar_ = -mass_ * momToRad;
  }

  void IPHelix::position(Vec4 &pos) const
  {
    double cDip = cosDip();
    double phi00 = phi0();
    double l = CLHEP::c_light * beta() * (pos.T() - t0()) * cDip;
    double ang = phi00 + l * omega();
    double cang = cos(ang);
    double sang = sin(ang);
    double sphi0 = sin(phi00);
    double cphi0 = cos(phi00);

    pos.SetPx((sang - sphi0) / omega() - d0() * sphi0);
    pos.SetPy(-(cang - cphi0) / omega() + d0() * cphi0);
    pos.SetPz(z0() + l * tanDip());
  }

  Vec3 IPHelix::position(double t) const
  {
    double cDip = cosDip();
    double phi00 = phi0();
    double l = CLHEP::c_light * beta() * (t - t0()) * cDip;
    double ang = phi00 + l * omega();
    double cang = cos(ang);
    double sang = sin(ang);
    double sphi0 = sin(phi00);
    double cphi0 = cos(phi00);

    return Vec3((sang - sphi0) / omega() - d0() * sphi0, -(cang - cphi0) / omega() + d0() * cphi0, z0() + l * tanDip());
  }

  Mom4 IPHelix::momentum(double tval) const
  {
    double l = beta() * CLHEP::c_light * (tval - t0()) * cosDip();
    return Mom4( Q() / omega() * cos(phi0() + omega() * l),
	Q() / omega() * sin(phi0() + omega() * l),
	Q() / omega() * tanDip(),
	mass_);
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
    double phival = phi(time);                   // azimuth at this point
    switch ( mdir ) {
      case LocalBasis::perpdir:
	return Vec3(-sinval * cos(phival), -sinval * sin(phival), cosval);
	break;
      case LocalBasis::phidir:
	return Vec3(-sin(phival), cos(phival), 0.0);
	break;
      case LocalBasis::momdir:
	return momentum(time).Vect().Unit();
	break;
      default:
	throw std::invalid_argument("Invalid direction");
    }
  }

  void IPHelix::momDeriv(double time, LocalBasis::LocDir mdir, DVEC &pder,Vec3& unit) const
  {
    // FIXME: these formulas need to be verified
    // compute some useful quantities
    double tanval = tanDip();
    double cosval = cosDip();
    double omval = omega();
    double l = translen(CLHEP::c_light * beta() * (time - t0()));
    double d0val = d0();
    unit = direction(time,mdir);
    // cases
    switch ( mdir ) {
      case LocalBasis::perpdir:
	// polar bending: only momentum and position are unchanged
	pder[d0_] = tanval*(1-cos(omval*l))/omval;
	pder[phi0_] = -tanval * sin(omval * l) / (1 + omval * d0val);
	pder[omega_] = omval * tanval;
	pder[z0_] = -l - tanval * tanval * sin(omval * l) / (omval * (1 + omval * d0val));
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
