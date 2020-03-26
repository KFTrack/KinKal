#include "KinKal/HHelix.hh"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  vector<string> HHelix::paramTitles_ = {
    "Distance of closest approach",
    "Angle in the xy plane at closest approach",
    "xy plane curvature of the track",
    "Distance from the closest approach to the origin",
    "Tangent of the track dip angle in the rho-z projection",
    "Time at Z=0 Plane"};
  vector<string> HHelix::paramNames_ = {
  "d0","phi0","omega","z0","tanDip","Time0"};
  std::vector<std::string> const& HHelix::paramNames() { return paramNames_; }
  std::vector<std::string> const& HHelix::paramTitles() { return paramTitles_; }
  std::string const& HHelix::paramName(paramIndex index) { return paramNames_[static_cast<size_t>(index)];}
  std::string const& HHelix::paramTitle(paramIndex index) { return paramTitles_[static_cast<size_t>(index)];}

  HHelix::HHelix(Vec4 const &pos, Mom4 const &mom, int charge, Context const &context,
                 TRange const &range) : TTraj(range), KInter(mom.M(), charge)
  {
    double momToRad = 1000.0 / (charge * context.bNom() * CLHEP::c_light);
    mbar_ = -mass_*momToRad;

    double pt = sqrt(mom.perp2());
    double radius = fabs(pt*momToRad);

    double lambda = -mom.z()*momToRad;
    float amsign = copysign(1.0, -charge * context.bNom());

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

  double HHelix::deltaPhi(double &phi, double refphi) const
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

  HHelix::HHelix(PDATA::DVEC const &pvec, PDATA::DMAT const &pcov, double mass, int charge, Context const &context,
                 TRange const &range) : TTraj(range), KInter(mass, charge), pars_(pvec, pcov)
  {
    double momToRad = 1000.0/(charge_*context.bNom()*CLHEP::c_light);
    mbar_ = -mass_*momToRad;
  }

  void HHelix::position(Vec4 &pos) const
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

  void HHelix::position(double t, Vec3 &pos) const
  {
    double cDip = cosDip();
    double phi00 = phi0();
    double l = CLHEP::c_light * beta() * (t - t0()) * cDip;
    double ang = phi00 + l * omega();
    double cang = cos(ang);
    double sang = sin(ang);
    double sphi0 = sin(phi00);
    double cphi0 = cos(phi00);

    pos.SetX((sang - sphi0) / omega() - d0() * sphi0);
    pos.SetY(-(cang - cphi0) / omega() + d0() * cphi0);
    pos.SetZ(z0() + l * tanDip());
  }

  void HHelix::momentum(double tval, Mom4 &mom) const
  {
    double l = beta() * CLHEP::c_light * (tval - t0()) * cosDip();
    mom.SetPx(Q() / omega() * cos(phi0() + omega() * l));
    mom.SetPy(Q() / omega() * sin(phi0() + omega() * l));
    mom.SetPz(Q() / omega() * tanDip());
    mom.SetM(mass_);
  }

  double HHelix::angle(const double &f) const
  {
    return phi0() + arc(f);
  }

  void HHelix::velocity(double tval, Vec3 &vel) const
  {
    Mom4 mom;
    momentum(tval, mom);
    vel = mom.Vect() * (CLHEP::c_light * fabs(Q() / ebar()));
  }

  void HHelix::direction(double tval, Vec3 &dir) const
  {
    Mom4 mom;
    momentum(tval, mom);
    dir = mom.Vect().Unit();
  }

  void HHelix::dirVector(MDir dir, double tval, Vec3 &unit) const
  {
    // FIXME: these formulas need to be verified
    double phival = phi(tval);                   // azimuth at this point
    // double norm = 1.0 / copysign(pbar(), mbar_); // sign matters!
    double sinval = sinDip();
    double cosval = cosDip();

    switch (dir)
    {
    case theta1:
      unit.SetX(-sinval * cos(phival));
      unit.SetY(-sinval * sin(phival));
      unit.SetZ(cosval);
      // unit *= norm;
      break;
    case theta2: // purely transverse
      unit.SetX(-sin(phival));
      unit.SetY(cos(phival));
      unit.SetZ(0.0);
      break;
    case momdir: // along momentum: sign matters!
      direction(tval, unit);
      break;
    default:
      throw std::invalid_argument("Invalid direction");
    }
  }

  void HHelix::momDeriv(MDir dir, double time, PDer& dermat) const {
    // FIXME: these formulas need to be verified
    // compute some useful quantities
    double tanval = tanDip();
    double cosval = cosDip();
    double omval = omega();
    double l = translen(CLHEP::c_light * beta() * (time - t0()));
    double d0val = d0();
    // cases
    switch ( dir ) {
      case theta1:
        // polar bending: only momentum and position are unchanged
        dermat[d0_][0] = tanval*(1-cos(omval*l))/omval;
        dermat[phi0_][0] = -tanval*sin(omval*l)/(1+omval*d0val);
        dermat[omega_][0] = omval*tanval;
        dermat[z0_][0] = -l-tanval*tanval*sin(omval*l)/(omval*(1+omval*d0val));
        dermat[tanDip_][0] = 1/(cosval*cosval);
        dermat[t0_][0] = dermat[z0_][0] / vz() + dermat[tanDip_][0] * (time - t0()) * cosval * cosval / tanval;
        break;
      case theta2:
        // Azimuthal bending: R, Lambda, t0 are unchanged
        dermat[d0_][0] = -sin(omval*l)/(omval*cosval);
        dermat[phi0_][0] = cos(omval*l)/(cosval*(1+omval*d0val));
        dermat[omega_][0] = 0;
        dermat[z0_][0] = -tanval/(omval*cosval)*(1-cos(omval*l)/(1+omval*d0val));
        dermat[tanDip_][0] = 0;
        dermat[t0_][0] = dermat[z0_][0]/vz();
        break;
      case momdir:
        // fractional momentum change: position and direction are unchanged
        dermat[d0_][0] = -(1-cos(omval*l))/omval;
        dermat[phi0_][0] = sin(omval*l)/(1+omval*d0val);
        dermat[omega_][0] = -omval;
        dermat[z0_][0] = -tanval*(l-sin(omval*l)/(omval*(1+omval*d0val)));
        dermat[tanDip_][0] = 0;
        dermat[t0_][0] = dermat[z0_][0] / vz();
        break;
      default:
        throw std::invalid_argument("Invalid direction");
    }
  }

  std::ostream& operator <<(std::ostream& ost, HHelix const& hhel) {
    ost << " HHelix parameters: ";
    for(size_t ipar=0;ipar < HHelix::npars_;ipar++){
      ost << HHelix::paramName(static_cast<HHelix::paramIndex>(ipar) ) << " : " << hhel.param(ipar);
      if(ipar < HHelix::npars_-1) ost << " , ";
    }
    return ost;
  }

} // KinKal namespace
