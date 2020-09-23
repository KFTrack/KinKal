#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Detector/BFieldUtils.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  const vector<string> CentralHelix::paramTitles_ = {
    "Distance of closest approach d_{0}",
    "Angle in the xy plane at closest approach #phi_{0}",
    "xy plane curvature of the track #omega",
    "Distance from the closest approach to the origin z_{0}",
    "Tangent of the track dip angle in the #rho - z projection tan#lambda",
    "Time at Z=0 Plane"};
  const vector<string> CentralHelix::paramNames_ = {
  "d_{0}","#phi_{0}","#omega","z_{0}","tan#lambda","t_{0}"};
  const vector<string> CentralHelix::paramUnits_ = {
      "mm", "rad", "rad", "mm", "", "ns"};
  const string CentralHelix::trajName_("CentralHelix");
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
    // std::cout << "Testpos " << testpos << std::endl;
    auto testpos = position3(pos0.T());
    auto testmom = momentum3(pos0.T());
    auto dp = testpos - pos0.Vect();
    auto dm = testmom - mom0.Vect();
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
    PSMAT dpds = dPardState(pstate.stateVector().time());
    pars_.covariance() = ROOT::Math::Similarity(dpds,pstate.stateCovariance());
  }

  VEC4 CentralHelix::position4(double time) const
 {
    VEC3 pos3 = position3(time);
    return VEC4( pos3.X(), pos3.Y(), pos3.Z(), time);
  }

  VEC3 CentralHelix::position3(double time) const
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

  VEC3 CentralHelix::momentum3(double time) const
  {
    return direction(time)* betaGamma()*mass_;
  }

  MOM4 CentralHelix::momentum4(double time) const
  {
    VEC3 mom3 = momentum3(time);
    return MOM4(mom3.X(), mom3.Y(), mom3.Z(), mass_);
  }

  double CentralHelix::angle(const double &f) const
  {
    return phi0() + arc(f);
  }

  VEC3 CentralHelix::velocity(double time) const
  {
    MOM4 mom = momentum4(time);
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

  DPDV CentralHelix::dPardMLoc(double time) const {
    // euclidean space is column, parameter space is row
    double phi00 = phi0();
    double cDip = cosDip();
    double l = CLHEP::c_light * beta() * (time - t0()) * cDip;
    double rho = 1.0/omega();
    double sphi0 = sin(phi00);
    double cphi0 = cos(phi00);
    double ang = phi00 + l * omega();
    double cang = cos(ang);
    double sang = sin(ang);
    double x = (sang - sphi0) / omega() - d0() * sphi0;
    double y = -(cang - cphi0) / omega() + d0() * cphi0;
    double px = Q()/omega()*cang;
    double py = Q()/omega()*sang;
    double phip = atan2(py,px);
    double cx = x-rho*sin(phip);
    double cy = y+rho*cos(phip);

    double omval = omega();
    double tanval = tanDip();

    // SVEC3 dz0_dM(0,0,0);
    SVEC3 dt0_dM (CLHEP::c_light * beta() * cDip * Q() * Q() * sang,
                  -CLHEP::c_light * beta() * cDip * Q() * Q() * cang,
                  0);
    SVEC3 domega_dM (-omval*omval*cang, -omval*omval*sang, 0);
    SVEC3 dtanDip_dM (-omval*cang*tanval, -omval*sang*tanval, omval);
    SVEC3 dphi0_dM (cx/(cx*cx+cy*cy), cy/(cx*cx+cy*cy), 0);

    dphi0_dM /= Q();
    domega_dM /= Q();
    dtanDip_dM /= Q();

    if( (cy*cphi0-cx*sphi0)*rho < 0.0 ){
    // wrong angular momentum: fix
      phi00 = phi00+CLHEP::pi;
      sphi0 = -sphi0;
      cphi0 = -cphi0;
    }

    double drho_dPx = omval*px/(Q()*Q());
    double drho_dPy = omval*py/(Q()*Q());

    SVEC3 dd0_dM(1./sphi0 * (cx*cphi0/sphi0 * dphi0_dM[0]) - drho_dPx,
                 1./sphi0 * (cx*cphi0/sphi0 * dphi0_dM[1] + 1./Q()) - drho_dPy,
                 0);
    if(fabs(sphi0)<=0.5)
      SVEC3 dd0_dM(1./cphi0 * (cy*sphi0/cphi0 * dphi0_dM[0] + 1./Q()) - drho_dPx,
                   1./cphi0 * (cy*sphi0/cphi0 * dphi0_dM[1]) - drho_dPy,
                   0);

    // SVEC3 dz0_dM(- dtanDip_dM[0]*(phip-phi00)/omval
    //              + tanDip()*domega_dM[0]*(phip-phi00)/(omval*omval)
    //              - tanDip()*(-omval/Q()*sang-dphi0_dM[0])/omval,
    //              - dtanDip_dM[1]*(phip-phi00)/omval
    //              + tanDip()*domega_dM[1]*(phip-phi00)/(omval*omval)
    //              - tanDip()*(omval/Q()*cang-dphi0_dM[1])/omval,
    //              dtanDip_dM[2]*(phip-phi00)/omval);

    SVEC3 dz0_dM(0,0,0);

    DPDV dPdM;
    dPdM.Place_in_row(dd0_dM,d0_,0);
    dPdM.Place_in_row(dphi0_dM,phi0_,0);
    dPdM.Place_in_row(domega_dM,omega_,0);
    dPdM.Place_in_row(dtanDip_dM,tanDip_,0);
    dPdM.Place_in_row(dz0_dM,z0_,0); // TODO
    dPdM.Place_in_row(dt0_dM,t0_,0); // TODO
    return dPdM;
  }

  DPDV CentralHelix::dPardM(double time) const {
    // now rotate these into local space
    RMAT g2lmat;
    g2l_.GetRotationMatrix(g2lmat);
    return dPardMLoc(time)*g2lmat;
  }

  DVDP CentralHelix::dMdPar(double time) const {
    double phi00 = phi0();
    double cDip = cosDip();
    double l = CLHEP::c_light * beta() * (time - t0()) * cDip;
    double omval = omega();
    double factor = Q()/omval;
    double ang = phi00 + l * omval;
    double cang = cos(ang);
    double sang = sin(ang);
    SVEC3 dM_dd0 (0,0,0);
    SVEC3 dM_dphi0 (-factor*sang, factor*cang, 0);
    SVEC3 dM_domega (-factor/omval*(l*omval*sang+cang),
                     factor/omval*(l*omval*cang-sang),
                     -factor/omval);
    SVEC3 dM_dtanDip (Q() * l * tanDip() * sang,
                      -Q() * l * tanDip() * cang,
                      factor);
    SVEC3 dM_dz0 (0,0,0);
    SVEC3 dM_dt0 (l/(time-t0()) * Q() * sang,
                  -l/(time-t0()) * Q() * cang,
                  0);
    DVDP dMdP;
    dMdP.Place_in_col(dM_dd0,0,d0_);
    dMdP.Place_in_col(dM_dphi0,0,phi0_);
    dMdP.Place_in_col(dM_domega,0,omega_);
    dMdP.Place_in_col(dM_dtanDip,0,tanDip_);
    dMdP.Place_in_col(dM_dz0,0,z0_);
    dMdP.Place_in_col(dM_dt0,0,t0_);
    // now rotate these into global space
    RMAT l2gmat;
    l2g_.GetRotationMatrix(l2gmat);
    return l2gmat*dMdP;
  }

  DPDV CentralHelix::dPardXLoc(double time) const {
    // euclidean space is column, parameter space is row
    SVEC3 dd0_dX (0,0,0);
    SVEC3 dz0_dX (0,0,0);
    SVEC3 dphi0_dX (0,0,0);
    SVEC3 dt0_dX (0,0,0);
    SVEC3 domega_dX (0,0,0);
    SVEC3 dtanDip_dX (0,0,0);
    DPDV dPdX;
    dPdX.Place_in_row(dd0_dX,d0_,0); // TODO
    dPdX.Place_in_row(dphi0_dX,phi0_,0); // TODO
    dPdX.Place_in_row(domega_dX,omega_,0); // TODO
    dPdX.Place_in_row(dz0_dX,z0_,0); // TODO
    dPdX.Place_in_row(dtanDip_dX,tanDip_,0); // TODO
    dPdX.Place_in_row(dt0_dX,t0_,0); // TODO
    return dPdX;
  }

  DVDP CentralHelix::dXdPar(double time) const {
    // first find the derivatives wrt local cartesian coordinates
    // euclidean space is row, parameter space is column
    double phi00 = phi0();
    double omval = omega();
    double d0val = d0();
    double cDip = cosDip();
    double sDip = sinDip();
    double l = CLHEP::c_light * beta() * (time - t0()) * cDip;
    double sang = sin(phi00+omval*l);
    double cang = cos(phi00+omval*l);
    SVEC3 dX_dd0 (-sin(phi00), cos(phi00), 0);
    SVEC3 dX_dphi0 (1./omval*cang - (1./omval+d0val)*cang,
                    1./omval*sang - (1./omval+d0val)*sang,
                    0);
    SVEC3 dX_domega ((l*omval*cang - sang + sin(phi00))/omval/omval,
                     (l*omval*sang + cang - cos(phi00))/omval/omval,
                     0);
    SVEC3 dX_dz0 (0,0,1);
    SVEC3 dX_dtanDip (-l*sDip*cDip*cang,
                      -l*sDip*cDip*sang,
                      l*cDip*cDip);
    SVEC3 dX_dt0 (-l/(time-t0())*cang,
                  -l/(time-t0())*sang,
                  -l/(time-t0())*tanDip());
    DVDP dXdP;
    dXdP.Place_in_col(dX_dd0,0,d0_);
    dXdP.Place_in_col(dX_dphi0,0,phi0_);
    dXdP.Place_in_col(dX_domega,0,omega_);
    dXdP.Place_in_col(dX_dz0,0,z0_);
    dXdP.Place_in_col(dX_dtanDip,0,tanDip_);
    dXdP.Place_in_col(dX_dt0,0,t0_);
    // now rotate these into global space
    RMAT l2gmat;
    l2g_.GetRotationMatrix(l2gmat);
    return l2gmat*dXdP;
  }

  DPDV CentralHelix::dPardX(double time) const {
    // rotate into local space
    RMAT g2lmat;
    g2l_.GetRotationMatrix(g2lmat);
    return dPardXLoc(time)*g2lmat;
  }

  DVEC CentralHelix::momDeriv(double time, MomBasis::Direction mdir) const
  {
    typedef ROOT::Math::SVector<double,3> SVEC3;
    DPDV dPdM = dPardM(time);
    auto dir = direction(time,mdir);
    double mommag = momentum(time);
    return mommag*(dPdM*SVEC3(dir.X(), dir.Y(), dir.Z()));
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
