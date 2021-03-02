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
    // kinematic to geometric conversion
    double radToMom = BFieldUtils::cbar()*charge_*bnom_.R();
    double momToRad = 1.0/radToMom;
    mbar_ = -mass_ * momToRad;
    // caches
    double pt = sqrt(mom.perp2());
    vt_ = CLHEP::c_light * pt / mom.E();
    double radius = fabs(pt*momToRad);
    double amsign = copysign(1.0,charge_*bnom_.Z());
    param(omega_) = -amsign/radius;
    param(tanDip_) = mom.Z()/pt; 
// vector pointing from the center to the measurement point; this is perp to the transverse momentum
    double phimom = atan2(mom.Y(),mom.X());
    double phirm = phimom + amsign*M_PI_2;
    auto relpos = radius*VEC3(cos(phirm),sin(phirm),0.0);
// center of the circle
    lcent_ = pos.Vect()  - relpos;
    double rcent = sqrt(lcent_.perp2());
    // central helix undefined for small center radius
    if(rcent < 1.0) throw invalid_argument("Central helix undefined for center at origin");
    double phicent = atan2(lcent_.Y(),lcent_.X());
    param(phi0_) = phicent + amsign*M_PI_2;
    param(d0_) = -amsign*(rcent-radius);
    // preliminary z0: this doesn't have winding correction yet
    double z0 = pos.Z() + amsign*tanDip()*(phimom-phi0())/omega();
    // Z change for 1 revolution; sign is irrelevant
    double deltaz = 2*M_PI*tanDip()/omega();
    // compute the winding; it should minimize |z0|
    nwind_ = round(z0/deltaz);
    param(z0_) = z0 - nwind_*deltaz;
    // t0, also correcting for winding
    param(t0_) = pos.T() +amsign*(phimom-phi0() + 2*M_PI*nwind_)/Omega();
    // test
    auto testpos = position3(pos0.T());
    auto testmom = momentum3(pos0.T());
    auto dp = testpos - pos0.Vect();
    auto dm = testmom - mom0.Vect();
    if(dp.R() > 1.0e-5 || dm.R() > 1.0e-5)throw invalid_argument("Rotation Error");
    // check
    auto lpos = localPosition(pos0.T());
    auto lmom = localMomentum(pos0.T());
    VEC3 lcent(lpos.X()-lmom.Y()/Q(),lpos.Y()+lmom.X()/Q(),0.0);
    if(fabs(lcent_.phi()-lcent.phi())>1e-5){
      cout << "center " << lcent_ << " test center " << lcent << endl;
    }
    if(fabs(tan(phi0()) +1.0/tan(lcent.phi())) > 1e-5){
      cout << "phi0 " << phi0() << " test phi0 " << -1.0/tan(lcent.phi()) << endl; 
    }
    double d0t = sign()*lcent.R()-sqrt(lmom.perp2())/Q();
    if(fabs(d0t - d0()) > 1e-5){
      cout  << " d0 " << d0() << " d0 test " << d0t << endl;
    }
  }

  void CentralHelix::setBNom(double time, VEC3 const& bnom) {
    // adjust the parameters for the change in bnom
    mbar_ *= bnom_.R()/bnom.R();
    pars_.parameters() += dPardB(time,bnom);
    bnom_ = bnom;
    // adjust rotations to global space
    g2l_ = Rotation3D(AxisAngle(VEC3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    l2g_ = g2l_.Inverse();
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

  CentralHelix::CentralHelix(ParticleStateEstimate const& pstate, int charge, VEC3 const& bnom, TimeRange const& range) :
  CentralHelix(pstate.stateVector(),charge,bnom,range) {
  // derive the parameter space covariance from the global state space covariance
    PSMAT dpds = dPardState(pstate.stateVector().time());
    pars_.covariance() = ROOT::Math::Similarity(dpds,pstate.stateCovariance());
  }

  double CentralHelix::momentumVar(double time) const {
    DVEC dMomdP(0.0,  0.0, -1.0/omega() , 0.0 , sinDip()*cosDip() , 0.0);
    dMomdP *= momentum(time);
    return ROOT::Math::Similarity(dMomdP,params().covariance());
  }

  VEC4 CentralHelix::position4(double time) const
 {
    VEC3 pos3 = position3(time);
    return VEC4( pos3.X(), pos3.Y(), pos3.Z(), time);
  }

  VEC3 CentralHelix::position3(double time) const {
    return l2g_(localPosition(time));
  } 

  VEC3 CentralHelix::localPosition(double time) const
  {
    double phit = phi(time);
    double cphit = cos(phit);
    double sphit = sin(phit);
    double sphi0 = sin(phi0());
    double cphi0 = cos(phi0());
    double rho = 1.0/omega();

    return rho*VEC3((sphit - sphi0),  -(cphit - cphi0), tanDip()*(phit-phi0())) +
	VEC3(- d0()*sphi0, d0()*cphi0, z0());
  }

  VEC3 CentralHelix::momentum3(double time) const
  {
    return direction(time)*momentum();
  }

  MOM4 CentralHelix::momentum4(double time) const
  {
    VEC3 mom3 = momentum3(time);
    return MOM4(mom3.X(), mom3.Y(), mom3.Z(), mass_);
  }

  VEC3 CentralHelix::velocity(double time) const
  {
    return CLHEP::c_light * beta()*direction(time,MomBasis::momdir_);
  }

  VEC3 CentralHelix::localMomentum(double time) const{
    return betaGamma()*mass()*localDirection(time);
  }

  VEC3 CentralHelix::direction(double time, MomBasis::Direction mdir) const {
    return l2g_(localDirection(time,mdir));
  }

  VEC3 CentralHelix::localDirection(double time,MomBasis::Direction mdir) const
  {
    double cosdip = cosDip();
    double sindip = sinDip();
    double phit = phi(time);
    double cphit = cos(phit);
    double sphit = sin(phit);

    switch ( mdir ) {
      case MomBasis::perpdir_:
        return VEC3(-sindip * cphit, -sindip * sphit, cosdip);
      case MomBasis::phidir_:
        return VEC3(-sphit, cphit, 0.0);
      case MomBasis::momdir_:
        return VEC3(cosdip * cphit, cosdip* sphit, sindip);
      default:
        throw std::invalid_argument("Invalid direction");
    }
  }
  
  DPDV CentralHelix::dPardMLoc(double time) const {
    auto locmom = localMomentum(time);
    double pt2 = locmom.perp2();
    double pt = sqrt(pt2);
    double fx = locmom.X()/pt2;
    double fy = locmom.Y()/pt2;
    double invqval = 1.0/Q();
    double cphi0 = cos(phi0());
    double sphi0 = sin(phi0());
    double invrc = 1.0/sqrt(lcent_.perp2());
    double omval = Omega();
//    double deltaphi = dphi(time);
    double inve = 1.0/energy();
    double dt = time-t0();
    double invc = 1.0/CLHEP::c_light;

    SVEC3 domega_dM (-omega()*fx, -omega()*fy,0.0);
    SVEC3 dtanDip_dM (-tanDip()*fx, -tanDip()*fy,1.0/pt);
    SVEC3 dphi0_dM = sign()*invrc*invqval*SVEC3(-sphi0,cphi0,0.0);
    SVEC3 dd0_dM = invqval*(sign()*invrc*SVEC3( lcent_.Y(), -lcent_.X(), 0.0) -
     (1.0/pt)*SVEC3(locmom.X(),locmom.Y(), 0.0));
    SVEC3 dt0_dM = -sign()*invqval*inve*invc*dphi(time)*SVEC3(locmom.X(), locmom.Y(), locmom.Z()) -
      (sign()/omval)*(SVEC3(-fy,fx,0) - dphi0_dM);
    SVEC3 dz0_dM = locmom.Z()*inve*inve*inve*CLHEP::c_light*dt*SVEC3(locmom.X(), locmom.Y(), locmom.Z()) + 
      locmom.Z()*inve*CLHEP::c_light*dt0_dM -
      inve*CLHEP::c_light*dt*SVEC3(0.0,0.0,1.0);

    // euclidean space is column, parameter space is row
//    double phi00 = phi0();
//    double cDip = cosDip();
//    double l = CLHEP::c_light * beta() * (time - t0()) * cDip;
//    double rho = 1.0/omega();
//    double sphi0 = sin(phi00);
//    double cphi0 = cos(phi00);
//    double ang = phi00 + l * omega();
//    double cang = cos(ang);
//    double sang = sin(ang);
//    double x = (sang - sphi0) / omega() - d0() * sphi0;
//    double y = -(cang - cphi0) / omega() + d0() * cphi0;
//    double px = Q()/omega()*cang;
//    double py = Q()/omega()*sang;
//    double phip = atan2(py,px);
//    double cx = x-rho*sin(phip);
//    double cy = y+rho*cos(phip);
//
//    double omval = omega();
//    double tanval = tanDip();
//
//    // SVEC3 dz0_dM(0,0,0);
//    SVEC3 dt0_dM (CLHEP::c_light * beta() * cDip * Q() * Q() * sang,
//                  -CLHEP::c_light * beta() * cDip * Q() * Q() * cang,
//                  0);
//    SVEC3 domega_dM (-omval*omval*cang, -omval*omval*sang, 0);
//    SVEC3 dtanDip_dM (-omval*cang*tanval, -omval*sang*tanval, omval);
//    SVEC3 dphi0_dM (cx/(cx*cx+cy*cy), cy/(cx*cx+cy*cy), 0);
//
//    dphi0_dM /= Q();
//    domega_dM /= Q();
//    dtanDip_dM /= Q();
//
//    if( (cy*cphi0-cx*sphi0)*rho < 0.0 ){
//    // wrong angular momentum: fix
//      phi00 = phi00+CLHEP::pi;
//      sphi0 = -sphi0;
//      cphi0 = -cphi0;
//    }
//
//    double drho_dPx = omega()*locmom.X()/(Q()*Q());
//    double drho_dPy = omega()*locmom.Y()/(Q()*Q());
//
//    SVEC3 dd0_dM(1./sphi0 * (cx*cphi0/sphi0 * dphi0_dM[0]) - drho_dPx,
//                 1./sphi0 * (cx*cphi0/sphi0 * dphi0_dM[1] + 1./Q()) - drho_dPy,
//                 0);
//    if(fabs(sphi0)<=0.5)
//      SVEC3 dd0_dM(1./cphi0 * (cy*sphi0/cphi0 * dphi0_dM[0] + 1./Q()) - drho_dPx,
//                   1./cphi0 * (cy*sphi0/cphi0 * dphi0_dM[1]) - drho_dPy,
//                   0);

    // SVEC3 dz0_dM(- dtanDip_dM[0]*(phip-phi00)/omval
    //              + tanDip()*domega_dM[0]*(phip-phi00)/(omval*omval)
    //              - tanDip()*(-omval/Q()*sang-dphi0_dM[0])/omval,
    //              - dtanDip_dM[1]*(phip-phi00)/omval
    //              + tanDip()*domega_dM[1]*(phip-phi00)/(omval*omval)
    //              - tanDip()*(omval/Q()*cang-dphi0_dM[1])/omval,
    //              dtanDip_dM[2]*(phip-phi00)/omval);

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
    double invrc = 1.0/sqrt(lcent_.perp2());
    double cphi0 = cos(phi0());
    double sphi0 = sin(phi0());

    SVEC3 domega_dX (0,0,0);
    SVEC3 dtanDip_dX (0,0,0);
    SVEC3 dphi0_dX = (sign()*invrc)*SVEC3(-cphi0,sphi0,0.0);
    // euclidean space is column, parameter space is row
    SVEC3 dd0_dX = (sign()*invrc)*SVEC3(lcent_.X(), lcent_.Y(), 0.0);
    SVEC3 dt0_dX = (-sign()/Omega())*dphi0_dX;
    SVEC3 dz0_dX = SVEC3(0.0,0.0,1.0) + (tanDip()*Omega()/omega())*dt0_dX;
    DPDV dPdX;
    dPdX.Place_in_row(dd0_dX,d0_,0);
    dPdX.Place_in_row(dphi0_dX,phi0_,0);
    dPdX.Place_in_row(domega_dX,omega_,0);
    dPdX.Place_in_row(dz0_dX,z0_,0);
    dPdX.Place_in_row(dtanDip_dX,tanDip_,0);
    dPdX.Place_in_row(dt0_dX,t0_,0);
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

  PSMAT CentralHelix::dPardState(double time) const{
  // aggregate state from separate X and M derivatives; parameter space is row
    DPDV dPdX = dPardX(time);
    DPDV dPdM = dPardM(time);
    PSMAT dpds;
    dpds.Place_at(dPdX,0,0);
    dpds.Place_at(dPdM,0,3);
    return dpds;
  }

  PSMAT CentralHelix::dStatedPar(double time) const {
  // aggregate state from separate X and M derivatives; parameter space is column
    DVDP dXdP = dXdPar(time);
    DVDP dMdP = dMdPar(time);
    PSMAT dsdp;
    dsdp.Place_at(dXdP,0,0);
    dsdp.Place_at(dMdP,3,0);
    return dsdp;
  }

  ParticleStateEstimate CentralHelix::stateEstimate(double time) const {
  // express the parameter space covariance in global state space
    PSMAT dsdp = dStatedPar(time);
    return ParticleStateEstimate(state(time),ROOT::Math::Similarity(dsdp,pars_.covariance()));
  }

  void CentralHelix::print(std::ostream& ost, int detail) const { 
    ost << " CentralHelix parameters: ";
    for(size_t ipar=0;ipar < CentralHelix::npars_;ipar++){
      ost << CentralHelix::paramName(static_cast<CentralHelix::ParamIndex>(ipar) ) << " : " << paramVal(ipar);
      if(ipar < CentralHelix::npars_-1) ost << " , ";
    }
    ost << std::endl;
  }

  std::ostream& operator <<(std::ostream& ost, CentralHelix const& hhel) {
    hhel.print(ost,0);
    return ost;
  }

} // KinKal namespace
