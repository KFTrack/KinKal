#include "KinKal/Trajectory/CentralHelix.hh"
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
    "D0","Phi0","Omega","Z0","TanDip","Time0"};
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
    double radToMom = BFieldMap::cbar()*charge_*bnom_.R();
    double momToRad = 1.0/radToMom;
    mbar_ = -mass_ * momToRad;
    // caches
    double pt = sqrt(mom.perp2());
    double radius = fabs(pt*momToRad);
    double amsign = sign();
    param(omega_) = amsign/radius;
    param(tanDip_) = mom.Z()/pt;
    // vector pointing to the circle center from the measurement point; this is perp to the transverse momentum
    double phimom = atan2(mom.Y(),mom.X());
    double phirm = phimom + amsign*M_PI_2;
    auto relpos = radius*VEC3(cos(phirm),sin(phirm),0.0);
    // center of the circle
    auto lcent = pos.Vect() + relpos;
    double rcent = sqrt(lcent.perp2());
    // central helix undefined for small center radius
    if(rcent < 1.0) throw invalid_argument("Central helix undefined for center at origin");
    double phicent = atan2(lcent.Y(),lcent.X());
    param(phi0_) = phicent - amsign*M_PI_2;
    // force phi0 in the range [-pi,pi];
    double nrot = round(phi0()/(2*M_PI));
    param(phi0_) -= nrot*2*M_PI;
    param(d0_) = amsign*(rcent-radius);
    // preliminary z0: this doesn't have winding correction yet
    double dphi = phimom-phi0();
    double z0 = pos.Z() - tanDip()*dphi/omega();
    // Z change for 1 revolution; sign is irrelevant
    double deltaz = 2*M_PI*tanDip()/omega();
    // compute the winding; it should minimize |z0|
    double nwind = round(z0/deltaz);
    param(z0_) = z0 - nwind*deltaz;
    // t0, also correcting for winding
    param(t0_) = pos.T() -(dphi + 2*M_PI*nwind)/Omega();
    // test
//    auto testpos = position3(pos0.T());
//    auto testmom = momentum3(pos0.T());
//    auto dp = testpos - pos0.Vect();
//    auto dm = testmom - mom0.Vect();
//    if(dp.R() > 1.0e-5 || dm.R() > 1.0e-5)throw invalid_argument("Construction Test Failure");
//    // check
//    auto lmom = localMomentum(pos0.T());
//    auto tcent = center();
//    if(fabs(lcent.phi()-tcent.phi())>1e-5 || fabs(lcent.perp2()-tcent.perp2()) > 1e-5){
//      cout << "center " << lcent << " test center " << tcent << endl;
//    }
//    if(fabs(tan(phi0()) +1.0/tan(lcent.phi())) > 1e-5){
//      cout << "phi0 " << phi0() << " test phi0 " << -1.0/tan(lcent.phi()) << endl;
//    }
//    double d0t = sign()*sqrt(lcent.perp2())-sqrt(lmom.perp2())/Q();
//    if(fabs(d0t - d0()) > 1e-5){
//      cout  << " d0 " << d0() << " d0 test " << d0t << endl;
//    }
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

  CentralHelix::CentralHelix(Parameters const &pdata, double mass, int charge, double bnom, TimeRange const& range) : trange_(range),  pars_(pdata), mass_(mass), charge_(charge), bnom_(VEC3(0.0,0.0,bnom)){
    // compute kinematic cache
    double momToRad = 1.0/(BFieldMap::cbar()*charge_*bnom);
    mbar_ = -mass_ * momToRad;
  }

  CentralHelix::CentralHelix(Parameters const &pdata, CentralHelix const& other) : CentralHelix(other) {
    pars_ = pdata;
  }

  CentralHelix::CentralHelix(ParticleState const& pstate, VEC3 const& bnom, TimeRange const& range) :
    CentralHelix(pstate.position4(),pstate.momentum4(),pstate.charge(),bnom,range)
  {}

  CentralHelix::CentralHelix(ParticleStateEstimate const& pstate, VEC3 const& bnom, TimeRange const& range) :
    CentralHelix((ParticleState)pstate,bnom,range) {
      // derive the parameter space covariance from the global state space covariance
      PSMAT dpds = dPardState(pstate.time());
      pars_.covariance() = ROOT::Math::Similarity(dpds,pstate.stateCovariance());
    }

  double CentralHelix::momentumVariance(double time) const {
    DVEC dMomdP(0.0,  0.0, -1.0/omega() , 0.0 , sinDip()*cosDip() , 0.0);
    dMomdP *= momentum(time);
    return ROOT::Math::Similarity(dMomdP,params().covariance());
  }

  double CentralHelix::positionVariance(double time, MomBasis::Direction mdir) const {
    auto dxdpvec = dXdPar(time);
    auto momdir = direction(time);
    auto posdir = MomBasis::direction(mdir, momdir);
    SVEC3 sdir(posdir.X(),posdir.Y(),posdir.Z());
    DVEC dxdp = sdir*dxdpvec;
    return ROOT::Math::Similarity(dxdp,params().covariance());
  }

  PMAT CentralHelix::planeCovariance(double time,Plane const& plane) const {
    // project covariance onto the U, V direction of the given plane
    // particle direction cannot be orthogonal to the plane normal
    auto momdir = direction(time);
    if(fabs(plane.normal().Dot(momdir)) < 1.0e-10)throw invalid_argument("Momentum direction lies in the plane");
    auto dxdpvec = dXdPar(time);
    SVEC3 uvec(plane.uDirection().X(),plane.uDirection().Y(),plane.uDirection().Z());
    SVEC3 vvec(plane.vDirection().X(),plane.vDirection().Y(),plane.vDirection().Z());
    PPMAT dPlanedPar;
    dPlanedPar.Place_in_row(uvec*dxdpvec,0,0);
    dPlanedPar.Place_in_row(vvec*dxdpvec,1,0);
    return ROOT::Math::Similarity(dPlanedPar,params().covariance());
  }

  VEC4 CentralHelix::position4(double time) const {
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

  VEC3 CentralHelix::acceleration(double time) const {
    auto phival = phi(time);
    auto locacc = acceleration()*VEC3(-sin(phival),cos(phival),0.0);
    return l2g_(locacc);
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
    double invrc = 1.0/sqrt(center().perp2());
    double omval = Omega();
    double inve = 1.0/energy();
    double dt = time-t0();
    double invc = 1.0/CLHEP::c_light;

    SVEC3 domega_dM (-omega()*fx, -omega()*fy,0.0);
    SVEC3 dtanDip_dM (-tanDip()*fx, -tanDip()*fy,1.0/pt);
    SVEC3 dphi0_dM = sign()*invrc*invqval*SVEC3(-sphi0,cphi0,0.0);
    SVEC3 dd0_dM = invqval*(sign()*invrc*SVEC3( center().Y(), -center().X(), 0.0) -
        (1.0/pt)*SVEC3(locmom.X(),locmom.Y(), 0.0));
    SVEC3 dt0_dM = -invqval*inve*invc*dphi(time)*SVEC3(locmom.X(), locmom.Y(), locmom.Z()) -
      (1.0/omval)*(SVEC3(-fy,fx,0) - dphi0_dM);
    SVEC3 dz0_dM = locmom.Z()*inve*inve*inve*CLHEP::c_light*dt*SVEC3(locmom.X(), locmom.Y(), locmom.Z()) +
      locmom.Z()*inve*CLHEP::c_light*dt0_dM -
      inve*CLHEP::c_light*dt*SVEC3(0.0,0.0,1.0);

    // euclidean space is column, parameter space is row
    DPDV dPdM;
    dPdM.Place_in_row(dd0_dM,d0_,0);
    dPdM.Place_in_row(dphi0_dM,phi0_,0);
    dPdM.Place_in_row(domega_dM,omega_,0);
    dPdM.Place_in_row(dtanDip_dM,tanDip_,0);
    dPdM.Place_in_row(dz0_dM,z0_,0);
    dPdM.Place_in_row(dt0_dM,t0_,0);
    return dPdM;
  }

  DPDV CentralHelix::dPardM(double time) const {
    // now rotate these into local space
    RMAT g2lmat;
    g2l_.GetRotationMatrix(g2lmat);
    return dPardMLoc(time)*g2lmat;
  }

  DVDP CentralHelix::dMdPar(double time) const {
    double cDip = cosDip();
    //    double sDip = tanDip()*cDip;
    double factor = Q()/omega();
    double dp = dphi(time);
    double ang = phi0()+dp;
    double cang = cos(ang);
    double sang = sin(ang);
    double bta = beta();
    auto lmom = localMomentum(time);
    SVEC3 momv(lmom.X(),lmom.Y(),lmom.Z());
    SVEC3 momperpv(-lmom.Y(), lmom.X(),0.0);

    SVEC3 dM_dd0 (0,0,0);
    SVEC3 dM_dphi0 (-factor*sang, factor*cang, 0);
    SVEC3 dM_domega = (1.0/omega())*(-momv + bta*bta*dp*momperpv);
    SVEC3 dM_dtanDip = (Q()/omega())*(SVEC3(0.0,0.0,1.0)- dp*bta*bta*cDip*cDip*tanDip()*SVEC3(-sang,cang,0.0));
    SVEC3 dM_dz0 (0,0,0);
    SVEC3 dM_dt0 = -bta*cDip*omega()*CLHEP::c_light*momperpv;
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
    double invrc = 1.0/sqrt(center().perp2());
    double cphi0 = cos(phi0());
    double sphi0 = sin(phi0());

    SVEC3 domega_dX (0,0,0);
    SVEC3 dtanDip_dX (0,0,0);
    SVEC3 dphi0_dX = (sign()*invrc)*SVEC3(-cphi0,-sphi0,0.0);
    // euclidean space is column, parameter space is row
    SVEC3 dd0_dX = (sign()*invrc)*SVEC3(center().X(), center().Y(), 0.0);
    SVEC3 dt0_dX = (1.0/Omega())*dphi0_dX;
    SVEC3 dz0_dX = (tanDip()/omega())*dphi0_dX + SVEC3(0.0,0.0,1.0);
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
    double dp = dphi(time);
    double phit = dp+phi0();
    double cDip = cosDip();
    double sphi = sin(phit);
    double cphi = cos(phit);
    double sphi0 = sin(phi0());
    double cphi0 = cos(phi0());
    double bta = beta();
    double invom = 1.0/omega();

    SVEC3 dX_dd0 (-sphi0, cphi0, 0);
    SVEC3 dX_dphi0 = invom*SVEC3( cphi - cphi0, sphi - sphi0, 0.0) - d0()*SVEC3( cphi0, sphi0, 0.0);
    SVEC3 dX_domega = invom*invom*(-SVEC3(sphi-sphi0,-cphi+cphi0,dp*tanDip()) + dp*bta*bta*SVEC3(cphi,sphi,tanDip()));
    SVEC3 dX_dz0 (0,0,1);
    SVEC3 dX_dtanDip = dp*invom*( -bta*bta*cDip*cDip*tanDip()*SVEC3(cphi,sphi,tanDip())  + SVEC3(0.0,0.0,1.0));
    SVEC3 dX_dt0 = -CLHEP::c_light*beta()*cDip*SVEC3(cphi, sphi, tanDip());
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

  PSMAT CentralHelix::dPardStateLoc(double time) const{
    // aggregate state from separate X and M derivatives; parameter space is row
    DPDV dPdX = dPardXLoc(time);
    DPDV dPdM = dPardMLoc(time);
    PSMAT dpds;
    dpds.Place_at(dPdX,0,0);
    dpds.Place_at(dPdM,0,3);
    return dpds;
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

  DVEC CentralHelix::dPardB(double time) const {
    auto lpos =localPosition(time);
    auto pdir = localDirection(time,MomBasis::phidir_);
    auto mperp = VEC3(pdir.Y(),-pdir.X(),0.0);
    auto cent = center();
    double rcval = rc();
    double rad = bendRadius();
    DVEC retval;
    retval[omega_] = omega();
    retval[tanDip_] = 0.0;
    retval[d0_] = sign()*( rad + cent.Dot(pdir)*rad/rcval);
    retval[phi0_] = -sign()*rad*cent.Dot(mperp)/(rcval*rcval);
    retval[z0_] = lpos.Z()-z0() + tanDip()*retval[phi0_]/omega();
    retval[t0_] = time-t0() + retval[phi0_]/Omega();
    return (1.0/bnom_.R())*retval;
  }

  DVEC CentralHelix::dPardB(double time, VEC3 const& BPrime) const {
    // rotate Bfield difference into local coordinate system
    VEC3 dB = g2l_(BPrime-bnom_);
    // find the parameter change due to BField magnitude change using component parallel to the local nominal Bfield (always along z)
    DVEC retval = dPardB(time)*dB.Z();
    // find the change in (local) position and momentum due to the rotation implied by the B direction change
    // work in local coordinate system to avoid additional matrix mulitplications
    auto xvec = localPosition(time);
    auto mvec = localMomentum(time);
    VEC3 BxdB =VEC3(0.0,0.0,1.0).Cross(dB)/bnomR();
    VEC3 dx = xvec.Cross(BxdB);
    VEC3 dm = mvec.Cross(BxdB);
    // convert these to a full state vector change
    ParticleState dstate(dx,dm,time,mass(),charge());
    // convert the change in (local) state due to rotation to parameter space
    retval += dPardStateLoc(time)*dstate.state();
    return retval;
  }

  VEC3 CentralHelix::center(double time) const {
    // local transverse position is at the center.  Use the time to define the Z position
    VEC3 cpos = center();
    cpos.SetZ(z0()+tanDip()*dphi(time)/omega());
    // transform to global coordinates
    auto gcpos = l2g_(cpos);
    return gcpos;
  }

  Ray CentralHelix::axis(double time) const {
    // direction is along Bnom, signed by pz
    VEC3 adir = bnom_.Unit();
    auto pzsign = sinDip();
    if(pzsign*adir.Z() < 0) adir.SetZ(-adir.Z());
    return Ray(adir,center(time));
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
