#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Detector/BFieldUtils.hh"
#include "Math/AxisAngle.h"
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  const vector<string> LoopHelix::paramTitles_ = {
    "Transverse Radius",
    "Longitudinal Wavelength",
    "Cylinder Center X",
    "Cylinder Center Y",
    "Azimuth at Z=0 Plane",
    "Time at Z=0 Plane"}; 
  const vector<string> LoopHelix::paramNames_ = {
    "Radius","Lambda","CenterX","CenterY","Phi0","Time0"};
  const vector<string> LoopHelix::paramUnits_ = {
    "mm","mm","mm","mm","radians","ns"};
  const string LoopHelix::trajName_("LoopHelix");  
  vector<string> const& LoopHelix::paramNames() { return paramNames_; }
  vector<string> const& LoopHelix::paramUnits() { return paramUnits_; }
  vector<string> const& LoopHelix::paramTitles() { return paramTitles_; }
  string const& LoopHelix::paramName(ParamIndex index) { return paramNames_[static_cast<size_t>(index)];}
  string const& LoopHelix::paramUnit(ParamIndex index) { return paramUnits_[static_cast<size_t>(index)];}
  string const& LoopHelix::paramTitle(ParamIndex index) { return paramTitles_[static_cast<size_t>(index)];}
  string const& LoopHelix::trajName() { return trajName_; }

  LoopHelix::LoopHelix( VEC4 const& pos0, MOM4 const& mom0, int charge, double bnom, TimeRange const& range) : LoopHelix(pos0,mom0,charge,VEC3(0.0,0.0,bnom),range) {}
  LoopHelix::LoopHelix( VEC4 const& pos0, MOM4 const& mom0, int charge, VEC3 const& bnom, TimeRange const& trange) : trange_(trange), mass_(mom0.M()), charge_(charge), bnom_(bnom) {
    static double twopi = 2*M_PI;
    // Transform into the system where Z is along the Bfield, which is the implicit coordinate system of the parameterization.
    // The transform is a pure rotation about the origin
    VEC4 pos(pos0);
    MOM4 mom(mom0);
    g2l_ = Rotation3D(AxisAngle(VEC3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    if(fabs(g2l_(bnom_).Theta()) > 1.0e-6)throw invalid_argument("Rotation Error");
    // to convert global vectors into parameters they must first be rotated into the local system.
    pos = g2l_(pos);
    mom = g2l_(mom);
    // create inverse rotation; this moves back into the global coordinate system
    l2g_ = g2l_.Inverse();
    // compute some simple useful parameters
    double pt = mom.Pt(); 
    double phibar = mom.Phi();
    // translation factor from MeV/c to curvature radius in mm, B in Tesla; signed by the charge!!!
    double momToRad = 1.0/(BFieldUtils::cbar()*charge_*bnom_.R());
    // reduced mass; note sign convention!
    mbar_ = -mass()*momToRad;
    // transverse radius of the helix
    param(rad_) = -pt*momToRad;
    // longitudinal wavelength
    param(lam_) = -mom.Z()*momToRad;
    // time at z=0
    double om = omega();
    param(t0_) = pos.T() - pos.Z()/(om*lam());
    // compute winding that puts phi0 in the range -pi,pi
    double nwind = rint((pos.Z()/lam() - phibar)/twopi);
    //  cout << "winding number = " << nwind << endl;
    // azimuth at z=0
    param(phi0_) = phibar - om*(pos.T()-t0()) + twopi*nwind;
    // circle center
    param(cx_) = pos.X() + mom.Y()*momToRad;
    param(cy_) = pos.Y() - mom.X()*momToRad;
    // test position and momentum function
//    auto testpos = position3(pos0.T());
//    auto testmom = momentum3(pos0.T());
//    auto dp = testpos - pos0.Vect();
//    auto dm = testmom - mom0.Vect();
//    if(dp.R() > 1.0e-5 || dm.R() > 1.0e-5)throw invalid_argument("Rotation Error");
  }

  void LoopHelix::setBNom(double time, VEC3 const& bnom) {
    // adjust the parameters for the change in bnom
    mbar_ *= bnom_.R()/bnom.R();
    pars_.parameters() += dPardB(time,bnom);
    bnom_ = bnom;
    // adjust rotations to global space
    g2l_ = Rotation3D(AxisAngle(VEC3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    l2g_ = g2l_.Inverse();
  }

  LoopHelix::LoopHelix(LoopHelix const& other, VEC3 const& bnom, double tref) : LoopHelix(other) {
    setBNom(tref,bnom);
  }

  LoopHelix::LoopHelix( Parameters const& pdata, LoopHelix const& other) : LoopHelix(other) {
    pars_ = pdata;
  }

  LoopHelix::LoopHelix( Parameters const& pars, TimeRange const& trange, double mass, int charge, VEC3 const& bnom) : 
  trange_(trange), pars_(pars), mass_(mass), charge_(charge), bnom_(bnom) {
    double momToRad = 1.0/(BFieldUtils::cbar()*charge_*bnom_.R());
    // set reduced mass
    mbar_ = -mass_*momToRad;
    // set the transforms
    g2l_ = Rotation3D(AxisAngle(VEC3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    l2g_ = g2l_.Inverse();
  }

  LoopHelix::LoopHelix(ParticleState const& pstate, int charge, VEC3 const& bnom, TimeRange const& range) :
    LoopHelix(pstate.position4(),pstate.momentum4(),charge,bnom,range) 
  {}

  LoopHelix::LoopHelix(ParticleStateEstimate const& pstate, int charge, VEC3 const& bnom, TimeRange const& range) :
  LoopHelix(pstate.stateVector(),charge,bnom,range) {
  // derive the parameter space covariance from the global state space covariance
    PSMAT dpds = dPardState(pstate.stateVector().time());
    pars_.covariance() = ROOT::Math::Similarity(dpds,pstate.stateCovariance());
  }

  double LoopHelix::momentumVar(double time) const {
    DVEC dMomdP(rad(), lam(),  0.0, 0.0 ,0.0 , 0.0);
    dMomdP *= mass()/(pbar()*mbar());
    return ROOT::Math::Similarity(dMomdP,params().covariance());
  }

  VEC4 LoopHelix::position4(double time) const {
    VEC3 temp = position3(time);
    return VEC4(temp.X(),temp.Y(),temp.Z(),time);
  }

  VEC3 LoopHelix::position3(double time) const {
    return l2g_(localPosition(time));
  } 

  MOM4 LoopHelix::momentum4(double time) const{
    VEC3 mom3 = momentum3(time);
    return MOM4(mom3.X(), mom3.Y(), mom3.Z(), mass());
  }

  VEC3 LoopHelix::momentum3(double time) const{
    return direction(time)*momentum();
  }

  VEC3 LoopHelix::velocity(double time) const{
    return direction(time)*speed(time); 
  }

  VEC3 LoopHelix::localDirection(double time, MomBasis::Direction mdir) const {
    double phival = phi(time);
    double invpb = sign()/pbar(); // need to sign
    switch ( mdir ) {
      case MomBasis::perpdir_:
	return VEC3( lam()*cos(phival)*invpb,lam()*sin(phival)*invpb,-rad()*invpb);
      case MomBasis::phidir_:
	return VEC3(-sin(phival),cos(phival),0.0);
      case MomBasis::momdir_:
	return VEC3( rad()*cos(phival)*invpb,rad()*sin(phival)*invpb,lam()*invpb);
      default:
	throw invalid_argument("Invalid direction");
    }
  }

  VEC3 LoopHelix::localMomentum(double time) const{
    return localDirection(time)*momentum();
  }

  VEC3 LoopHelix::localPosition(double time) const {
    double df = dphi(time);
    double phival = df + phi0();
    return VEC3(cx() + rad()*sin(phival), cy() - rad()*cos(phival), df*lam());
  } 

  VEC3 LoopHelix::direction(double time, MomBasis::Direction mdir) const {
    return l2g_(localDirection(time,mdir));
  }

  // derivatives of parameters WRT momentum projected along the given momentum basis direction 
  DVEC LoopHelix::momDeriv(double time, MomBasis::Direction mdir) const {
    DPDV dPdM = dPardM(time);
    auto dir = direction(time,mdir);
    double mommag = momentum(time);
    return mommag*(dPdM*SVEC3(dir.X(), dir.Y(), dir.Z())); // normalize to fractional change
  }

  DPDV LoopHelix::dPardXLoc(double time) const {
    // euclidean space is column, parameter space is row
    double omval = omega();
    SVEC3 zdir(0.0,0.0,1.0);
    SVEC3 dCx_dX (1.0,0.0,0.0);
    SVEC3 dCy_dX (0.0,1.0,0.0);
    SVEC3 dphi0_dX = -zdir/lam();
    SVEC3 dt0_dX = -zdir/(omval*lam());
    DPDV dPdX;
    dPdX.Place_in_row(dCx_dX,cx_,0);
    dPdX.Place_in_row(dCy_dX,cy_,0);
    dPdX.Place_in_row(dphi0_dX,phi0_,0);
    dPdX.Place_in_row(dt0_dX,t0_,0);
    return dPdX;
  }

  DPDV LoopHelix::dPardMLoc(double time) const {
    // euclidean space is column, parameter space is row
    double omval = omega();
    double dt = time-t0();
    double dphi = omval*dt;
    double phival = dphi + phi0();
    double sphi = sin(phival);
    double cphi = cos(phival);
    double inve2 = 1.0/ebar2();
    SVEC3 T2(-sphi,cphi,0.0);
    SVEC3 dR_dM(cphi,sphi,0.0);
    SVEC3 dL_dM(0.0,0.0,1.0);
    SVEC3 mdir = rad()*dR_dM + lam()*dL_dM; 
    SVEC3 dCx_dM (0.0,-1.0,0.0);
    SVEC3 dCy_dM (1.0,0.0,0.0);
    SVEC3 dphi0_dM = T2/rad() + (dphi/lam())*dL_dM;
    SVEC3 dt0_dM = -dt*(inve2*mdir - dL_dM/lam());
    DPDV dPdM;
    dPdM.Place_in_row(dR_dM,rad_,0);
    dPdM.Place_in_row(dL_dM,lam_,0);
    dPdM.Place_in_row(dCx_dM,cx_,0);
    dPdM.Place_in_row(dCy_dM,cy_,0);
    dPdM.Place_in_row(dphi0_dM,phi0_,0);
    dPdM.Place_in_row(dt0_dM,t0_,0);
    dPdM *= 1.0/Q();
    return dPdM;
  }

  DPDV LoopHelix::dPardX(double time) const {
// rotate into local space
    RMAT g2lmat;
    g2l_.GetRotationMatrix(g2lmat);
    return dPardXLoc(time)*g2lmat;
  }

  DPDV LoopHelix::dPardM(double time) const {
// now rotate these into local space
    RMAT g2lmat;
    g2l_.GetRotationMatrix(g2lmat);
    return dPardMLoc(time)*g2lmat;
  }

  DVEC LoopHelix::dPardB(double time) const {
    double phival = phi(time);
    DVEC retval;
    retval[rad_] = -rad();
    retval[lam_] = -lam();
    retval[cx_] = rad()*sin(phival);
    retval[cy_] = -rad()*cos(phival);
    retval[phi0_] = -dphi(time); 
    retval[t0_] = 0.0;
    return (1.0/bnom_.R())*retval;
  }

  DVEC LoopHelix::dPardB(double time, VEC3 const& BPrime) const {
  // rotate new B field difference into local coordinate system
    VEC3 dB = g2l_(BPrime-bnom_);
    // find the parameter change due to BField magnitude change usng component parallel to the local nominal Bfield (always along z)
    DVEC retval = dPardB(time)*dB.Z();
    // find the change in (local) position and momentum due to the rotation implied by the B direction change
    // work in local coordinate system to avoid additional matrix mulitplications
    auto xvec = localPosition(time);
    auto mvec = localMomentum(time);
    VEC3 BxdB =VEC3(0.0,0.0,1.0).Cross(dB)/bnomR();
    VEC3 dx = xvec.Cross(BxdB);
    VEC3 dm = mvec.Cross(BxdB);
    // convert these to a full state vector change
    ParticleState dstate(dx,dm,time,mass());
    // convert the change in (local) state due to rotation to parameter space
    retval += dPardStateLoc(time)*dstate.state();
    return retval;
  }

  DVDP LoopHelix::dXdPar(double time) const {
    // first find the derivatives wrt local cartesian coordinates
    // euclidean space is row, parameter space is column
    double omval = omega();
    double dt = time-t0();
    double dphi = omval*dt;
    double phival = dphi + phi0();
    double sphi = sin(phival);
    double cphi = cos(phival);
    double inve2 = 1.0/ebar2();
    SVEC3 T2(-sphi,cphi,0.0);
    SVEC3 T3(cphi,sphi,0.0);
    SVEC3 zdir(0.0,0.0,1.0);
    SVEC3 mdir = rad()*T3 + lam()*zdir;
    SVEC3 dX_dR = -T2  -rad()*dphi*inve2*mdir;
    SVEC3 dX_dL = dphi*zdir  -lam()*dphi*inve2*mdir;
    SVEC3 dX_dCx (1.0,0.0,0.0); // along X
    SVEC3 dX_dCy (0.0,1.0,0.0); // along Y
    SVEC3 dX_dphi0 = rad()*T3;
    SVEC3 dX_dt0 = -omval*mdir;
    DVDP dXdP;
    dXdP.Place_in_col(dX_dR,0,rad_);
    dXdP.Place_in_col(dX_dL,0,lam_);
    dXdP.Place_in_col(dX_dCx,0,cx_);
    dXdP.Place_in_col(dX_dCy,0,cy_);
    dXdP.Place_in_col(dX_dphi0,0,phi0_);
    dXdP.Place_in_col(dX_dt0,0,t0_);
// now rotate these into global space
    RMAT l2gmat;
    l2g_.GetRotationMatrix(l2gmat);
    return l2gmat*dXdP;
  }

  DVDP LoopHelix::dMdPar(double time) const {
    double omval = omega();
    double dt = time-t0();
    double dphi = omval*dt;
    double phival = dphi + phi0();
    double sphi = sin(phival);
    double cphi = cos(phival);
    double inve2 = 1.0/ebar2();
    SVEC3 T2(-sphi,cphi,0.0);
    SVEC3 T3(cphi,sphi,0.0);
    SVEC3 zdir(0.0,0.0,1.0);
    SVEC3 dM_dphi0 = rad()*T2;
    SVEC3 dM_dR = T3 -rad()*dphi*inve2*dM_dphi0;
    SVEC3 dM_dL = zdir  -lam()*dphi*inve2*dM_dphi0;
    SVEC3 dM_dt0 = -omval*dM_dphi0;
    DVDP dMdP;
    dMdP.Place_in_col(dM_dR,0,rad_);
    dMdP.Place_in_col(dM_dL,0,lam_);
    dMdP.Place_in_col(dM_dphi0,0,phi0_);
    dMdP.Place_in_col(dM_dt0,0,t0_);
    dMdP *= Q(); // scale to momentum
// now rotate these into global space
    RMAT l2gmat;
    l2g_.GetRotationMatrix(l2gmat);
    return l2gmat*dMdP;
  }

  PSMAT LoopHelix::dPardStateLoc(double time) const{
  // aggregate state from separate X and M derivatives; parameter space is row
    DPDV dPdX = dPardXLoc(time);
    DPDV dPdM = dPardMLoc(time);
    PSMAT dpds;
    dpds.Place_at(dPdX,0,0);
    dpds.Place_at(dPdM,0,3);
    return dpds;
  }

  PSMAT LoopHelix::dPardState(double time) const{
  // aggregate state from separate X and M derivatives; parameter space is row
    DPDV dPdX = dPardX(time);
    DPDV dPdM = dPardM(time);
    PSMAT dpds;
    dpds.Place_at(dPdX,0,0);
    dpds.Place_at(dPdM,0,3);
    return dpds;
  }

  PSMAT LoopHelix::dStatedPar(double time) const {
  // aggregate state from separate X and M derivatives; parameter space is column
    DVDP dXdP = dXdPar(time);
    DVDP dMdP = dMdPar(time);
    PSMAT dsdp;
    dsdp.Place_at(dXdP,0,0);
    dsdp.Place_at(dMdP,3,0);
    return dsdp;
  }

  ParticleStateEstimate LoopHelix::stateEstimate(double time) const {
  // express the parameter space covariance in global state space
    PSMAT dsdp = dStatedPar(time);
    return ParticleStateEstimate(state(time),ROOT::Math::Similarity(dsdp,pars_.covariance()));
  }

  void LoopHelix::print(ostream& ost, int detail) const {
    auto pvar = params().covariance().Diagonal();
    ost << " LoopHelix " << range() << " parameters: ";
    for(size_t ipar=0;ipar < LoopHelix::npars_;ipar++){
      ost << LoopHelix::paramName(static_cast<LoopHelix::ParamIndex>(ipar) ) << " " << paramVal(ipar) << " +- " << sqrt(pvar(ipar));
      if(ipar < LoopHelix::npars_-1) ost << " ";
    }
    ost << " with rotation around Bnom " << bnom_ << endl;
  }

  ostream& operator <<(ostream& ost, LoopHelix const& lhel) {
    lhel.print(ost,0);
    return ost;
  }

} // KinKal namespace
