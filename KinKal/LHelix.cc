#include "KinKal/LHelix.hh"
#include "KinKal/BField.hh"
#include "KinKal/BFieldUtils.hh"
#include "Math/AxisAngle.h"
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
    vector<string> LHelix::paramTitles_ = {
    "Transverse Radius",
    "Longitudinal Wavelength",
    "Cylinder Center X",
    "Cylinder Center Y",
    "Azimuth at Z=0 Plane",
    "Time at Z=0 Plane"}; 
  vector<string> LHelix::paramNames_ = {
    "Radius","Lambda","CenterX","CenterY","Phi0","Time0"};
  vector<string> LHelix::paramUnits_ = {
    "mm","mm","mm","mm","radians","ns"};
  string LHelix::trajName_("LHelix");  
  vector<string> const& LHelix::paramNames() { return paramNames_; }
  vector<string> const& LHelix::paramUnits() { return paramUnits_; }
  vector<string> const& LHelix::paramTitles() { return paramTitles_; }
  string const& LHelix::paramName(ParamIndex index) { return paramNames_[static_cast<size_t>(index)];}
  string const& LHelix::paramUnit(ParamIndex index) { return paramUnits_[static_cast<size_t>(index)];}
  string const& LHelix::paramTitle(ParamIndex index) { return paramTitles_[static_cast<size_t>(index)];}
  string const& LHelix::trajName() { return trajName_; }

  LHelix::LHelix( Vec4 const& pos0, Mom4 const& mom0, int charge, double bnom, TRange const& range) : LHelix(pos0,mom0,charge,Vec3(0.0,0.0,bnom),range) {}
  LHelix::LHelix( Vec4 const& pos0, Mom4 const& mom0, int charge, Vec3 const& bnom, TRange const& trange) : trange_(trange), mass_(mom0.M()), charge_(charge), bnom_(bnom) {
    static double twopi = 2*M_PI;
    // Transform into the system where Z is along the Bfield, which is the implicit coordinate system of the parameterization.
    // The transform is a pure rotation about the origin
    Vec4 pos(pos0);
    Mom4 mom(mom0);
    g2l_ = Rotation3D(AxisAngle(Vec3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
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
    mbar_ = -mass_*momToRad;
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
    Vec4 testpos(pos0);
    position(testpos);
    Mom4 testmom = momentum(testpos.T());
    auto dp = testpos.Vect() - pos0.Vect();
    auto dm = testmom.Vect() - mom0.Vect();
    if(dp.R() > 1.0e-5 || dm.R() > 1.0e-5)throw invalid_argument("Rotation Error");
  }

  void LHelix::setBNom(double time, Vec3 const& bnom) {
    // adjust the parameters for the change in bnom
    pars_.parameters() += dPardB(time,bnom);
    bnom_ = bnom;
    // adjust rotations to global space
    g2l_ = Rotation3D(AxisAngle(Vec3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    l2g_ = g2l_.Inverse();
  }

  LHelix::LHelix(LHelix const& other, Vec3 const& bnom, double trot) : LHelix(other) {
    mbar_ *= bnom_.R()/bnom.R();
    bnom_ = bnom;
    pars_.parameters() += other.dPardB(trot,bnom);
    g2l_ = Rotation3D(AxisAngle(Vec3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    l2g_ = g2l_.Inverse();
  }

  LHelix::LHelix( PData const& pdata, LHelix const& other) : LHelix(other) {
    pars_ = pdata;
  }

  LHelix::LHelix(StateVector const& pstate, double time, double mass, int charge, Vec3 const& bnom, TRange const& range) :
    LHelix(Vec4(pstate.position().X(),pstate.position().Y(),pstate.position().Z(),time),
	Mom4(pstate.momentum().X(),pstate.momentum().Y(),pstate.momentum().Z(),mass),
	charge,bnom,range) 
  {}

  LHelix::LHelix(StateVectorMeasurement const& pstate, double time, double mass, int charge, Vec3 const& bnom, TRange const& range) :
  LHelix(pstate.stateVector(),time,mass,charge,bnom,range) {
  // derive the parameter space covariance from the global state space covariance
    DPDS dpds = dPardState(time);
    pars_.covariance() = ROOT::Math::Similarity(dpds,pstate.stateCovariance());
  }

  double LHelix::momentumVar(double time) const {
    DVEC dMomdP(rad(), lam(),  0.0, 0.0 ,0.0 , 0.0);
    dMomdP *= mass()/(pbar()*mbar());
    return ROOT::Math::Similarity(dMomdP,params().covariance());
  }

  Vec4 LHelix::pos4(double time) const {
    Vec3 temp = position(time);
    return Vec4(temp.X(),temp.Y(),temp.Z(),time);
  }

  void LHelix::position(Vec4& pos) const {
    Vec3 temp = position(pos.T());
    pos.SetXYZT(temp.X(),temp.Y(),temp.Z(),pos.T());
  }

  Vec3 LHelix::position(double time) const {
    double df = dphi(time);
    double phival = df + phi0();
    return l2g_(Vec3(cx() + rad()*sin(phival), cy() - rad()*cos(phival), df*lam()));
  } 

  Mom4 LHelix::momentum(double time) const{
    Vec3 dir = direction(time);
    double bgm = betaGamma()*mass_;
    return Mom4(bgm*dir.X(), bgm*dir.Y(), bgm*dir.Z(), mass_);
  }

  Vec3 LHelix::velocity(double time) const{
    return direction(time)*speed(time); 
  }

  Vec3 LHelix::localDirection(double time, MomBasis::Direction mdir) const {
    double phival = phi(time);
    double invpb = sign()/pbar(); // need to sign
    switch ( mdir ) {
      case MomBasis::perpdir_:
	return Vec3( lam()*cos(phival)*invpb,lam()*sin(phival)*invpb,-rad()*invpb);
      case MomBasis::phidir_:
	return Vec3(-sin(phival),cos(phival),0.0);
      case MomBasis::momdir_:
	return Vec3( rad()*cos(phival)*invpb,rad()*sin(phival)*invpb,lam()*invpb);
      default:
	throw invalid_argument("Invalid direction");
    }
  }

  Vec3 LHelix::localMomentum(double time) const{
    return betaGamma()*mass_*localDirection(time);
  }

  Vec3 LHelix::localPosition(double time) const {
    double df = dphi(time);
    double phival = df + phi0();
    return Vec3(cx() + rad()*sin(phival), cy() - rad()*cos(phival), df*lam());
  } 

  Vec3 LHelix::direction(double time, MomBasis::Direction mdir) const {
    return l2g_(localDirection(time,mdir));
  }

  // derivatives of parameters WRT momentum projected along the given momentum basis direction 
  DVEC LHelix::momDeriv(double time, MomBasis::Direction mdir) const {
    DPDV dPdM = dPardM(time);
    auto dir = direction(time,mdir);
    double mommag = momentumMag(time);
    return mommag*(dPdM*SVec3(dir.X(), dir.Y(), dir.Z())); // normalize to fractional change
  }

  DPDV LHelix::dPardXLoc(double time) const {
    // euclidean space is column, parameter space is row
    double omval = omega();
    SVec3 zdir(0.0,0.0,1.0);
    SVec3 dCx_dX (1.0,0.0,0.0);
    SVec3 dCy_dX (0.0,1.0,0.0);
    SVec3 dphi0_dX = -zdir/lam();
    SVec3 dt0_dX = -zdir/(omval*lam());
    DPDV dPdX;
    dPdX.Place_in_row(dCx_dX,cx_,0);
    dPdX.Place_in_row(dCy_dX,cy_,0);
    dPdX.Place_in_row(dphi0_dX,phi0_,0);
    dPdX.Place_in_row(dt0_dX,t0_,0);
    return dPdX;
  }

  DPDV LHelix::dPardMLoc(double time) const {
    // euclidean space is column, parameter space is row
    double omval = omega();
    double dt = time-t0();
    double dphi = omval*dt;
    double phival = dphi + phi0();
    double sphi = sin(phival);
    double cphi = cos(phival);
    double inve2 = 1.0/ebar2();
    SVec3 T2(-sphi,cphi,0.0);
    SVec3 T3(cphi,sphi,0.0);
    SVec3 zdir(0.0,0.0,1.0);
    SVec3 mdir = rad()*T3 + lam()*zdir; 
    SVec3 dR_dM = T3; 
    SVec3 dL_dM = zdir; 
    SVec3 dCx_dM (0.0,-1.0,0.0);
    SVec3 dCy_dM (1.0,0.0,0.0);
    SVec3 dphi0_dM = T2/rad() + (dphi/lam())*zdir;
    SVec3 dt0_dM = -dt*(inve2*mdir - zdir/lam());
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

  DPDV LHelix::dPardX(double time) const {
// rotate into local space
    RMAT g2lmat;
    g2l_.GetRotationMatrix(g2lmat);
    return dPardXLoc(time)*g2lmat;
  }

  DPDV LHelix::dPardM(double time) const {
// now rotate these into local space
    RMAT g2lmat;
    g2l_.GetRotationMatrix(g2lmat);
    return dPardMLoc(time)*g2lmat;
  }

  DVEC LHelix::dPardB(double time) const {
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

  DVEC LHelix::dPardB(double time, Vec3 const& BPrime) const {
  // rotate new B field difference into local coordinate system
    Vec3 dB = g2l_(BPrime-bnom_);
    // find the parameter change due to BField magnitude change usng component parallel to the local nominal Bfield (always along z)
    DVEC retval = dPardB(time)*dB.Z();
    // find the change in (local) position and momentum due to the rotation implied by the B direction change
    // work in local coordinate system to avoid additional matrix mulitplications
    auto xvec = localPosition(time);
    auto mvec = localMomentum(time);
    Vec3 BxdB =Vec3(0.0,0.0,1.0).Cross(dB)/bnomR();
    Vec3 dx = xvec.Cross(BxdB);
    Vec3 dm = mvec.Cross(BxdB);
    // convert these to a full state vector change
    StateVector dstate(dx,dm);
    // convert the change in (local) state due to rotation to parameter space
    retval += dPardStateLoc(time)*dstate.state();
    return retval;
  }

  DVDP LHelix::dXdPar(double time) const {
    // first find the derivatives wrt local cartesian coordinates
    // euclidean space is row, parameter space is column
    double omval = omega();
    double dt = time-t0();
    double dphi = omval*dt;
    double phival = dphi + phi0();
    double sphi = sin(phival);
    double cphi = cos(phival);
    double inve2 = 1.0/ebar2();
    SVec3 T2(-sphi,cphi,0.0);
    SVec3 T3(cphi,sphi,0.0);
    SVec3 zdir(0.0,0.0,1.0);
    SVec3 mdir = rad()*T3 + lam()*zdir;
    SVec3 dX_dR = -T2  -rad()*dphi*inve2*mdir;
    SVec3 dX_dL = dphi*zdir  -lam()*dphi*inve2*mdir;
    SVec3 dX_dCx (1.0,0.0,0.0); // along X
    SVec3 dX_dCy (0.0,1.0,0.0); // along Y
    SVec3 dX_dphi0 = rad()*T3;
    SVec3 dX_dt0 = -omval*mdir;
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

  DVDP LHelix::dMdPar(double time) const {
    double omval = omega();
    double dt = time-t0();
    double dphi = omval*dt;
    double phival = dphi + phi0();
    double sphi = sin(phival);
    double cphi = cos(phival);
    double inve2 = 1.0/ebar2();
    SVec3 T2(-sphi,cphi,0.0);
    SVec3 T3(cphi,sphi,0.0);
    SVec3 zdir(0.0,0.0,1.0);
    SVec3 dM_dphi0 = rad()*T2;
    SVec3 dM_dR = T3 -rad()*dphi*inve2*dM_dphi0;
    SVec3 dM_dL = zdir  -lam()*dphi*inve2*dM_dphi0;
    SVec3 dM_dt0 = -omval*dM_dphi0;
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

  DSDP LHelix::dPardStateLoc(double time) const{
  // aggregate state from separate X and M derivatives; parameter space is row
    DPDV dPdX = dPardXLoc(time);
    DPDV dPdM = dPardMLoc(time);
    DPDS dpds;
    dpds.Place_at(dPdX,0,0);
    dpds.Place_at(dPdM,0,3);
    return dpds;
  }

  DSDP LHelix::dPardState(double time) const{
  // aggregate state from separate X and M derivatives; parameter space is row
    DPDV dPdX = dPardX(time);
    DPDV dPdM = dPardM(time);
    DPDS dpds;
    dpds.Place_at(dPdX,0,0);
    dpds.Place_at(dPdM,0,3);
    return dpds;
  }

  DPDS LHelix::dStatedPar(double time) const {
  // aggregate state from separate X and M derivatives; parameter space is column
    DVDP dXdP = dXdPar(time);
    DVDP dMdP = dMdPar(time);
    DSDP dsdp;
    dsdp.Place_at(dXdP,0,0);
    dsdp.Place_at(dMdP,3,0);
    return dsdp;
  }

  StateVector LHelix::state(double time) const {
    return StateVector(position(time),momentum(time).Vect());
  }

  StateVectorMeasurement LHelix::measurementState(double time) const {
  // express the parameter space covariance in global state space
    DSDP dsdp = dStatedPar(time);
    return StateVectorMeasurement(state(time),ROOT::Math::Similarity(dsdp,pars_.covariance()));
  }

  void LHelix::print(ostream& ost, int detail) const {
    auto pvar = params().covariance().Diagonal();
    ost << " LHelix " << range() << " parameters: ";
    for(size_t ipar=0;ipar < LHelix::npars_;ipar++){
      ost << LHelix::paramName(static_cast<LHelix::ParamIndex>(ipar) ) << " " << paramVal(ipar) << " +- " << sqrt(pvar(ipar));
      if(ipar < LHelix::npars_-1) ost << " ";
    }
    ost << " with rotation around Bnom " << bnom_ << endl;
  }

  ostream& operator <<(ostream& ost, LHelix const& lhel) {
    lhel.print(ost,0);
    return ost;
  }

} // KinKal namespace
