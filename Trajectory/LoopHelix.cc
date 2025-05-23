#include "KinKal/Trajectory/LoopHelix.hh"
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

  LoopHelix::LoopHelix() : mass_(0.0), charge_(0) {}
  LoopHelix::LoopHelix( VEC4 const& gpos, MOM4 const& gmom, int charge, double bnom, TimeRange const& range) : LoopHelix(gpos,gmom,charge,VEC3(0.0,0.0,bnom),range) {}
  LoopHelix::LoopHelix( VEC4 const& gpos, MOM4 const& gmom, int charge, VEC3 const& bnom, TimeRange const& trange) : trange_(trange), mass_(gmom.M()), charge_(charge), bnom_(bnom) {
    static double twopi = 2*M_PI;
    // Transform position and momentum into the system where Z is along the Bfield, which is the implicit coordinate system of the parameterization.
    // This is a pure rotation about the origin
    setTransforms();
    auto lpos = g2l_(gpos);
    auto lmom = g2l_(gmom);
    // compute some simple useful parameters
    double pt = lmom.Pt();
    double lmomphi = lmom.Phi();
    // translation factor from MeV/c to curvature radius in mm, B in Tesla; signed by the charge!!!
    double invq = 1.0/Q();
    // transverse radius of the helix
    param(rad_) = pt*invq;
    // longitudinal wavelength
    param(lam_) = lmom.Z()*invq;
    // time when particle reaches local z=0
    double zbar = lpos.Z()/lam();
    param(t0_) = lpos.T() - zbar/omega();
    // compute winding that puts phi0 in the range -pi,pi
    double nwind = round((zbar - lmomphi)/twopi);
    // particle momentum azimuth at z=0
    param(phi0_) = lmomphi - zbar + twopi*nwind;
    // circle center
    param(cx_) = lpos.X() - lmom.Y()*invq;
    param(cy_) = lpos.Y() + lmom.X()*invq;
  }

  void LoopHelix::syncPhi0(LoopHelix const& other) {
// adjust the phi0 of this traj to agree with the reference, keeping its value (mod 2pi) the same.
    static double twopi = 2*M_PI;
    int nloop = static_cast<int>(round( (other.phi0() - phi0())/twopi));
    if(nloop != 0) pars_.parameters()[phi0_] += nloop*twopi;
  }

  void LoopHelix::setTransforms() {
    // adjust rotations to global space
    g2l_ = Rotation3D(AxisAngle(VEC3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    l2g_ = g2l_.Inverse();
  }

  void LoopHelix::setBNom(double time, VEC3 const& newbnom) {
    auto db = newbnom - bnom_;
//    PSMAT dpdpdb =ROOT::Math::SMatrixIdentity();
//    PSMAT dpdpdb = dPardPardB(time,db);
//    std::cout << "dpdpdb = " << dpdpdb << std::endl;
//    pars_.covariance() = ROOT::Math::Similarity(dpdpdb,pars_.covariance());
    pars_.parameters() += dPardB(time,db);
    // rotate covariance: for now, just for the magnitude change.  Rotation is still TODO
    resetBNom(newbnom);
  }

  void LoopHelix::resetBNom(VEC3 const& newbnom) {
    bnom_ = newbnom;
    setTransforms();
  }

  LoopHelix::LoopHelix(LoopHelix const& other, VEC3 const& bnom, double tref) : LoopHelix(other) {
    setBNom(tref,bnom);
  }

  LoopHelix::LoopHelix( Parameters const& pdata, LoopHelix const& other) : LoopHelix(other) {
    pars_ = pdata;
  }

  LoopHelix::LoopHelix( Parameters const& pars, double mass, int charge, VEC3 const& bnom, TimeRange const& trange ) :
    trange_(trange), pars_(pars), mass_(mass), charge_(charge), bnom_(bnom) {
    setTransforms();
  }

  LoopHelix::LoopHelix(ParticleState const& pstate, VEC3 const& bnom, TimeRange const& range) :
    LoopHelix(pstate.position4(),pstate.momentum4(),pstate.charge(),bnom,range)
  {}

  LoopHelix::LoopHelix(ParticleStateEstimate const& pstate, VEC3 const& bnom, TimeRange const& range) :
    LoopHelix((ParticleState)pstate,bnom,range) {
      // derive the parameter space covariance from the global state space covariance
      PSMAT dpds = dPardState(pstate.time());
      pars_.covariance() = ROOT::Math::Similarity(dpds,pstate.stateCovariance());
    }

  double LoopHelix::momentumVariance(double time) const {
    DVEC dMomdP(rad(), lam(),  0.0, 0.0 ,0.0 , 0.0);
    dMomdP *= Q()/pbar();
    return ROOT::Math::Similarity(dMomdP,params().covariance());
  }

  double LoopHelix::positionVariance(double time, MomBasis::Direction mdir) const {
    auto dxdpvec = dXdPar(time);
    auto momdir = direction(time);
    auto posdir = MomBasis::direction(mdir, momdir);
    SVEC3 sdir(posdir.X(),posdir.Y(),posdir.Z());
    DVEC dxdp = sdir*dxdpvec;
    return ROOT::Math::Similarity(dxdp,params().covariance());
  }

  PMAT LoopHelix::planeCovariance(double time,Plane const& plane) const {
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

  VEC3 LoopHelix::acceleration(double time) const {
    auto phival = phi(time);
    auto locacc = acceleration()*VEC3(-sin(phival),cos(phival),0.0);
    return l2g_(locacc);
  }

  VEC3 LoopHelix::localDirection(double time, MomBasis::Direction mdir) const {
    double phival = phi(time);
    double invpb = -sign()/pbar(); // need to sign
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

  PSMAT LoopHelix::dPardPardB(double time,VEC3 const& db) const {
    auto newbnom = bnom_ + db;
    double bfrac = bnom_.R()/newbnom.R();
    auto xvec = localPosition(time);
    double phival = phi(time);
    double cphi = cos(phival);
    double sphi = sin(phival);
    // compute the rows
    DVEC dpdpdb_rad; dpdpdb_rad[rad_] = bfrac;
    DVEC dpdpdb_lam; dpdpdb_lam[lam_] = bfrac;
    DVEC dpdpdb_cx; dpdpdb_cx[cx_] = 1.0; dpdpdb_cx[rad_] = sphi*(1-bfrac);
    DVEC dpdpdb_cy; dpdpdb_cy[cy_] = 1.0; dpdpdb_cy[rad_] = -cphi*(1-bfrac);
    DVEC dpdpdb_phi0; dpdpdb_phi0[phi0_] = 1.0;  dpdpdb_phi0[lam_] = -xvec.Z()*(1-1.0/bfrac)/(lam()*lam());
    DVEC dpdpdb_t0; dpdpdb_t0[t0_] = 1.0;
    // build the matrix from these rows
    PSMAT dpdpdb;
    dpdpdb.Place_in_row(dpdpdb_rad,rad_,0);
    dpdpdb.Place_in_row(dpdpdb_lam,lam_,0);
    dpdpdb.Place_in_row(dpdpdb_cx,cx_,0);
    dpdpdb.Place_in_row(dpdpdb_cy,cy_,0);
    dpdpdb.Place_in_row(dpdpdb_phi0,phi0_,0);
    dpdpdb.Place_in_row(dpdpdb_t0,t0_,0);
    return dpdpdb;
  }

  DVEC LoopHelix::dPardB(double time, VEC3 const& db) const {
    // record position and momentum in local coordinates; these are constant
    auto xvec = localPosition(time);
    auto mvec = localMomentum(time);
    // translate newbnom to local coordinates
    VEC3 dbloc = g2l_(db);
    // find the parameter change due to bnom magnitude change.  These are exact
    double br = bnom_.R();
    auto zhat = VEC3(0.0,0.0,1.0);
    auto bloc = br*zhat;
    auto newbloc = dbloc + bloc;
    double nbr = newbloc.R();
    DVEC retval;
    double dbf  = (br - nbr)/nbr;
    retval[rad_] = rad()*dbf;
    retval[lam_] = lam()*dbf;
    retval[phi0_] = xvec.Z()*(br -nbr)/br/lam();
    retval[t0_] = 0.0;
    retval[cx_] = -sin(mvec.Phi())*retval[rad_];
    retval[cy_] = cos(mvec.Phi())*retval[rad_];

    // find the change in local position and momentum due to the rotation implied by the B direction change
    // work in local coordinate system to avoid additional matrix mulitplications
    VEC3 BxdB = zhat.Cross(dbloc)/br;
    VEC3 dx = xvec.Cross(BxdB);
    VEC3 dm = mvec.Cross(BxdB);
    // convert these to the state vector change
    ParticleState dstate(dx,dm,time,mass(),charge());
    // convert the (local) state due to rotation to parameter space
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

  VEC3 LoopHelix::center(double time) const {
    // local transverse position is at the center.  Use the time to define the Z position
    VEC3 cpos(cx(),cy(),dphi(time)*lam());
    // transform to global coordinates
    auto gcpos = l2g_(cpos);
    return gcpos;
  }

  Ray LoopHelix::axis(double time) const {
    // direction is along Bnom, signed by pz.  Note Bnom is in global coordinates
    VEC3 adir = bnom_.Unit();
    auto pzsign = -lam()*sign();
    if(pzsign*adir.Z() < 0) adir*= -1.0;
    return Ray(adir,center(time));
  }

  double LoopHelix::sagitta(double trange) const {
    double tlen = fabs(trange*transverseSpeed());
    double brad = bendRadius();
    if(tlen < M_PI*brad){
      double drunit = (1.0-cos(0.5*tlen/brad)); // unit circle
      return 0.125*brad*drunit*drunit;
    }
    return brad; // maximum possible sagitta
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
