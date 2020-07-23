#include "KinKal/LHelix.hh"
#include "KinKal/BField.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  typedef ROOT::Math::SVector<double,3> SVec3;
  typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepStd<double,3,3> > RMAT; // algebraic rotation matrix
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
    // Transform into the system where Z is along the Bfield.  This is a pure rotation about the origin
    Vec4 pos(pos0);
    Mom4 mom(mom0);
    g2l_ = Rotation3D(AxisAngle(Vec3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
    if(fabs(g2l_(bnom_).Theta()) > 1.0e-6)throw invalid_argument("Rotation Error");
    pos = g2l_(pos);
    mom = g2l_(mom);
    // create inverse rotation; this moves back into the original coordinate system
    l2g_ = g2l_.Inverse();
    // compute some simple useful parameters
    double pt = mom.Pt(); 
    double phibar = mom.Phi();
    // translation factor from MeV/c to curvature radius in mm, B in Tesla; signed by the charge!!!
    double momToRad = 1.0/(BField::cbar()*charge_*bnom_.R());
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

  LHelix::LHelix( PDATA const& pdata, LHelix const& other) : LHelix(other) {
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
    PDATA::DVEC dMomdP(rad(), lam(),  0.0, 0.0 ,0.0 , 0.0);
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

  Vec3 LHelix::direction(double time, LocalBasis::LocDir mdir) const {
    double phival = phi(time);
    double invpb = sign()/pbar(); // need to sign
    switch ( mdir ) {
      case LocalBasis::perpdir:
	return l2g_(Vec3( lam()*cos(phival)*invpb,lam()*sin(phival)*invpb,-rad()*invpb));
      case LocalBasis::phidir:
	return l2g_(Vec3(-sin(phival),cos(phival),0.0));
      case LocalBasis::momdir:
	return l2g_(Vec3( rad()*cos(phival)*invpb,rad()*sin(phival)*invpb,lam()*invpb));
      default:
	throw invalid_argument("Invalid direction");
    }
  }


  // derivatives of momentum projected along the given basis WRT the 6 parameters, and the physical direction associated with that
  LHelix::DVEC LHelix::momDeriv(double time, LocalBasis::LocDir mdir) const {
    typedef ROOT::Math::SVector<double,3> SVec3;
    DPDV dPdM = dPardM(time);
    auto dir = direction(time,mdir);
    double mommag = momentumMag(time);
    return mommag*(dPdM*SVec3(dir.X(), dir.Y(), dir.Z()));
  }

  VMAT LHelix::momCovar(double time) const {
    return VMAT(); //TODO!
  }
  VMAT LHelix::posCovar(double time) const {
    return VMAT(); //TODO!
  }
  LHelix::DPDV LHelix::dPardX(double time) const {
    // euclidean space is column, parameter space is row
    double omval = omega();
    SVec3 zdir(0.0,0.0,1.0);
    SVec3 dCx_dX (1.0,0.0,0.0);
    SVec3 dCy_dX (0.0,1.0,0.0);
    SVec3 dphi0_dX = -zdir/lam();
    SVec3 dt0_dX = -zdir/(omval*lam());
    LHelix::DPDV dPdX;
    dPdX.Place_in_row(dCx_dX,cx_,0);
    dPdX.Place_in_row(dCy_dX,cy_,0);
    dPdX.Place_in_row(dphi0_dX,phi0_,0);
    dPdX.Place_in_row(dt0_dX,t0_,0);
// now rotate these into local space
    RMAT g2lmat;
    g2l_.GetRotationMatrix(g2lmat);
    return dPdX*g2lmat;
  }

  LHelix::DPDV LHelix::dPardM(double time) const {
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
    LHelix::DPDV dPdM;
    dPdM.Place_in_row(dR_dM,rad_,0);
    dPdM.Place_in_row(dL_dM,lam_,0);
    dPdM.Place_in_row(dCx_dM,cx_,0);
    dPdM.Place_in_row(dCy_dM,cy_,0);
    dPdM.Place_in_row(dphi0_dM,phi0_,0);
    dPdM.Place_in_row(dt0_dM,t0_,0);
    dPdM *= 1.0/Q();
// now rotate these into local space
    RMAT g2lmat;
    g2l_.GetRotationMatrix(g2lmat);
    return dPdM*g2lmat;
  }

  LHelix::DVDP LHelix::dXdPar(double time) const {
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
    LHelix::DVDP dXdP;
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

  LHelix::DVDP LHelix::dMdPar(double time) const {
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
    LHelix::DVDP dMdP;
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

  DSDP LHelix::dPardState(double time) const{
  // aggregate state from separate X and M derivatives; parameter space is row
    LHelix::DPDV dPdX = dPardX(time);
    LHelix::DPDV dPdM = dPardM(time);
    DPDS dpds;
    dpds.Place_at(dPdX,0,0);
    dpds.Place_at(dPdM,0,3);
    return dpds;
  }

  DPDS LHelix::dStatedPar(double time) const {
  // aggregate state from separate X and M derivatives; parameter space is column
    LHelix::DVDP dXdP = dXdPar(time);
    LHelix::DVDP dMdP = dMdPar(time);
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

  void LHelix::rangeInTolerance(TRange& drange, BField const& bfield, double tol) const {
    // compute scaling factor
    double bn = bnom_.R();
    double spd = speed(drange.low());
    double sfac = spd*spd/(bn*pbar());
    // estimate step size from initial BField difference
    Vec3 tpos = position(drange.low());
    Vec3 bvec = bfield.fieldVect(tpos);
    auto db = (bvec - bnom_).R();
    double tstep(0.1);
    // this next part should have the hard-coded numbers replaced by parameters.  Some calculations should move to BField FIXME!
    if(db > 1e-4) tstep = 0.2*sqrt(tol/(sfac*db)); // step increment from difference from nominal
    Vec3 dBdt = bfield.fieldDeriv(tpos,velocity(drange.low()));
    tstep = std::min(tstep, 0.5*cbrt(tol/(sfac*dBdt.R())));
    //
    // loop over the trajectory in fixed steps to compute integrals and domains.
    // step size is defined by momentum direction tolerance.
    drange.high() = drange.low();
    double dx(0.0);
    // advance till spatial distortion exceeds position tolerance or we reach the range limit
    do{
      // increment the range
      drange.high() += tstep;
      tpos = position(drange.high());
      bvec = bfield.fieldVect(tpos);
      // BField diff with nominal
      auto db = (bvec - bnom_).R();
      // spatial distortion accumulation
      dx += sfac*drange.range()*tstep*db;
    } while(fabs(dx) < tol && drange.high() < range().high());
    //    std::cout << "tstep " << tstep << " trange " << drange.range() << std::endl;
  }

  void LHelix::print(ostream& ost, int detail) const {
    auto perr = params().diagonal(); 
    ost << " LHelix " << range() << " parameters: ";
    for(size_t ipar=0;ipar < LHelix::npars_;ipar++){
      ost << LHelix::paramName(static_cast<LHelix::ParamIndex>(ipar) ) << " " << paramVal(ipar) << " +- " << perr(ipar);
      if(ipar < LHelix::npars_-1) ost << " ";
    }
    ost << " with rotation around Bnom " << bnom_ << endl;
  }

  ostream& operator <<(ostream& ost, LHelix const& lhel) {
    lhel.print(ost,0);
    return ost;
  }

} // KinKal namespace
