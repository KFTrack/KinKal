#include "KinKal/LHelix.hh"
#include "KinKal/BField.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  vector<string> LHelix::paramTitles_ = {
    "Transverse Radius",
    "Longiduinal Wavelength",
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
  LHelix::LHelix( Vec4 const& pos0, Mom4 const& mom0, int charge, Vec3 const& bnom, TRange const& trange) : KInter(mom0.M(),charge), trange_(trange), bnom_(bnom), needsrot_(false) {
    static double twopi = 2*M_PI; // FIXME
    // Transform into the system where Z is along the Bfield.  This is a pure rotation about the origin
    Vec4 pos(pos0);
    Mom4 mom(mom0);
    if(fabs(bnom_.Theta()) >1.0e-6){
      needsrot_ = true;
      Rotation3D rot(AxisAngle(Vec3(sin(bnom_.Phi()),-cos(bnom_.Phi()),0.0),bnom_.Theta()));
      pos = rot(pos);
      mom = rot(mom);
      // create inverse rotation; this moves back into the original coordinate system
      brot_ = rot.Inverse();
      // check; in the internal coordinate system, B should be along the Z axis.
      auto test = rot(bnom_);
      if(fabs(test.Theta()) > 1.0e-6)throw invalid_argument("Rotation Error");
    }
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
    if(needsrot_){
     // test position and momentum function
      Vec4 testpos; testpos.SetE(pos0.T());
      Mom4 testmom;
      position(testpos);
      momentum(testpos.T(),testmom);
      auto dp = testpos.Vect() - pos0.Vect();
      auto dm = testmom.Vect() - mom0.Vect();
      if(dp.R() > 1.0e-3 || dm.R() > 1.0e-3)throw invalid_argument("Rotation Error");
    }
  }

  LHelix::LHelix( PDATA const& pdata, double mass, int charge, double bnom, TRange const& range) : LHelix(pdata,mass,charge,Vec3(0.0,0.0,bnom),range) {}
  LHelix::LHelix( PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& trange) : 
     KInter(mass,charge), trange_(trange), pars_(pdata), bnom_(bnom) {
      double momToRad = 1000.0/(charge_*bnom_.R()*CLHEP::c_light);
      // reduced mass; note sign convention!
      mbar_ = -mass_*momToRad;
    }
  
  double LHelix::momentumVar(float time) const {
    PDATA::DVEC dMomdP(rad(), lam(),  0.0, 0.0 ,0.0 , 0.0);
    dMomdP *= mass()/(pbar()*mbar());
    return ROOT::Math::Similarity(dMomdP,params().covariance());
  }

  void LHelix::position(Vec4& pos) const {
    Vec3 temp;
    position(pos.T(),temp);
    pos.SetXYZT(temp.X(),temp.Y(),temp.Z(),pos.T());
  }
  
  void LHelix::position(float time, Vec3& pos) const {
    // compute azimuthal angle
    double df = dphi(time);
    double phival = df + phi0();
    // now compute position
    pos.SetX(cx() + rad()*sin(phival));
    pos.SetY(cy() - rad()*cos(phival));
    pos.SetZ(df*lam());
    if(needsrot_) pos = brot_(pos);
 } 

  void LHelix::momentum(float time,Mom4& mom) const{
    Vec3 dir;
    direction(time,dir);
    double bgm = betaGamma()*mass_;
    mom.SetPx(bgm*dir.X());
    mom.SetPy(bgm*dir.Y());
    mom.SetPz(bgm*dir.Z());
    mom.SetM(mass_);
  }

 void LHelix::velocity(float time,Vec3& vel) const{
    direction(time,vel);
    vel *= speed(time); 
  }

  void LHelix::direction(float time,Vec3& dir) const{
    dir = direction(time);
  }

  Vec3 LHelix::direction(float time) const {
    return momDir(KInter::momdir,time);
  }
  
  Vec3 LHelix::rawMomDir(MDir mdir, float time) const {
    double phival = phi(time);
    double invpb = sign()/pbar(); // need to sign
    switch ( mdir ) {
      case theta1:
	return Vec3( lam()*cos(phival)*invpb,lam()*sin(phival)*invpb,-rad()*invpb);
      case theta2:
	return Vec3(-sin(phival),cos(phival),0.0);
      case momdir:
	return Vec3( rad()*cos(phival)*invpb,rad()*sin(phival)*invpb,lam()*invpb);
      default:
	throw invalid_argument("Invalid direction");
    }
  }
  
  Vec3 LHelix::momDir(MDir mdir, float time) const {
    if(needsrot_)
      return brot_(rawMomDir(mdir,time));
    else
      return rawMomDir(mdir,time);
  }

// derivatives of momentum projected along the given basis WRT the 6 parameters, and the physical direction associated with that
  void LHelix::momDeriv(MDir mdir, float time, DVEC& pder, Vec3& unit) const {
    // compute some useful quantities
    double bval = beta();
    double omval = omega();
    double pb = pbar()*sign(); // need to sign
    double dt = time-t0();
    double phival = omval*dt + phi0();
    // set unit vector.  remove this eventually FIXME!
    unit = momDir(mdir,time);
    // cases
    switch ( mdir ) {
      case theta1:
	// polar bending: only momentum and position are unchanged
	pder[rad_] = lam();
	pder[lam_] = -rad();
	pder[t0_] = -dt*rad()/lam();
	pder[phi0_] = -omval*dt*rad()/lam();
	pder[cx_] = -lam()*sin(phival);
	pder[cy_] = lam()*cos(phival);
	break;
      case theta2:
	// Azimuthal bending: R, Lambda, t0 are unchanged
	pder[rad_] = 0.0;
	pder[lam_] = 0.0;
	pder[t0_] = 0.0;
	pder[phi0_] = pb/rad();
	pder[cx_] = -pb*cos(phival);
	pder[cy_] = -pb*sin(phival);
	break;
      case momdir:
	// fractional momentum change: position and direction are unchanged
	pder[rad_] = rad();
	pder[lam_] = lam();
	pder[t0_] = dt*(1.0-bval*bval);
	pder[phi0_] = omval*dt;
	pder[cx_] = -rad()*sin(phival);
	pder[cy_] = +rad()*cos(phival);
	break;
      default:
	throw invalid_argument("Invalid direction");
    }
  }

  void LHelix::rangeInTolerance(TRange& drange, BField const& bfield, float dtol, float ptol) const {
    // compute scaling factor
    float bn = bnom_.R();
    float spd = speed(drange.low());
    float sfac = spd*spd/(bn*pbar());
    // loop over the trajectory in fixed steps to compute integrals and domains.
    // step size is defined by momentum direction tolerance.
    float tstep = dtol*ebar()/CLHEP::c_light;
    drange.high() = drange.low();
    float dx(0.0);
    // advance till spatial distortion exceeds position tolerance or we reach the range limit
    do{
      // increment the range
      drange.high() += tstep;
      Vec3 tpos, bvec;
      position(drange.high(),tpos);
      bfield.fieldVect(tpos,bvec);
      // BField diff with nominal
      auto dbvec = bvec - bnom_;
      // spatial distortion accumulation
      dx += sfac*drange.range()*tstep*dbvec.R();
    } while(fabs(dx) < ptol && drange.high() < range().high());
//     std::cout << "tstep " << tstep << " trange " << drange.range() << std::endl;
 }
 
  void LHelix::print(ostream& ost, int detail) const {
    auto perr = params().diagonal(); 
    ost << " LHelix " << range() << " parameters: ";
    for(size_t ipar=0;ipar < LHelix::npars_;ipar++){
      ost << LHelix::paramName(static_cast<LHelix::ParamIndex>(ipar) ) << " " << paramVal(ipar) << " +- " << perr(ipar);
      if(ipar < LHelix::npars_-1) ost << " ";
    }
    if(needsrot_) ost << " with rotation around Bnom " << bnom_;
    ost << endl;
  }

  ostream& operator <<(ostream& ost, LHelix const& lhel) {
    lhel.print(ost,0);
    return ost;
  }

} // KinKal namespace
