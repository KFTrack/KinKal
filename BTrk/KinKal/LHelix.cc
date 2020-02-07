#include "BTrk/KinKal/LHelix.hh"
#include <math.h>

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
  std::vector<std::string> const& LHelix::paramNames() { return paramNames_; }
  std::vector<std::string> const& LHelix::paramTitles() { return paramTitles_; }
  std::string const& LHelix::paramName(paramIndex index) { return paramNames_[static_cast<size_t>(index)];}
  std::string const& LHelix::paramTitle(paramIndex index) { return paramTitles_[static_cast<size_t>(index)];}

  LHelix::LHelix( Vec4 const& pos, Mom4 const& mom, int charge, Context const& context) : KTraj(mom.M(),charge) {
    static double twopi = M_PI*M_PI;
    // compute some simple useful parameters
    double pt = mom.Pt(); 
    double phibar = mom.Phi();
    // translation factor from MeV/c to curvature radius in mm
    double momToRad = 1000.0/(charge_*context.Bz_*c_);
    // reduced mass; note sign convention!
    mbar_ = -mass_*momToRad;
    // transverse radius of the helix
    pars_.vec()[rad_] = -pt*momToRad;
    //tan dip
    pars_.vec()[lam_] = -mom.Z()*momToRad;
    // time at z=0
    double om = omega();
    pars_.vec()[t0_] = pos.T() - pos.Z()/(om*lam());
    // compute winding that miminimizes z1
    double nwind = rint((pos.Z()/lam() - phibar)/twopi);
    //  cout << "winding number = " << nwind << endl;
    // azimuth at z=0
    pars_.vec()[phi0_] = phibar - om*(pos.T()-t0()) + twopi*nwind;
    // circle center
    pars_.vec()[cx_] = pos.x() + mom.Y()*momToRad;
    pars_.vec()[cy_] = pos.X() - mom.X()*momToRad;
  }

  LHelix::LHelix( TDATA::DVec const& pvec, TDATA::DMat const& pcov, double mass, int charge, Context const& context) : KTraj(mass,charge), pars_(pvec,pcov) {
    double momToRad = 1000.0/(charge_*context.Bz_*c_);
    mbar_ = -mass_*momToRad;
  }

  void LHelix::position(Vec4& pos) const {
    // compute azimuthal angle
    double df = dphi(pos.T());
    double phival = df + phi0();
    // now compute position
    pos.SetPx(cx() + rad()*sin(phival));
    pos.SetPy(cy() - rad()*cos(phival));
    pos.SetPz(df*lam());
  }

  void LHelix::position(double t, Vec3& pos) const {
    // compute azimuthal angle
    double df = dphi(t);
    double phival = df + phi0();
    // now compute position
    pos.SetX(cx() + rad()*sin(phival));
    pos.SetY(cy() - rad()*cos(phival));
    pos.SetZ(df*lam());
 } 

  void LHelix::momentum(double tval,Mom4& mom) const{
    double phival = phi(tval);
    double factor = mass_/mbar_;
    mom.SetPx(factor * rad() * cos(phival));
    mom.SetPy(factor * rad() * sin(phival));
    mom.SetPz(factor * lam());
    mom.SetM(mass_);;
  }

  void LHelix::velocity(double tval,Vec3& vel) const{
    double phival = phi(tval);
    double factor = c_/ebar();
    vel.SetX(factor * rad() * cos(phival));
    vel.SetY(factor * rad() * sin(phival));
    vel.SetZ(factor * lam());
  }

  void LHelix::direction(double tval,Vec3& dir) const{
    double phival = phi(tval);
    double factor = 1.0/pbar();
    dir.SetX(factor * rad() * cos(phival));
    dir.SetY(factor * rad() * sin(phival));
    dir.SetZ(factor * lam());
  }

  void LHelix::dirVector(trajdir dir,double tval,Vec3& unit) const {
    double phival = phi(tval); // azimuth at this point
    double invpmm = 1.0/pbar(); 
    switch ( dir ) {
      case theta1:
	unit.SetX(lam()*cos(phival)*invpmm);
	unit.SetY(lam()*sin(phival)*invpmm);
	unit.SetZ(rad()*invpmm);
	break;
      case theta2:
	unit.SetX(sin(phival));
	unit.SetY(-cos(phival));
	unit.SetZ(0.0);
	break;
      case momdir:
	unit.SetX(rad()*cos(phival)*invpmm);
	unit.SetY(rad()*sin(phival)*invpmm);
	unit.SetZ(-lam()*invpmm);
	break;
      default:
      // should throw here FIXME!
	cout << "Error: unknown direction!" << dir << endl;
    }

  }

  void LHelix::momDeriv(trajdir dir, double time, PDer& dermat) const {
    // compute some useful quantities
    double bval = beta();
    double omval = omega();
    double pb = pbar();
    double phival = omval*(time - t0()) + phi0();
    // cases
    switch ( dir ) {
      case theta1:
	// polar bending: only momentum and position are unchanged
	dermat[rad_][0] = lam();
	dermat[lam_][0] = -rad();
	dermat[t0_][0] = (time-t0())/lam();
	dermat[phi0_][0] = omval*(time-t0())/lam();
	dermat[cx_][0] = -lam()*sin(phival);
	dermat[cy_][0] = lam()*cos(phival);
	break;
      case theta2:
	// Azimuthal bending: R, Lambda, t0 are unchanged
	dermat[rad_][0] = 0.0;
	dermat[lam_][0] = 0.0;
	dermat[t0_][0] = 0.0;
	dermat[phi0_][0] = copysign(1.0,omval)*pb/rad();
	dermat[cx_][0] = -copysign(1.0,omval)*pb*cos(phival);
	dermat[cy_][0] = -copysign(1.0,omval)*pb*sin(phival);
	break;
      case momdir:
	// fractional momentum change: position and direction are unchanged
	dermat[rad_][0] = rad();
	dermat[lam_][0] = lam();
	dermat[t0_][0] = (time-t0())*(1.0-bval*bval);
	dermat[phi0_][0] = omval*(time-t0());
	dermat[cx_][0] = -rad()*sin(phival);
	dermat[cy_][0] = +rad()*cos(phival);
	break;
      default:
	// should throw here FIXME!
	cout << "Error: unknown direction!" << dir << endl;
    }
  }
  std::ostream& operator <<(std::ostream& ost, LHelix const& lhel) {
    ost << " LHelix parameters: ";
    for(size_t ipar=0;ipar < LHelix::npars_;ipar++){
      ost << LHelix::paramName(static_cast<LHelix::paramIndex>(ipar) ) << " : " << lhel.param(ipar);
      if(ipar < LHelix::npars_-1) ost << " , ";
    }
    return ost;
  }

} // KinKal namespace
