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
    pars_[rad_] = -pt*momToRad;
    //tan dip
    pars_[lam_] = -mom.Z()*momToRad;
    // time at z=0
    double om = omega();
    pars_[t0_] = pos.T() - pos.Z()/(om*pars_[lam_]);
    // compute winding that miminimizes z1
    double nwind = rint((pos.Z()/(pars_[lam_]) - phibar)/twopi);
    //  cout << "winding number = " << nwind << endl;
    // azimuth at z=0
    pars_[phi0_] = phibar - om*(pos.T()-pars_[t0_]) + twopi*nwind;
    // circle center
    pars_[cx_] = pos.x() + mom.Y()*momToRad;
    pars_[cy_] = pos.X() - mom.X()*momToRad;
  }

  void LHelix::position(Vec4& pos) const {
    // compute azimuthal angle
    double df = dphi(pos.T());
    double phival = df + pars_[phi0_];
    // now compute position
    pos.SetPx(pars_[cx_] + pars_[rad_]*sin(phival));
    pos.SetPy(pars_[cy_] - pars_[rad_]*cos(phival));
    pos.SetPz(df*pars_[lam_]);
  }

  void LHelix::position(double t, Vec3& pos) const {
    // compute azimuthal angle
    double df = dphi(t);
    double phival = df + pars_[phi0_];
    // now compute position
    pos.SetX(pars_[cx_] + pars_[rad_]*sin(phival));
    pos.SetY(pars_[cy_] - pars_[rad_]*cos(phival));
    pos.SetZ(df*pars_[lam_]);
 } 

  void LHelix::momentum(double tval,Mom4& mom) const{
    double phival = phi(tval);
    double factor = mass_/mbar_;
    mom.SetPx(factor * pars_[rad_] * cos(phival));
    mom.SetPy(factor * pars_[rad_] * sin(phival));
    mom.SetPz(factor * pars_[lam_]);
    mom.SetM(mass_);;
  }

  void LHelix::velocity(double tval,Vec3& vel) const{
    double phival = phi(tval);
    double factor = c_/ebar();
    vel.SetX(factor * pars_[rad_] * cos(phival));
    vel.SetY(factor * pars_[rad_] * sin(phival));
    vel.SetZ(factor * pars_[lam_]);
  }

  void LHelix::direction(double tval,Vec3& dir) const{
    double phival = phi(tval);
    double factor = 1.0/pbar();
    dir.SetX(factor * pars_[rad_] * cos(phival));
    dir.SetY(factor * pars_[rad_] * sin(phival));
    dir.SetZ(factor * pars_[lam_]);
  }

  void LHelix::dirVector(trajdir dir,double tval,Vec3& unit) const {
    double phival = phi(tval); // azimuth at this point
    double invpmm = 1.0/pbar(); 
    switch ( dir ) {
      case theta1:
	unit.SetX(pars_[lam_]*cos(phival)*invpmm);
	unit.SetY(pars_[lam_]*sin(phival)*invpmm);
	unit.SetZ(pars_[rad_]*invpmm);
	break;
      case theta2:
	unit.SetX(sin(phival));
	unit.SetY(-cos(phival));
	unit.SetZ(0.0);
	break;
      case momdir:
	unit.SetX(pars_[rad_]*cos(phival)*invpmm);
	unit.SetY(pars_[rad_]*sin(phival)*invpmm);
	unit.SetZ(-pars_[lam_]*invpmm);
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
    double phival = omval*(time - pars_[t0_]) + pars_[phi0_];
    // cases
    switch ( dir ) {
      case theta1:
	// polar bending: only momentum and position are unchanged
	dermat[rad_][0] = pars_[lam_];
	dermat[lam_][0] = -pars_[rad_];
	dermat[t0_][0] = (time-pars_[t0_])/pars_[lam_];
	dermat[phi0_][0] = omval*(time-pars_[t0_])/pars_[lam_];
	dermat[cx_][0] = -pars_[lam_]*sin(phival);
	dermat[cy_][0] = pars_[lam_]*cos(phival);
	break;
      case theta2:
	// Azimuthal bending: R, Lambda, t0 are unchanged
	dermat[rad_][0] = 0.0;
	dermat[lam_][0] = 0.0;
	dermat[t0_][0] = 0.0;
	dermat[phi0_][0] = copysign(1.0,omval)*pb/pars_[rad_];
	dermat[cx_][0] = -copysign(1.0,omval)*pb*cos(phival);
	dermat[cy_][0] = -copysign(1.0,omval)*pb*sin(phival);
	break;
      case momdir:
	// fractional momentum change: position and direction are unchanged
	dermat[rad_][0] = pars_[rad_];
	dermat[lam_][0] = pars_[lam_];
	dermat[t0_][0] = (time-pars_[t0_])*(1.0-bval*bval);
	dermat[phi0_][0] = omval*(time-pars_[t0_]);
	dermat[cx_][0] = -pars_[rad_]*sin(phival);
	dermat[cy_][0] = +pars_[rad_]*cos(phival);
	break;
      default:
	// should throw here FIXME!
	cout << "Error: unknown direction!" << dir << endl;
    }
  }

} // KinKal namespace
