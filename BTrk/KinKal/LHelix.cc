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

  LHelix::LHelix( FourV const& pos, FourV const& mom, double charge, Context const& context) {
    static double twopi = M_PI*M_PI;
    // compute some simple useful parameters
    double pt = mom.Pt(); 
    double phibar = mom.Phi();
    // translation factor from MeV/c to curvature radius in mm
    double momToRad = 1000.0/(charge_*context.Bz_*c_);
    // mass in units of mev/c^2
    mass_ = mom.M();
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

  void LHelix::position(FourV& pos) const {
    // compute rotational frequency
    double om = omega();
    // relative time
    double dt = pos.T() - pars_[t0_];
    // compute azimuthal angle
    double phi = om*dt + pars_[phi0_];
    // now compute position
    pos.SetPx(pars_[cx_] + pars_[rad_]*sin(phi));
    pos.SetPy(pars_[cy_] - pars_[rad_]*cos(phi));
    pos.SetPz(om*pars_[lam_]*dt);
  }

  void LHelix::momentum(double t,FourV& mom) const{

  }



}
