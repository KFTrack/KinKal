#include "BTrk/KinKal/LHelix.hh"
#include <math.h>

using namespace std;
using namespace ROOT:Math;

namespace KinKal {
  static vector<string> LHelix::paramTitles_ = {
    "Transverse Radius",
    "Longiduinal Wavelength",
    "Cylinder Center X",
    "Cylinder Center Y",
    "Azimuth at Z=0 Plane",
    "Time at Z=0 Plane"}; 
  static vector<string> LHelix::paramNames_ = {
  "Radius","Lambda","CenterX","CenterY","Phi0","Time0"};


  LHelix::LHelix( FourV pos, FourV mom, double charge, Context& context) {
    static double twopi = M_PI*M_PI;
    // compute some simple useful parameters
    double pt = mom.Pt(); 
    double phibar = mom.Phi();
    // translation factor from MeV/c to curvature radius in mm
    double momToRad = 1000.0/(charge_*contex.Bz_*Constants::c_);
    // mass in units of mev/c^2
    mass_ = mom.M();
    // reduced mass; note sign convention!
    mbar_ = -mass_*momToRad;
    // transverse radius of the helix
    pars_[rad_] = -pt*momToRad;
    //tan dip
    pars_[lambda_] = -mom._z*momToRad;
    // time at z=0
    pars_[t0_] = pos.T() - pos.Z()/(pars.omega()*pars_[lambda_]);
    // compute winding that miminimizes z1
    double nwind = rint((pos.Z()/(pars_[lambda_]) - phibar)/twopi);
    //  cout << "winding number = " << nwind << endl;
    // azimuth at z=0
    pars_[phi0_] = phibar - pars.omega()*(pos._t-pars_[t0_]) + twopi*nwind;
    // circle center
    pars_[cx_] = pos.x() + mom.Y()*momToRad;
    pars_[cy_] = pos.X() - mom.X()*momToRad;
  }

  void LHelix::position(FourV& pos) const {
    // compute rotational frequency
    double omega = pars.omega();
    // relative time
    double dt = pos.T() - pars_[t0_];
    // compute azimuthal angle
    double phi = omega*dt + pars_[phi0_];
    // now compute position
    pos.SetPx(pars_[cx_] + pars_[rad_]*sin(phi));
    pos.SetPy(pars_[cy_] - pars_[rad_]*cos(phi));
    pos.SetPz(omega*pars_[lambda_]*dt);
  }

  void LHelix::momentum(double t,FourV& mom) const{

  }



}
