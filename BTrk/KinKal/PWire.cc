#include "BTrk/KinKal/PWire.hh"
#include "BTrk/KinKal/Constants.hh"
#include <iostream>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  vector<string> PWire::paramTitles_ = {
    "Transverse DOCA to Z Axis",
    "Azimuth of POCA"
    "Z at POCA",
    "Cos Theta",
    "Time at POCA"}; 
  vector<string> PWire::paramNames_ = {
  "D0","Phi0","Z0","CTheta","Time0"};
  std::vector<std::string> const& PWire::paramNames() { return paramNames_; }
  std::vector<std::string> const& PWire::paramTitles() { return paramTitles_; }
  std::string const& PWire::paramName(paramIndex index) { return paramNames_[static_cast<size_t>(index)];}
  std::string const& PWire::paramTitle(paramIndex index) { return paramTitles_[static_cast<size_t>(index)];}

  PWire::PWire(Vec3 const& p0, Vec3 const& svel, double tmeas) {
  // check direction: must be clost to perpendicular to z axis
    static const Vec3 zdir(0.0,0.0,1.0);
    auto sdir = svel.Unit();
    static const double tol =0.02; // within 2% of perpendicular
    double zsdot = zdir.Dot(sdir);
    if(fabs(zsdot) > tol)
      cout << "PWire direction not perpendicular!" << endl; // should throw FIXME!

    double stheta2 = (1.0 -zsdot*zsdot);
    pars_[cost_] = zsdot;
    // find the POCA with the z axis
    double psdot = p0.Dot(sdir);
    double slen = (p0.Z()*zsdot - psdot)/stheta2;
    auto poca = p0 + sdir*slen;
    pars_[d0_] = poca.Rho();
    pars_[phi0_] = atan2(poca.Y(),poca.X());
    pars_[z0_] = poca.Z();
    // check
    if(fabs(poca.Z()+(psdot*zsdot - p0.Z())/stheta2) > 1e-5)
      cout << "POCA calculation failed!" << endl; // should throw FIXME!
    // sign velocity according to the counter-clockwise direction convention
    vel_ = sqrt(svel.Mag2()); // FIXME!
    // move the time to POCA
    pars_[t0_] = tmeas + slen*vel_;

  }

  void PWire::position(Vec4& pos) const {
    Vec3 pos3;
    position(pos.T(),pos3);
    pos.SetXYZT(pos3.X(),pos3.Y(),pos3.Z(),pos.T());
  }

  void PWire::pos0( Vec3& pos) const {
    pos = Vec3(pars_[d0_]*cos(pars_[phi0_]),
	pars_[d0_]*sin(pars_[phi0_]),
	pars_[z0_]);
  }

  void PWire::position(double time, Vec3& pos) const {
    Vec3 p0; pos0(p0);
    Vec3 dir; direction(time,dir);
    pos = p0 + vel_*(time-pars_[t0_])*dir;
  }


  void PWire::velocity(double time, Vec3& vel) const {
    direction(time,vel);
    vel *= vel_;
  }

  void PWire::direction(double time, Vec3& dir) const {
    double sint = sinTheta();
    dir.SetXYZ(-sint*sin(pars_[phi0_]), sint*cos(pars_[phi0_]), pars_[cost_]);
  }

}
