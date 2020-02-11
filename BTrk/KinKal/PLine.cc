#include "BTrk/KinKal/PLine.hh"
#include "BTrk/KinKal/Constants.hh"
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  vector<string> PLine::paramTitles_ = {
    "Transverse DOCA to Z Axis",
    "Azimuth of POCA"
    "Z at POCA",
    "Cos Theta",
    "Time at POCA"}; 
  vector<string> PLine::paramNames_ = {
  "D0","Phi0","Z0","CTheta","Time0"};
  std::vector<std::string> const& PLine::paramNames() { return paramNames_; }
  std::vector<std::string> const& PLine::paramTitles() { return paramTitles_; }
  std::string const& PLine::paramName(paramIndex index) { return paramNames_[static_cast<size_t>(index)];}
  std::string const& PLine::paramTitle(paramIndex index) { return paramTitles_[static_cast<size_t>(index)];}

  PLine::PLine(Vec3 const& p0, Vec3 const& svel, double tmeas, TRange const& range) : TTraj(range) {
  // check direction: must be clost to perpendicular to z axis
    static const Vec3 zdir(0.0,0.0,1.0);
    static const double tol =0.02; // within 2% of perpendicular
    auto sdir = svel.Unit();
    double zsdot = zdir.Dot(sdir);
    if(fabs(zsdot) > tol)
      throw std::invalid_argument("PLine direction not perpendicular!");

    double stheta2 = (1.0 -zsdot*zsdot);
    pars_.vec()[cost_] = zsdot;
    // find the POCA with the z axis; this defines the reference point
    double psdot = p0.Dot(sdir);
    double slen = (p0.Z()*zsdot - psdot)/stheta2;
    auto poca = p0 + sdir*slen;
    pars_.vec()[d0_] = poca.Rho();
    pars_.vec()[phi0_] = atan2(poca.Y(),poca.X());
    pars_.vec()[z0_] = poca.Z();
    // check
    if(fabs(poca.Z()+(psdot*zsdot - p0.Z())/stheta2) > 1e-5)
      throw std::range_error("POCA calculation failed!");
    vel_ = sqrt(svel.Mag2());
    // move the time to POCA
    pars_.vec()[t0_] = tmeas - slen/vel_;
  }

  void PLine::position(Vec4& pos) const {
    Vec3 pos3;
    position(pos.T(),pos3);
    pos.SetXYZT(pos3.X(),pos3.Y(),pos3.Z(),pos.T());
  }

  void PLine::pos0( Vec3& pos) const {
    pos = Vec3(d0()*cos(phi0()), d0()*sin(phi0()), z0());
  }

  void PLine::position(double time, Vec3& pos) const {
    Vec3 p0; pos0(p0);
    Vec3 vel; velocity(time,vel);
    pos = p0 + (time-t0())*vel;
  }


  void PLine::velocity(double time, Vec3& vel) const {
    direction(time,vel);
    vel *= vel_;
  }

  void PLine::direction(double time, Vec3& dir) const {
    double sint = sinTheta();
    dir.SetXYZ(sint*sin(phi0()), -sint*cos(phi0()), cost());
  }

}
